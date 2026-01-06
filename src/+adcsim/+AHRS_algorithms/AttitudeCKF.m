classdef AttitudeCKF < handle
%AttitudeCKF  Cubature Kalman Filter (CKF) for attitude estimation.
%
% Structure mirrors AttitudeUKF:
%   - Predict(gyro)  : high-rate propagation using gyro
%   - Correct(y_meas): low-rate correction using acc/mag/(stars)
%
% CKF specifics (Appendix A idea):
%   - Uses 2n cubature points (n = error-state dimension = 3)
%   - Points:  xi_i = ±sqrt(n) * e_i
%   - Weights: w_i = 1/(2n)
%
% State representation:
%   - Nominal attitude: quaternion (scalar-first) q = [w x y z]'
%   - Error state: 3x1 small rotation vector (delta-theta)
%   - Covariance stored for error state only: P_est (3x3)

    %% Public Properties
    properties (Access = public)
        % --- Core Filter Properties ---
        dt        double    % Sample period for the high-frequency loop (gyro).
        x_est     double    % Estimated attitude quaternion [w, x, y, z]'.
        P_est     double    % Estimated attitude ERROR covariance (3x3).
        gyro_bias double

        % --- Noise Covariances ---
        Q         double    % Process noise covariance for error-state (3x3).
        R         double    % Measurement noise covariance (6+3*Nstars x 6+3*Nstars).

        % --- Earth Reference Vectors ---
        g_ref     double
        m_ref     double
        stars_ref double

        imu
        starTracker
    end

    properties (Access = private)
        % --- CKF internal ---
        num_error_states  double   % n = 3
        num_cub_points    double   % 2n
        weight            double   % 1/(2n)
    end

    %% Public Methods
    methods (Access = public)
        function obj = AttitudeCKF(ckf_params, sample_time, imu_obj, st_obj)
            %AttitudeCKF Constructor.
            %
            % ckf_params fields (suggested, same as UKF noise params):
            %   .gyro_std (3x1)
            %   .acc_std  (3x1)
            %   .mag_std  (3x1)
            %   .star_std (3x1)  (per star)
            %
            obj.imu = imu_obj;

            % Star tracker handling
            if nargin > 3 && ~isempty(st_obj)
                obj.starTracker = st_obj;
            else
                obj.starTracker = struct('enable', false);
            end

            obj.dt    = sample_time;
            obj.x_est = [1; 0; 0; 0];
            obj.P_est = eye(3) * 0.1;
            obj.gyro_bias = zeros(3,1);

            % Error-state dimension
            obj.num_error_states = 3;          % attitude error only
            obj.num_cub_points   = 2*obj.num_error_states;
            obj.weight           = 1/obj.num_cub_points;

            % Process noise (gyro-driven error noise)
            obj.Q = diag(ckf_params.gyro_std.^2);

            % Measurement noise (acc/mag/stars stacked)
            R_acc = diag(ckf_params.acc_std.^2);
            R_mag = diag(ckf_params.mag_std.^2);

            if obj.starTracker.enable
                R_star_single = diag(ckf_params.star_std.^2);
                R_star = kron(eye(obj.starTracker.numberOfStars), R_star_single);
                obj.stars_ref = obj.starTracker.stars_I;
            else
                R_star = [];
                obj.stars_ref = [];
            end

            obj.R = blkdiag(R_acc, R_mag, R_star);

            % References in inertial frame
            obj.g_ref = obj.imu.g_I;
            obj.m_ref = obj.imu.m_I;
        end

        function Predict(obj, gyro)
            %PREDICT  CKF time update using gyro.
            %
            % Appendix A (concept):
            %   1) Build cubature points around (x_est, P_est)
            %   2) Propagate each point through process model
            %   3) Recover mean quaternion + covariance

            % 1) Cubature quaternion points
            q_points = obj.generate_cubature_points(obj.x_est, obj.P_est);

            % 2) Propagate each point through process model
            q_prop = zeros(4, obj.num_cub_points);
            for i = 1:obj.num_cub_points
                q_prop(:,i) = obj.process_model(q_points(:,i), gyro, obj.dt);
            end

            % 3) Recover predicted mean & covariance in error space
            [q_pred, P_pred] = obj.recover_statistics_quat_ckf(q_prop, obj.Q);

            obj.x_est = q_pred;
            obj.P_est = P_pred;
        end

        function Correct(obj, y_meas)
            %CORRECT  CKF measurement update using vector observations.
            %
            % y_meas stacking:
            %   [acc; mag; stars]  where stars is 3*Nstars stacked (if enabled)

            % Predicted state
            q_pred = obj.x_est;
            P_pred = obj.P_est;

            % Split / normalize measurements (recommended)
            acc = y_meas(1:3);
            mag = y_meas(4:6);

            if obj.starTracker.enable
                star = y_meas(7:end);
                numberOfStars = obj.starTracker.numberOfStars;
            else
                star = [];
                numberOfStars = 0;
            end

            % (Optional) normalize measured vectors to reduce scale issues
            acc = acc / max(norm(acc), eps);
            mag = mag / max(norm(mag), eps);
            if numberOfStars > 0
                for k = 1:numberOfStars
                    s = star(3*k-2:3*k);
                    star(3*k-2:3*k) = s / max(norm(s), eps);
                end
            end

            z_actual = [acc; mag; star];
            num_meas = 6 + 3*numberOfStars;

            % 1) Cubature points around predicted
            q_points = obj.generate_cubature_points(q_pred, P_pred);

            % 2) Push points through measurement model
            Z_points = zeros(num_meas, obj.num_cub_points);
            for i = 1:obj.num_cub_points
                Z_points(:,i) = obj.measurement_model(q_points(:,i));
            end

            % 3) Predicted measurement mean
            z_pred = obj.weight * sum(Z_points, 2);

            % 4) Innovation covariance Pzz
            Pzz = obj.R;
            for i = 1:obj.num_cub_points
                dz = Z_points(:,i) - z_pred;
                Pzz = Pzz + obj.weight * (dz*dz');
            end
            Pzz = (Pzz + Pzz')/2;

            % 5) Cross-covariance Pxz (in error-state space)
            Pxz = zeros(obj.num_error_states, num_meas);
            for i = 1:obj.num_cub_points
                % error vector from point quaternion relative to mean quaternion
                q_err = quatmultiply(q_points(:,i)', quatinv(q_pred'))';
                dx = obj.quat_to_error_vec(q_err);   % 3x1
                dz = Z_points(:,i) - z_pred;         % mx1
                Pxz = Pxz + obj.weight * (dx * dz');
            end

            % 6) Kalman gain & correction
            K = Pxz / Pzz;
            innovation = z_actual - z_pred;

            dtheta = K * innovation;                 % 3x1
            dq = obj.error_vec_to_quat(dtheta);      % 4x1

            obj.x_est = quatmultiply(dq', q_pred')';
            obj.x_est = obj.x_est / norm(obj.x_est);

            obj.P_est = P_pred - K * Pzz * K';
            obj.P_est = (obj.P_est + obj.P_est')/2;
        end
    end

    %% Private Methods
    methods (Access = private)
        function q_points = generate_cubature_points(obj, q_mean, P_error)
            % Generate 2n cubature quaternion points from mean quaternion and 3x3 error covariance.
            n = obj.num_error_states;

            % Cholesky (regularize if needed)
            try
                S = chol(P_error, 'lower');
            catch
                S = chol(P_error + 1e-9*eye(n), 'lower');
            end

            % Cubature directions: ±sqrt(n) e_i
            Xi = sqrt(n) * [eye(n), -eye(n)];  % 3 x 2n

            q_points = zeros(4, obj.num_cub_points);
            for i = 1:obj.num_cub_points
                dtheta = S * Xi(:,i);                % 3x1
                dq = obj.error_vec_to_quat(dtheta);  % 4x1
                q_points(:,i) = quatmultiply(dq', q_mean')';
                q_points(:,i) = q_points(:,i) / norm(q_points(:,i));
            end
        end

        function [mean_q, P_error] = recover_statistics_quat_ckf(obj, q_points, q_noise_cov)
            % Recover quaternion mean and 3x3 error covariance from 2n propagated quaternion points.
            %
            % Similar to UKF's quaternion averaging, but weights are uniform and no central point.

            q_avg = q_points(:,1);
            error_vectors = zeros(obj.num_error_states, obj.num_cub_points);

            % Iterative quaternion mean (small number of iterations is enough)
            for iter = 1:3
                for i = 1:obj.num_cub_points
                    q_err = quatmultiply(q_points(:,i)', quatinv(q_avg'))';
                    error_vectors(:,i) = obj.quat_to_error_vec(q_err);
                end
                mean_err = obj.weight * sum(error_vectors, 2);
                dq_mean = obj.error_vec_to_quat(mean_err);
                q_avg = quatmultiply(dq_mean', q_avg')';
                q_avg = q_avg / norm(q_avg);
            end

            mean_q = q_avg;

            % Covariance in error space
            P_error = q_noise_cov;
            for i = 1:obj.num_cub_points
                de = error_vectors(:,i) - mean_err;
                P_error = P_error + obj.weight * (de*de');
            end
            P_error = (P_error + P_error')/2;
        end

        function z_pred = measurement_model(obj, q)
            % Predicts stacked measurement vectors in BODY frame given quaternion q.
            R_i2b = quat2dcm(q');     % MATLAB Aerospace Toolbox (expects [w x y z])

            acc_pred = R_i2b * obj.g_ref;
            mag_pred = R_i2b * obj.m_ref;

            if obj.starTracker.enable
                Ns = obj.starTracker.numberOfStars;
                star_pred = zeros(3*Ns,1);
                for k = 1:Ns
                    star_pred(3*k-2:3*k) = R_i2b * obj.stars_ref(:,k);
                end
            else
                star_pred = [];
            end

            % Normalize predicted vectors (recommended)
            acc_pred = acc_pred / max(norm(acc_pred), eps);
            mag_pred = mag_pred / max(norm(mag_pred), eps);
            if ~isempty(star_pred)
                for k = 1:(numel(star_pred)/3)
                    s = star_pred(3*k-2:3*k);
                    star_pred(3*k-2:3*k) = s / max(norm(s), eps);
                end
            end

            z_pred = [acc_pred; mag_pred; star_pred];
        end
    end

    %% Static Helpers (same as your UKF)
    methods (Static, Access = private)
        function q_next = process_model(q, gyro, dt)
            omega = [0; gyro(:)];
            q_dot = 0.5 * quatmultiply(q', omega')';
            q_next = q + q_dot * dt;
            q_next = q_next / norm(q_next);
        end

        function err_vec = quat_to_error_vec(q_err)
            q_err = q_err / norm(q_err);
            angle = 2 * acos(q_err(1));
            if angle > pi, angle = angle - 2*pi; end

            if abs(angle) > 1e-9
                axis = q_err(2:4) / sin(angle/2);
                err_vec = angle * axis;
            else
                err_vec = zeros(3,1);
            end
        end

        function q_err = error_vec_to_quat(err_vec)
            angle = norm(err_vec);
            if angle > 1e-9
                axis = err_vec / angle;
                q_err = [cos(angle/2); axis * sin(angle/2)];
            else
                q_err = [1; 0; 0; 0];
            end
        end
    end
end
