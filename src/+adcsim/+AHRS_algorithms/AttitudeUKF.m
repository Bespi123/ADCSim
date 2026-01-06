classdef AttitudeUKF < handle
%AttitudeUKF An Unscented Kalman Filter with separate Predict and Correct steps.
%
%   This version is structured for real-time applications where sensors
%   operate at different sampling rates. The Predict step should be called
%   at the high frequency of the gyroscope, while the Correct step is called
%   at the lower frequency of the attitude sensors (accelerometer/magnetometer).
%
%   Date          Author          Notes
%   08/10/2025    bespi123        Update code to add star Sensor references.

    %% Public Properties
    properties (Access = public)
        % --- Core Filter Properties ---
        dt        double    % Sample period for the high-frequency loop (gyro).
        x_est     double    % The estimated state (attitude) quaternion [w, x, y, z]'.
        P_est     double    % The estimated state ERROR covariance matrix (3x3).
        
        % --- UKF Parameters ---
        alpha     double
        beta      double
        kappa     double
        
        % --- Noise Covariances ---
        Q         double    % Process noise covariance matrix (for gyro error, 3x3).
        R         double    % Measurement noise covariance matrix (6x6).
        
        % --- Earth Reference Vectors ---
        g_ref     double
        m_ref     double
        stars_ref double

        imu
        starTracker
    end
    
    properties (Access = private)
        % --- Internal UKF Calculation Variables ---
        lambda_param      double
        weights_m         double
        weights_c         double
        num_error_states  double
        num_sigma_points  double
    end

    %% Public Methods
    methods (Access = public)
        function obj = AttitudeUKF(ukf_params, sample_time, imu_obj, st_obj)
            %AttitudeUKF Constructor for the UKF class.            
            obj.imu = imu_obj;
            % Check if a StarTracker was provided
            if nargin > 3 && ~isempty(st_obj)
                obj.starTracker = st_obj;
            else
                obj.starTracker = struct('enable', false); % Default to disabled
            end
            obj.dt = sample_time;
            obj.x_est = [1; 0; 0; 0];
            obj.P_est = eye(3) * 0.1;
            
            obj.num_error_states = 3;
            obj.num_sigma_points = 2 * obj.num_error_states + 1;
            obj.alpha = ukf_params.alpha;
            obj.beta  = ukf_params.beta;
            obj.kappa = ukf_params.kappa;
            
            obj.Q = diag(ukf_params.gyro_std.^2);
            R_acc = diag(ukf_params.acc_std.^2);
            R_mag = diag(ukf_params.mag_std.^2);
            if obj.starTracker.enable
                R_star_single = diag(ukf_params.star_std.^2);
                R_star = kron(eye(obj.starTracker.numberOfStars), R_star_single); % Stack R for N stars
                obj.stars_ref = obj.starTracker.stars_I;
            else
                obj.stars_ref = [];
                R_star = [];
            end
            obj.R = blkdiag(R_acc, R_mag, R_star);
            
            obj.g_ref = obj.imu.g_I;
            obj.m_ref = obj.imu.m_I;
            
            obj.calculate_weights();
        end
        
        function Predict(obj, gyro)
            %PREDICT Propagates the state estimate using gyroscope data.
            %   This method should be called at a high frequency.
            
            % Generate sigma points from the current state and error covariance
            sigma_quats = obj.generate_sigma_points(obj.x_est, obj.P_est);

            % Propagate each sigma point through the process model
            propagated_quats = zeros(4, obj.num_sigma_points);
            for i = 1:obj.num_sigma_points
                propagated_quats(:, i) = obj.process_model(sigma_quats(:, i), gyro, obj.dt);
            end

            % Recover the new predicted mean and error covariance
            [x_pred, P_pred] = obj.recover_statistics_quat(propagated_quats, obj.Q);
            
            % Update the object's state and covariance with the predicted values
            obj.x_est = x_pred;
            obj.P_est = P_pred;
        end
        
        function Correct(obj, y_meas)
            %CORRECT Corrects the state estimate using accelerometer and magnetometer data.
            %   This method should be called at a lower frequency, when new
            %   attitude sensor measurements are available.
            
            % The current state obj.x_est is already the result of the last prediction
            x_pred = obj.x_est;
            P_pred = obj.P_est;
            
            % Normalize measurements
            acc = y_meas(1:3);
            mag = y_meas(4:6);
            
            if obj.starTracker.enable
                star = y_meas(7:end);
                numberOfStars = obj.starTracker.numberOfStars;
            else 
                star = [];
                numberOfStars = 0;
            end
                        
            % Generate a new set of sigma points around the predicted state
            sigma_quats = obj.generate_sigma_points(x_pred, P_pred);
            
            % Transform sigma points into measurement space
            num_meas = 6+3*numberOfStars;

            measurement_points = zeros(num_meas, obj.num_sigma_points);
            for i = 1:obj.num_sigma_points
                measurement_points(:, i) = obj.measurement_model(sigma_quats(:, i));
            end

            [z_pred, Pz] = obj.recover_statistics_vector(measurement_points, obj.R);
            
            Pxz = obj.calculate_cross_covariance(sigma_quats, x_pred, measurement_points, z_pred);
            
            K = Pxz / Pz; % Kalman Gain

            z_actual = [acc; mag; star];
            innovation = z_actual - z_pred;
            
            correction_error_vec = K * innovation; % 3x1 rotation vector
            
            % Apply the correction as a quaternion rotation
            correction_quat = obj.error_vec_to_quat(correction_error_vec);
            obj.x_est = quatmultiply(correction_quat', x_pred')';
            obj.x_est = obj.x_est / norm(obj.x_est);
            
            % Update the 3x3 error covariance
            obj.P_est = P_pred - K * Pz * K';
        end
    end
    
    %% Private Methods
    methods (Access = private)
        function calculate_weights(obj)
            n = obj.num_error_states;
            alpha_sq = obj.alpha^2;
            obj.lambda_param = alpha_sq * (n + obj.kappa) - n;
            lambda_n_sum = obj.lambda_param + n;
            obj.weights_m = zeros(1, obj.num_sigma_points);
            obj.weights_c = zeros(1, obj.num_sigma_points);
            obj.weights_m(1) = obj.lambda_param / lambda_n_sum;
            obj.weights_c(1) = obj.lambda_param / lambda_n_sum + (1 - alpha_sq + obj.beta);
            common_weight = 0.5 / lambda_n_sum;
            obj.weights_m(2:end) = common_weight;
            obj.weights_c(2:end) = common_weight;
        end

        function sigma_quats = generate_sigma_points(obj, q_mean, P_error)
            n = obj.num_error_states;
            scale = sqrt(n + obj.lambda_param);
            try
                L = scale * chol(P_error, 'lower');
            catch
                P_error_reg = P_error + eye(n)*1e-9;
                L = scale * chol(P_error_reg, 'lower');
            end
            sigma_quats = zeros(4, obj.num_sigma_points);
            sigma_quats(:, 1) = q_mean;
            error_vectors = [L, -L];
            for i = 1:(2*n)
                q_err = obj.error_vec_to_quat(error_vectors(:, i));
                sigma_quats(:, i + 1) = quatmultiply(q_err', q_mean')';
            end
        end
        
        function [mean_q, P_error] = recover_statistics_quat(obj, q_points, q_noise_cov)
            q_avg = q_points(:,1);
            error_vectors = zeros(obj.num_error_states, obj.num_sigma_points);
            for iter = 1:3
                for i = 1:obj.num_sigma_points
                    q_err = quatmultiply(q_points(:,i)', quatinv(q_avg'))';
                    error_vectors(:,i) = obj.quat_to_error_vec(q_err);
                end
                mean_error_vec = error_vectors * obj.weights_m';
                mean_error_quat = obj.error_vec_to_quat(mean_error_vec);
                q_avg = quatmultiply(mean_error_quat', q_avg')';
                q_avg = q_avg / norm(q_avg);
            end
            mean_q = q_avg;
            diff_error = error_vectors - mean_error_vec;
            P_error = (obj.weights_c .* diff_error) * diff_error' + q_noise_cov;
        end
        function [mean_val, cov_val] = recover_statistics_vector(obj, points, noise_cov)
            mean_val = points * obj.weights_m';
            diff = points - mean_val;
            cov_val = (obj.weights_c .* diff) * diff' + noise_cov;
        end
        function Pxz = calculate_cross_covariance(obj, q_points, q_mean, z_points, z_mean)
            num_meas = size(z_points, 1);
            Pxz = zeros(obj.num_error_states, num_meas);
            for i = 1:obj.num_sigma_points
                q_err = quatmultiply(q_points(:, i)', quatinv(q_mean'))';
                diff_x_error = obj.quat_to_error_vec(q_err);
                diff_z = z_points(:, i) - z_mean;
                Pxz = Pxz + obj.weights_c(i) * (diff_x_error * diff_z');
            end
        end
        function z_pred = measurement_model(obj, x)
            R_i2b = quat2dcm(x');
            acc_pred = R_i2b * obj.g_ref;
            mag_pred = R_i2b * obj.m_ref;
            if obj.starTracker.enable
                star_pred = zeros(3*obj.starTracker.numberOfStars,1);
                for k = 1:obj.starTracker.numberOfStars
                    star_pred(3*k-2:3*k,1) = R_i2b*obj.stars_ref(:,k);
                end
            else
                star_pred = [];
            end    
            z_pred = [acc_pred; mag_pred; star_pred];
        end
    end
    methods (Static, Access = private)
        function x_next = process_model(x, gyro, dt)
            omega = [0; gyro(:)];
            q = x(:);
            q_dot = 0.5 * quatmultiply(q', omega')';
            x_next = q + q_dot * dt;
            x_next = x_next / norm(x_next);
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