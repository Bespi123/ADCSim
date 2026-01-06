classdef REQUEST
    %REQUEST  Recursive QUEST (REQUEST) attitude estimator.
    %
    % This class implements a practical REQUEST loop:
    %   1) High-rate propagation using gyro:   K <- Phi*K*Phi'
    %   2) Low-rate measurement update using vector observations (acc/mag/stars):
    %        K <- (P*m/(P*m+dm))*K + (1/(P*m+dm))*dK
    %        m <- m + dm
    %   3) Attitude quaternion is extracted as the principal eigenvector of K.
    %
    % Conventions:
    % - Quaternion output is scalar-first: x_est = [q0; q1; q2; q3]
    % - Vector measurements b_i are in BODY frame (stacked in y_meas).
    % - Reference vectors r_i are in INERTIAL frame (columns of y_ref).
    %
    % Notes:
    % - REQUEST "pure" does NOT estimate gyro bias. gyro_bias is kept here
    %   only for interface compatibility; set to zeros unless you add a bias
    %   estimator externally.
    %
    % Expected y_meas stacking (3*n x 1):
    %   y_meas = [b1; b2; ...; bn]
    % where b_i corresponds to y_ref(:,i) and has weight a(i).

    properties
        % Outputs / state
        x_est      double   % Estimated attitude quaternion (4x1), scalar-first
        gyro_bias  double   % (Optional) gyro bias (3x1). Not estimated by REQUEST here.

        % REQUEST internal memory
        K_acc      double   % Accumulated Davenport K matrix (4x4), K_{k/k}
        m_acc      double   % Accumulated total weight m_k (scalar)
        P          double   % Fading memory factor (0<P<=1). Use 1 for no fading.

        % Measurement model
        a          double   % Weights a_i (1xn)
        y_ref      double   % Inertial reference vectors r_i (3xn)

        % Timing
        dt         double   % Sample time [s] used for propagation step
        gyro_dt    double   % Sample time [s] used for propagation step (gyro)

        % Sensor references (optional, used only to build y_ref/a conveniently)
        imu
        starTracker
    end

    methods
        function obj = REQUEST(request_params, sample_time, imu_obj, st_obj)
            %REQUEST Constructor.
            %
            % request_params fields (suggested):
            %   .P   : fading memory factor (0<P<=1)
            %   .K_1 : weight for gravity (accelerometer)
            %   .K_2 : weight for magnetic field (magnetometer)
            %   .K_3 : weight for star vectors (star tracker), applied per star
            %
            % imu_obj must provide:
            %   imu_obj.g_I (3x1) gravity direction in inertial frame
            %   imu_obj.m_I (3x1) magnetic field direction in inertial frame
            %
            % st_obj (optional) must provide:
            %   st_obj.enable (logical)
            %   st_obj.stars_I (3xNstars) star inertial directions
            %   st_obj.numberOfStars

            obj.dt = sample_time;

            if isfield(request_params, 'P')
                obj.P = request_params.P;
            else
                obj.P = 1.0; % default: no fading
            end

            % Save sensor references
            obj.imu = imu_obj;

            if nargin > 3 && ~isempty(st_obj)
                obj.starTracker = st_obj;
            else
                obj.starTracker = struct('enable', false);
            end

            % Build reference vectors and weights
            if obj.starTracker.enable
                obj.y_ref = [obj.imu.g_I, obj.imu.m_I, obj.starTracker.stars_I];
                obj.a     = [request_params.K_1, request_params.K_2, ...
                             request_params.K_3 * ones(1, obj.starTracker.numberOfStars)];
            else
                obj.y_ref = [obj.imu.g_I, obj.imu.m_I];
                obj.a     = [request_params.K_1, request_params.K_2];
            end

            % Initialize REQUEST memory
            obj.K_acc = eye(4);
            obj.m_acc = 0;

            % Initialize outputs
            obj.x_est = [1;0;0;0];
            obj.gyro_bias = zeros(3,1);
        end

        function obj = Predict(obj, Gyroscope)
            %PREDICT High-rate propagation using gyro (REQUEST step 3/4).
            %
            % Propagates the accumulated Davenport matrix:
            %   K <- Phi*K*Phi'
            %
            % Gyroscope is assumed to be the angular rate over the interval dt.
            % If you have an external bias estimate, subtract it before calling.

            omega = Gyroscope - obj.gyro_bias;

            % Angular rate matrix Omega (as per your paper/figure)
            Omega = [ ...
                 0        ,  omega(3), -omega(2),  omega(1);
                -omega(3) ,  0       ,  omega(1),  omega(2);
                 omega(2) , -omega(1),  0       ,  omega(3);
                -omega(1) , -omega(2), -omega(3),  0        ];

            Phi = expm(0.5 * Omega * obj.dt);

            % If not initialized with any measurements yet, keep K_acc as identity.
            % (You could also propagate identity, which stays identity.)
            if obj.m_acc > 0
                obj.K_acc = Phi * obj.K_acc * Phi';
                obj.K_acc = (obj.K_acc + obj.K_acc')/2;
            end

            % Optionally update x_est from propagated K (not strictly required each step)
            obj.x_est = obj.extractQuatFromK(obj.K_acc);
        end

        function obj = Update(obj, Gyroscope, y_meas)
            %UPDATE Low-rate measurement update using new vector measurements (REQUEST step 5).
            %
            % Inputs:
            %   Gyroscope : (3x1) rate (used only if you want to propagate here too;
            %              recommended to call Predict() separately at high rate)
            %   y_meas    : stacked body-frame measurement vectors (3*n x 1):
            %              [b1; b2; ...; bn]
            %
            % This method:
            %   - builds delta_K from the CURRENT measurement batch
            %   - applies fading-memory update to the accumulated K_acc
            %   - extracts quaternion x_est

            %#ok<NASGU> Gyroscope  % not used here if Predict() is called externally

            b_list = y_meas;     % stacked measurements b_i in BODY frame
            r_list = obj.y_ref;  % references r_i in INERTIAL frame
            a_list = obj.a;      % weights a_i

            n = size(r_list,2);

            % δm2 = sum(a_i) for this batch
            dm = sum(a_list);

            % Accumulators (UN-normalized) for delta terms
            Bsum = zeros(3,3);
            zsum = zeros(3,1);

            for i = 1:n
                b_i = b_list(3*(i-1)+1 : 3*(i-1)+3);
                r_i = r_list(:,i);

                % Optional normalization (recommended with real sensors)
                % nb = norm(b_i); if nb > 0, b_i = b_i/nb; end
                % nr = norm(r_i); if nr > 0, r_i = r_i/nr; end

                Bsum = Bsum + a_list(i) * (b_i * r_i');
                zsum = zsum + a_list(i) * cross(b_i, r_i);
            end

            % Build delta_K (incremental Davenport matrix) using unnormalized sums:
            dB = Bsum;
            dz = zsum;
            dS = dB + dB';
            dsigma = trace(dB);

            dK = [dS - dsigma*eye(3), dz;
                  dz',               dsigma];
            dK = (dK + dK')/2;

            % REQUEST fading-memory update
            if obj.m_acc == 0
                % First measurement batch initializes K_acc as normalized dK
                obj.K_acc = (1/dm) * dK;
                obj.m_acc = dm;
            else
                denom = obj.P*obj.m_acc + dm;
                obj.K_acc = (obj.P*obj.m_acc/denom) * obj.K_acc + (1/denom) * dK;
                obj.m_acc = obj.m_acc + dm;
            end

            obj.K_acc = (obj.K_acc + obj.K_acc')/2;

            % Extract quaternion estimate
            obj.x_est = obj.extractQuatFromK(obj.K_acc);
        end
    end

    methods (Access = private)
        function q = extractQuatFromK(~, Kmat)
            % Extract principal eigenvector of Davenport K.
            % Returns scalar-first quaternion [q0;q1;q2;q3].
            [V,D] = eig(Kmat);
            [~,idx] = max(diag(D));
            qv = V(:,idx);              % [q_vec; q_scalar]
            qv = qv / norm(qv);
            if qv(4) < 0
                qv = -qv;
            end
            q = [qv(4); qv(1:3)];
        end
    end
end
