classdef Mahony
    %MAHONY  Nonlinear Mahony complementary attitude filter.
    %
    % This class implements a Mahony-style nonlinear complementary filter
    % for attitude estimation using quaternions. The filter fuses gyroscope
    % measurements with vector observations (accelerometer, magnetometer,
    % and optionally star tracker) and includes integral estimation of the
    % gyroscope bias.
    %
    % The implementation follows the philosophy of:
    %   Mahony, Hamel, Pflimlin (2008),
    %   "Nonlinear Complementary Filters on the Special Orthogonal Group".
    %
    % IMPORTANT CONVENTION NOTES:
    % ------------------------------------------------------------------
    % - x_est is a unit quaternion representing the attitude estimate.
    % - quat2rot(q) is assumed to return R_IB:
    %       rotation from BODY frame to INERTIAL frame.
    % - Therefore, Rhat' = R_BI is used to map inertial reference vectors
    %   into the body frame for comparison with measurements.
    %
    % - All measured vectors y_meas are assumed to be expressed in the
    %   BODY frame.
    % - All reference vectors y_ref are expressed in the INERTIAL frame.
    %
    % Frame consistency is CRITICAL. If the filter diverges with perfect
    % data, the first thing to check is the rotation convention and the
    % transpose used when computing predicted vectors.
    %
    % STATE:
    %   x_est     (4x1) quaternion attitude estimate
    %   gyro_bias (3x1) estimated gyroscope bias [rad/s]
    %
    % MAIN METHODS:
    %   Predict() : high-rate propagation using gyro only
    %   Update()  : correction using vector measurements + bias update
    %
    % WARNING:
    %   Do NOT call Predict() and Update() with full integration at the
    %   same rate unless you know what you are doing, otherwise the gyro
    %   may be integrated twice per time step.

    properties
        x_est        double   % Estimated attitude quaternion (4x1)
        K_P          double   % Proportional gain (attitude correction)
        K_I          double   % Integral gain (gyro bias estimation)
        K            double   % Measurement reliability weights (1xn)
        y_ref        double   % Reference inertial vectors (3xn)
        dt           double   % Filter sampling time [s]
        gyro_bias    double   % Estimated gyroscope bias (3x1)

        % Sensor object references
        imu                 % IMU object (must provide g_I, m_I in inertial frame)
        starTracker         % StarTracker object (optional)
    end

    methods (Access = public)

        function obj = Mahony(mahony_params, sample_time, imu_obj, st_obj)
            %MAHONY Constructor.
            %
            % Inputs:
            %   mahony_params : struct containing filter gains
            %       - K_P : proportional gain
            %       - K_I : integral gain
            %       - K_1 : accelerometer weight
            %       - K_2 : magnetometer weight
            %       - K_3 : star tracker weight
            %
            %   sample_time  : filter sampling time [s]
            %   imu_obj      : IMU object providing inertial reference vectors
            %   st_obj       : StarTracker object (optional, may be empty)
            %
            % The constructor initializes the filter state, stores sensor
            % references, and builds the set of inertial reference vectors
            % and their corresponding confidence weights.

            % Save sensor references
            obj.imu = imu_obj;
            obj.dt  = sample_time;

            % Star tracker handling
            if nargin > 3 && ~isempty(st_obj)
                obj.starTracker = st_obj;
            else
                obj.starTracker = struct('enable', false);
            end

            % Initial state
            obj.x_est     = [1; 0; 0; 0];   % Identity quaternion
            obj.gyro_bias = zeros(3,1);     % Initial gyro bias estimate

            % Gains
            obj.K_P = mahony_params.K_P;
            obj.K_I = mahony_params.K_I;

            % Build reference vector set and weights
            if obj.starTracker.enable
                % Inertial references: gravity, magnetic field, star vectors
                obj.y_ref = [obj.imu.g_I, ...
                             obj.imu.m_I, ...
                             obj.starTracker.stars_I];

                % Corresponding weights
                obj.K = [mahony_params.K_1, ...
                         mahony_params.K_2, ...
                         mahony_params.K_3 * ones(1, obj.starTracker.numberOfStars)];
            else
                % Only gravity and magnetic field
                obj.y_ref = [obj.imu.g_I, obj.imu.m_I];
                obj.K     = [mahony_params.K_1, mahony_params.K_2];
            end
        end

        function obj = Predict(obj, Gyroscope)
            %PREDICT High-rate attitude propagation using gyro only.
            %
            % This method integrates the quaternion kinematics using the
            % bias-corrected gyroscope measurement. It should be called at
            % high frequency (IMU rate).
            %
            % Input:
            %   Gyroscope : measured angular velocity (3x1) [rad/s], body frame

            import adcsim.utils.*

            % Current attitude estimate
            q = obj.x_est;

            % Bias-corrected angular rate
            Gyroscope_corrected = Gyroscope - obj.gyro_bias;

            % Quaternion kinematics:
            %   q_dot = 0.5 * q ⊗ [0; omega]
            qDot = 0.5 * quaternProd(q', [0, Gyroscope_corrected']);

            % Euler integration
            q = q + qDot' * obj.dt;

            % Normalize to enforce unit quaternion constraint
            obj.x_est = q / norm(q);
        end

        function obj = Update(obj, Gyroscope, y_meas)
            %UPDATE Measurement-based correction and bias estimation.
            %
            % This method computes the Mahony attitude error using vector
            % measurements and applies proportional-integral feedback.
            %
            % Inputs:
            %   Gyroscope : measured angular velocity (3x1) [rad/s]
            %   y_meas    : stacked measurement vectors (3*n x 1), each
            %               vector expressed in the BODY frame
            %
            % The error is constructed as:
            %   S = Σ (k_i/2) * (v_i * v̂_i' - v̂_i * v_i')
            %   e = -vex(S)
            %
            % where:
            %   v_i     : measured direction (body frame)
            %   v̂_i     : predicted direction (body frame)
            %
            % NOTE:
            %   The minus sign in e = -vex(S) is consistent with the chosen
            %   quaternion kinematics and rotation convention.

            import adcsim.utils.*

            % Accumulate skew-symmetric error matrix
            S = zeros(3,3);

            % Convert quaternion to rotation matrix
            % quat2rot() is assumed to return R_IB (body -> inertial)
            Rhat = quat2rot(obj.x_est');

            % Number of reference vectors
            n = size(obj.y_ref, 2);

            for i = 1:n
                % Extract measured vector (body frame)
                v_i = y_meas(3*(i-1)+1 : 3*(i-1)+3);

                % Predicted vector:
                %   v̂_i = R_BI * v_I = Rhat' * y_ref
                hat_v_i = Rhat' * obj.y_ref(:,i);

                % Normalize vectors
                v_i     = v_i / norm(v_i);
                hat_v_i = hat_v_i / norm(hat_v_i);

                % Mahony error contribution (skew-symmetric form)
                S = S + (obj.K(i)/2) * (v_i*hat_v_i' - hat_v_i*v_i');
            end

            % Map skew-symmetric matrix to vector error
            e = -obj.vex(S);

            % Gyroscope bias estimator:
            %   b_dot = -K_I * e
            hat_b_dot = -obj.K_I * e;
            obj.gyro_bias = obj.gyro_bias + hat_b_dot * obj.dt;

            % Corrected angular velocity
            omega_corr = Gyroscope - obj.gyro_bias + obj.K_P * e;

            % Quaternion correction step
            hat_q_dot = 0.5 * quaternProd(obj.x_est', [0; omega_corr]');

            % Integrate and normalize
            obj.x_est = obj.x_est + hat_q_dot' * obj.dt;
            obj.x_est = obj.x_est / norm(obj.x_est);
        end
    end

    methods (Access = private)
        function x = vex(~, S)
            %VEX Convert a skew-symmetric matrix to a vector.
            %
            % For:
            %   S = [  0   -x3   x2;
            %          x3   0   -x1;
            %         -x2   x1   0 ]
            %
            % vex(S) = [x1; x2; x3]

            x = [ S(3,2);
                  S(1,3);
                  S(2,1) ];
        end
    end
end
