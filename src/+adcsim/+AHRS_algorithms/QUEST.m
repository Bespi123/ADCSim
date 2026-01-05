classdef QUEST
    %QUEST  Quaternion ESTimator (QUEST) attitude solver using Davenport's K.
    %
    % This class computes an attitude quaternion from a set of vector
    % observations by solving Wahba's problem. The implementation uses
    % Davenport's K-matrix and extracts the quaternion as the eigenvector
    % associated with the largest eigenvalue (a common and robust QUEST form).
    %
    % INPUT CONVENTIONS (CRITICAL):
    % -------------------------------------------------------------------------
    % - b_i (measurements) are vectors expressed in the BODY frame.
    % - r_i (references) are the corresponding known vectors expressed in the
    %   INERTIAL frame.
    % - Each (b_i, r_i) pair must represent the SAME physical direction.
    %   Examples:
    %     * Accelerometer (gravity direction): b_g in body, r_g in inertial
    %     * Magnetometer (mag field direction): b_m in body, r_m in inertial
    %     * Star tracker: body-frame star direction(s), inertial catalog direction(s)
    %
    % OUTPUT CONVENTION:
    % -------------------------------------------------------------------------
    % The eigenvector of Davenport's K is returned here as [q_vec; q_scalar]
    % (vector part first). This class converts it to scalar-first:
    %   x_est = [q0; q1; q2; q3]
    %
    % NOTE ON GYROSCOPE:
    % -------------------------------------------------------------------------
    % QUEST itself is a "measurement-only" attitude determination algorithm.
    % It does NOT require gyroscope measurements. In this class, Gyroscope is
    % accepted in Update() only to keep the interface consistent with other
    % estimators (e.g., Mahony/EKF). It is currently not used.
    %
    % NUMERICAL NOTE:
    % -------------------------------------------------------------------------
    % The K matrix is symmetrized as (K+K')/2 before eigen decomposition to
    % reduce numerical asymmetries.

    properties
        x_est        double   % Estimated attitude quaternion (4x1), scalar-first [q0;q1;q2;q3]
        dt           double   % Sampling time [s] (not used in pure QUEST update)
        gyro_bias    double   % Placeholder for gyro bias (3x1). Not used by QUEST update here.
        K            double   % Measurement reliability weights a_i (1xn)
        y_ref        double   % Reference inertial vectors r_i (3xn), columns are reference directions

        % Sensor object references
        imu                 % IMU object (must provide g_I, m_I in inertial frame)
        starTracker         % StarTracker object (optional; must provide enable, stars_I, numberOfStars)
    end

    methods
        function obj = QUEST(quest_params, sample_time, imu_obj, st_obj)
            %QUEST Constructor.
            %
            % Inputs:
            %   quest_params : struct with fields:
            %       K_1 : weight for accelerometer/gravity vector
            %       K_2 : weight for magnetometer vector
            %       K_3 : weight for star tracker vectors (applied to each star)
            %   sample_time : dt [s]
            %   imu_obj     : IMU object providing inertial reference vectors:
            %                imu.g_I, imu.m_I (3x1 each, should be unit or close)
            %   st_obj      : StarTracker object (optional). If empty, star tracker is disabled.
            %
            % This constructor builds:
            %   - y_ref : stacked inertial reference vectors r_i (3xn)
            %   - K     : associated weights a_i (1xn)
            % matching the measurement stacking expected in Update().

            % Save sensor references
            obj.imu = imu_obj;
            obj.dt  = sample_time;

            % Star tracker handling
            if nargin > 3 && ~isempty(st_obj)
                obj.starTracker = st_obj;
            else
                obj.starTracker = struct('enable', false);
            end

            % Initial state (identity attitude)
            obj.x_est     = [1; 0; 0; 0];
            obj.gyro_bias = zeros(3,1); % Not used by QUEST in this implementation

            % Build inertial reference vectors and weights
            if obj.starTracker.enable
                % Inertial references: gravity, magnetic field, and star catalog vectors
                obj.y_ref = [obj.imu.g_I, ...
                             obj.imu.m_I, ...
                             obj.starTracker.stars_I];

                % Corresponding weights (a_i)
                obj.K = [quest_params.K_1, ...
                         quest_params.K_2, ...
                         quest_params.K_3 * ones(1, obj.starTracker.numberOfStars)];
            else
                % Only gravity and magnetic field
                obj.y_ref = [obj.imu.g_I, obj.imu.m_I];
                obj.K     = [quest_params.K_1, quest_params.K_2];
            end
        end

        function obj = Update(obj, Gyroscope, y_meas)
            %UPDATE Solve Wahba's problem and update quaternion estimate.
            %
            % Inputs:
            %   Gyroscope : (3x1) angular rate [rad/s] (NOT used in this QUEST update)
            %   y_meas    : stacked body-frame measurement vectors (3*n x 1):
            %              [b_1; b_2; ...; b_n], each b_i is 3x1
            %
            % Assumptions:
            %   - Each b_i corresponds to r_i = y_ref(:,i)
            %   - y_meas is provided in BODY frame
            %   - y_ref is provided in INERTIAL frame
            %
            % Algorithm (Davenport K / QUEST via eigen):
            %   m = Σ a_i
            %   B = (1/m) Σ a_i b_i r_i^T
            %   z = (1/m) Σ a_i (b_i × r_i)
            %   S = B + B^T
            %   σ = tr(B)  (equivalently (1/m) Σ a_i b_i^T r_i)
            %   K = [S - σI, z; z^T, σ]
            %   q = eigenvector of largest eigenvalue of K
            %   convert to scalar-first and normalize

            %#ok<NASGU>  % Gyroscope currently unused (kept for interface consistency)

            b_list = y_meas;   % stacked measurements b_i in BODY frame
            r_list = obj.y_ref; % references r_i in INERTIAL frame
            a_list = obj.K;     % weights a_i

            n = size(r_list,2);
            m_k = sum(a_list);  % total weight

            % Accumulators for Wahba/QUEST terms
            Bsum = zeros(3,3);
            zsum = zeros(3,1);

            % (Optional) σ can be computed as trace(B) later; kept here for clarity
            for i = 1:n
                % Extract 3x1 measured vector b_i from stacked array
                b_i = b_list(3*(i-1)+1: 3*(i-1)+3);

                %%% Normalize b_i defensively (recommended in real data)
                %%%b_i = b_i / norm(b_i);

                % Reference vector r_i (often already unit; normalize defensively if needed)
                r_i = r_list(:,i);
                %%%r_i = r_i / norm(r_i);

                % Wahba accumulators
                Bsum = Bsum + a_list(i) * (b_i * r_i');      % Σ a_i b_i r_i^T
                zsum = zsum + a_list(i) * cross(b_i, r_i);   % Σ a_i (b_i × r_i)
            end

            % Normalize
            B = (1/m_k) * Bsum;
            z = (1/m_k) * zsum;

            % Symmetric part and trace
            S = B + B';
            sigma = trace(B);

            % Davenport K matrix (4x4)
            Kmat = [S - sigma*eye(3), z;
                    z',              sigma];

            % Enforce symmetry for numerical robustness
            Kmat = (Kmat + Kmat')/2;

            % Largest eigenvector corresponds to optimal quaternion
            [V,D] = eig(Kmat);
            [~,idx] = max(diag(D));
            q = V(:,idx);  % q = [q_vec; q_scalar] (vector-first in this construction)

            % Normalize and enforce sign convention (optional)
            q = q / norm(q);
            if q(4) < 0
                q = -q;
            end

            % Convert to scalar-first: [q0; q1; q2; q3]
            obj.x_est = [q(4); q(1:3)];
        end
    end
end
