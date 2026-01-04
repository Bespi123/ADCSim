classdef MadgwickAHRS < handle
%MADGWICKAHRS Implementation of Madgwick's IMU and AHRS algorithms
%
%   For more information see:
%   http://www.x-io.co.uk/node/8#open_source_ahrs_and_imu_algorithms
%
%   Date          Author          Notes
%   28/09/2011    SOH Madgwick    Initial release
%   06/10/2025    Bespi123        Extention to use more reference vectors 
    %% Public properties
    properties (Access = public)
        dt        double
        x_est     double
        gyro_bias double
        bias_gain double
        Beta      double
        % References to sensor objects
        imu
        starTracker
    end
    %% Public methods
    methods (Access = public)
        function obj = MadgwickAHRS(madwick_params, sample_time, imu_obj, st_obj)
            import adcsim.utils.*
            % Save references to sensor objects
            obj.imu = imu_obj;       % Assign IMU object reference
            obj.dt  = sample_time;
            % Check if a StarTracker was provided
            if nargin > 3 && ~isempty(st_obj)
                obj.starTracker = st_obj;
            else
                obj.starTracker = struct('enable', false); % Default to disabled
            end
            % Initialize state and covariance
            obj.x_est = [1, 0, 0, 0]';      % output quaternion describing the Earth relative to the sensor
            obj.gyro_bias = zeros(3,1);     % estimated gyroscope bias
            obj.bias_gain = madwick_params.bias_gain; % bias gain
            obj.Beta = madwick_params.beta;           % algorithm gain
        end

        function obj = Predict(obj, Gyroscope)
            % Predicts the next state using only the gyroscope reading.
            % This should be called at a high frequency.
            import adcsim.utils.*
            
            q = obj.x_est;
            
            % Remove the estimated bias from the gyro reading
            Gyroscope_corrected = Gyroscope - obj.gyro_bias;
            
            % Compute the rate of change of the quaternion from angular velocity
            qDot = 0.5 * quaternProd(q', [0, Gyroscope_corrected']);
            
            % Integrate to get the new, predicted orientation
            q = q + qDot' * obj.dt;
            obj.x_est = q / norm(q); % Normalize the quaternion
        end
        function obj = Update(obj, Gyroscope, y_meas)
            import adcsim.utils.*
            q = obj.x_est; % short name local variable for readability
            acc = y_meas(1:3);
            mag = y_meas(4:6);
            
            % Reference direction of Earth's magnetic feild
            h = quaternProd(q', quaternProd([0;mag]', quaternConj(q')));
            b = [0, norm([h(2) h(3)]), 0, h(4)];
            % Calculate the objective function (error) for the accelerometer and magnetometer
            F_acc = obj.computeObjectiveFunction(obj.imu.g_I, acc); % Gravity reference vector
            F_mag = obj.computeObjectiveFunction(b(2:4), mag);
           
            % Calculate the Jacobian for the accelerometer and magnetometer
            J_acc = obj.computeJacobian(obj.imu.g_I);
            J_mag = obj.computeJacobian(b(2:4));
            if  obj.starTracker.enable
                numberOfStars = obj.starTracker.numberOfStars;
                % Use star tracker measurements to correct the orientation
                starMeasurement = y_meas(7:end); % Get star tracker measurement
                % Initialize F_star and J_star
                F_star = NaN(3*numberOfStars, 1);
                J_star = NaN(3*numberOfStars, 4);
                
                for k = 1:numberOfStars
                    star_I = obj.starTracker.stars_I(:,k);
                    starMeasurement_k = starMeasurement(3*k-2:3*k);
                    F_star(3*k-2:3*k,:) = obj.computeObjectiveFunction(star_I, starMeasurement_k); 
                    J_star(3*k-2:3*k,:) = obj.computeJacobian(star_I);
                end
            else
                F_star = []; J_star=[];
            end
            F = [F_acc; F_mag; F_star]; % Combine errors into a single vector
            J = [J_acc; J_mag; J_star]; % Combine Jacobians into a single matrix
            step = (J'*F);
            step = step / norm(step);	% normalize step magnitude
            
            % Compute bias
            omega_e_t = quaternProd(2*quaternConj(q'),step');
            obj.gyro_bias = obj.gyro_bias + obj.bias_gain * omega_e_t(2:4)'*obj.dt; % Update gyroscope bias
            % Compute rate of change of quaternion
            qDot = 0.5 * quaternProd(q', [0, Gyroscope(1)-obj.gyro_bias(1),... 
                Gyroscope(2)-obj.gyro_bias(2), Gyroscope(3)-obj.gyro_bias(3)]) - obj.Beta * step';
            % Integrate to yield quaternion
            q = q + qDot' * obj.dt;  % 
            obj.x_est = q / norm(q); % normalise quaternion
        end
    end
    methods(Access=private)
        function f = computeObjectiveFunction(obj, d, s)
            % COMPUTEOBJECTIVEFUNCTION - Calculates the objective function (error) for a reference vector.
            %   This function computes the difference between the reference vector 'd' (in the
            %   Earth frame) rotated by the current orientation 'q', and the vector measured
            %   by the sensor 's' (in the body frame).
            
            q = obj.x_est;
            q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4);
            dx = d(1); dy = d(2); dz = d(3);
            sx = s(1); sy = s(2); sz = s(3);
            f = [ 2*(dx*(0.5 - q3^2 - q4^2) + dy*(q1*q4 + q2*q3) + dz*(q2*q4 - q1*q3)) - sx;
                  2*(dx*(q2*q3 - q1*q4) + dy*(0.5 - q2^2 - q4^2) + dz*(q1*q2 + q3*q4)) - sy;
                  2*(dx*(q1*q3 + q2*q4) + dy*(q3*q4 - q1*q2) + dz*(0.5 - q2^2 - q3^2)) - sz ];
      end
      
      function J = computeJacobian(obj, d)
        %Generates the Jacobian matrix (derivative of 'f' wrt 'q').
        q = obj.x_est;
        q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4);
        dx = d(1); dy = d(2); dz = d(3);
        J = zeros(3, 4);
        J(1, 1) = 2*(dy*q4 - dz*q3);
        J(2, 1) = 2*(-dx*q4 + dz*q2);
        J(3, 1) = 2*(dx*q3 - dy*q2);
        
        J(1, 2) = 2*(dy*q3 + dz*q4);
        J(2, 2) = 2*(dx*q3 - 2*dy*q2 + dz*q1);
        J(3, 2) = 2*(dx*q4 - dy*q1 - 2*dz*q2);
        
        J(1, 3) = 2*(-2*dx*q3 + dy*q2 - dz*q1);
        J(2, 3) = 2*(dx*q2 + dz*q4);
        J(3, 3) = 2*(dx*q1 + dy*q4 - 2*dz*q3);
        
        J(1, 4) = 2*(-2*dx*q4 + dy*q1 + dz*q2);
        J(2, 4) = 2*(-dx*q1 - 2*dy*q4 + 2*dz*q3);
        J(3, 4) = 2*(dx*q2 + dy*q3);
      end
    end
end