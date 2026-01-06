classdef IMU_Sophisticated < handle
%IMU_Sophisticated Encapsulates a sophisticated model of an Inertial Measurement Unit.
%   This class models an IMU with advanced error sources, including:
%   - Bias Instability (modeled as a random walk).
%   - Scale Factor errors.
%   - Axis Misalignment (non-orthogonality).
%   This provides a more realistic simulation of real-world MEMS sensors.
%   Date          Author          Notes
%   08/20/2025    Bespi123        Create simple IMU  model
    %% Public Properties
    properties
        Gyroscope       struct  % Models the gyroscope with advanced errors.
        Accelerometer   struct  % Models the accelerometer with advanced errors.
        Magnetometer    struct  % Models the magnetometer with advanced errors.
        
        g_I             double  % Gravity reference vector in the inertial frame.
        m_I             double  % Magnetic reference vector in the inertial frame.
    end
    
    methods
        function obj = IMU_Sophisticated(sensors)
            %IMU_Sophisticated Constructor for the advanced IMU model.
            %   Initializes sensor models with parameters for noise, bias,
            %   scale factor, misalignment, and bias instability.
            
            % --- Gyroscope Model Initialization ---
            % Build Scale Factor and Misalignment matrices from error vectors.
            gyro_sf_matrix = diag([1,1,1] + sensors.gyro.scale_factor_errors');
            gyro_ma_matrix = [1, sensors.gyro.misalignment_errors(1), sensors.gyro.misalignment_errors(2);
                              sensors.gyro.misalignment_errors(3), 1, sensors.gyro.misalignment_errors(4);
                              sensors.gyro.misalignment_errors(5), sensors.gyro.misalignment_errors(6), 1];
            
            obj.Gyroscope = struct(...
                'std', sensors.gyro.std, ...
                'bias', sensors.gyro.initial_bias, ... % This bias will now drift over time.
                'std_RW', sensors.gyro.std_RW, ...     % Std. dev. for the bias random walk.
                'scale_factor', gyro_sf_matrix, ...
                'misalignment', gyro_ma_matrix, ...
                'tau_gyro_filter', sensors.gyro.tau, ...
                'last_filtered_omega', zeros(3,1),...
                'ts', sensors.gyro.Ts ...
            );
            
            % --- Accelerometer Model Initialization ---
            accel_sf_matrix = diag([1,1,1] + sensors.acc.scale_factor_errors');
            accel_ma_matrix = [1, sensors.acc.misalignment_errors(1), sensors.acc.misalignment_errors(2);
                               sensors.acc.misalignment_errors(3), 1, sensors.acc.misalignment_errors(4);
                               sensors.acc.misalignment_errors(5), sensors.acc.misalignment_errors(6), 1];

            obj.Accelerometer = struct(...
                'std', sensors.acc.std, ...
                'bias', sensors.acc.bias, ... % Bias is constant for accelerometer.
                'scale_factor', accel_sf_matrix, ...
                'misalignment', accel_ma_matrix ...
            );
            
            % --- Magnetometer Model Initialization ---
            mag_sf_matrix = diag([1,1,1] + sensors.mag.scale_factor_errors');
            mag_ma_matrix = [1, sensors.mag.misalignment_errors(1), sensors.mag.misalignment_errors(2);
                             sensors.mag.misalignment_errors(3), 1, sensors.mag.misalignment_errors(4);
                             sensors.mag.misalignment_errors(5), sensors.mag.misalignment_errors(6), 1];

            obj.Magnetometer = struct(...
                'std', sensors.mag.std, ...
                'bias', sensors.mag.bias, ... % Bias is constant for magnetometer.
                'scale_factor', mag_sf_matrix, ...
                'misalignment', mag_ma_matrix ...
            );
            
            % Store the inertial frame reference vectors.
            obj.g_I = sensors.acc.g_I;
            obj.m_I = sensors.mag.m_I;
        end
        
        function reading = getGyroscopeReading(obj, omega, dt)
            %getGyroscopeReading Simulates a gyroscope reading with a sophisticated error model.
            %   Includes scale factor, misalignment, a time-varying bias (random walk), and white noise.
            
            % 1. Distort the true angular velocity with scale factor and misalignment.
            distorted_omega = obj.Gyroscope.misalignment * obj.Gyroscope.scale_factor * omega;
            
            % 2. Update the gyroscope bias using a random walk model.
            %    The bias drifts over time based on std_RW (rad/s/sqrt(Hz)).
            bias_drift = obj.Gyroscope.std_RW .* sqrt(dt) .* randn(3,1);
            obj.Gyroscope.bias = obj.Gyroscope.bias + bias_drift;
            
            % 3. Add the current (time-varying) bias and Gaussian white noise.
            white_noise = obj.Gyroscope.std .* randn(3, 1);
            reading = distorted_omega + obj.Gyroscope.bias + white_noise;
        end
        
        function reading_filtered = filterGyroscopeReading(obj, gyro_reading, Ts_gyro)
            %filterGyroscopeReading Applies a low-pass filter to the gyro measurement.
            % (This method remains unchanged).
            alpha = Ts_gyro / (obj.Gyroscope.tau_gyro_filter + Ts_gyro);
            reading_filtered = (1 - alpha) * obj.Gyroscope.last_filtered_omega + alpha * gyro_reading;
            obj.Gyroscope.last_filtered_omega = reading_filtered;
        end
        
        function reading = getAccelerometerReading(obj, R)
            %getAccelerometerReading Simulates an accelerometer reading with scale factor and misalignment.
            
            % 1. Get the "true" reading by rotating the reference vector.
            true_reading = R' * obj.g_I;
            
            % 2. Distort the true reading with scale factor and misalignment.
            distorted_reading = obj.Accelerometer.misalignment * obj.Accelerometer.scale_factor * true_reading;

            % 3. Add the constant bias and Gaussian white noise.
            white_noise = obj.Accelerometer.std .* randn(3, 1);
            reading = distorted_reading + obj.Accelerometer.bias + white_noise;
            
            % 4. Normalize the final vector.
            reading = reading / norm(reading);
        end
        
        function reading = getMagnetometerReading(obj, R)
            %getMagnetometerReading Simulates a magnetometer reading with scale factor and misalignment.

            % 1. Get the "true" reading by rotating the reference vector.
            true_reading = R' * obj.m_I;

            % 2. Distort the true reading with scale factor and misalignment.
            distorted_reading = obj.Magnetometer.misalignment * obj.Magnetometer.scale_factor * true_reading;

            % 3. Add the constant bias and Gaussian white noise.
            white_noise = obj.Magnetometer.std .* randn(3, 1);
            reading = distorted_reading + obj.Magnetometer.bias + white_noise;
            
            % 4. Normalize the final vector.
            reading = reading / norm(reading);
        end
    end
end