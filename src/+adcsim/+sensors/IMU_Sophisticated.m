classdef IMU_Sophisticated < handle
    % IMU_Sophisticated Encapsulates a high-fidelity model of an Inertial Measurement Unit.
    %   This class implements advanced error sources including bias instability (random walk), 
    %   scale factor errors, and axis misalignment as described in equations (8) to (12)[cite: 162, 163].
    %   It also features saturation and quantization to model hardware constraints as suggested by Reviewer #1.

    %% Public Properties
    properties (Access = public)
        Gyroscope       struct  % Contains Tma, Tsf, drift, and limits[cite: 33, 41].
        Accelerometer   struct  % Contains bias, noise, and limits[cite: 33, 41].
        Magnetometer    struct  % Contains bias, noise, and limits[cite: 33, 41].
        
        g_I             double  % Gravity reference vector in the inertial frame[cite: 156].
        m_I             double  % Magnetic field reference vector in the inertial frame[cite: 156].
    end
    
    %% Public Methods
    methods (Access = public)
        function obj = IMU_Sophisticated(sensors)
            % IMU_Sophisticated Constructor for the advanced IMU model[cite: 33].
            
            % --- Gyroscope Model Initialization ---
            gyro_sf_matrix = diag([1,1,1] + sensors.gyro.scale_factor_errors');
            gyro_ma_matrix = [1, sensors.gyro.misalignment_errors(1), sensors.gyro.misalignment_errors(2);
                              sensors.gyro.misalignment_errors(3), 1, sensors.gyro.misalignment_errors(4);
                              sensors.gyro.misalignment_errors(5), sensors.gyro.misalignment_errors(6), 1];
            
            obj.Gyroscope = struct(...
                'std', sensors.gyro.std, ...
                'bias', sensors.gyro.initial_bias, ... 
                'std_RW', sensors.gyro.std_RW, ...     
                'scale_factor', gyro_sf_matrix, ...
                'misalignment', gyro_ma_matrix, ...
                'tau_gyro_filter', sensors.gyro.tau, ...
                'last_filtered_omega', zeros(3,1),...
                'ts', sensors.gyro.Ts, ...
                'limit', sensors.gyro.limit, ...      
                'bits', sensors.gyro.bits ...          
            );
            
            % --- Accelerometer Model Initialization ---
            accel_sf_matrix = diag([1,1,1] + sensors.acc.scale_factor_errors');
            accel_ma_matrix = [1, sensors.acc.misalignment_errors(1), sensors.acc.misalignment_errors(2);
                               sensors.acc.misalignment_errors(3), 1, sensors.acc.misalignment_errors(4);
                               sensors.acc.misalignment_errors(5), sensors.acc.misalignment_errors(6), 1];
            obj.Accelerometer = struct(...
                'std', sensors.acc.std, ...
                'bias', sensors.acc.bias, ...
                'scale_factor', accel_sf_matrix, ...
                'misalignment', accel_ma_matrix, ...
                'limit', sensors.acc.limit, ...       
                'bits', sensors.acc.bits ...          
            );
            
            % --- Magnetometer Model Initialization ---
            mag_sf_matrix = diag([1,1,1] + sensors.mag.scale_factor_errors');
            mag_ma_matrix = [1, sensors.mag.misalignment_errors(1), sensors.mag.misalignment_errors(2);
                             sensors.mag.misalignment_errors(3), 1, sensors.mag.misalignment_errors(4);
                             sensors.mag.misalignment_errors(5), sensors.mag.misalignment_errors(6), 1];
            obj.Magnetometer = struct(...
                'std', sensors.mag.std, ...
                'bias', sensors.mag.bias, ...
                'scale_factor', mag_sf_matrix, ...
                'misalignment', mag_ma_matrix, ...
                'limit', sensors.mag.limit, ... 
                'bits', sensors.mag.bits ...          
            );
            
            obj.g_I = sensors.acc.g_I;
            obj.m_I = sensors.mag.m_I;
        end
        
        function reading = getGyroscopeReading(obj, omega)
            % getGyroscopeReading with added hardware constraints[cite: 168].
            dt = obj.Gyroscope.ts;
            distorted_omega = obj.Gyroscope.misalignment * obj.Gyroscope.scale_factor * omega;
            
            % Update drift (Random Walk)[cite: 171].
            bias_drift = obj.Gyroscope.std_RW .* sqrt(dt) .* randn(3,1);
            obj.Gyroscope.bias = obj.Gyroscope.bias + bias_drift;
            
            white_noise = obj.Gyroscope.std .* randn(3, 1);
            reading = distorted_omega + obj.Gyroscope.bias + white_noise;
            
            % Final Step: Apply Hardware Limits.
            reading = obj.applySaturationAndQuantization(reading, -obj.Gyroscope.limit, obj.Gyroscope.limit, obj.Gyroscope.bits);
        end

         function reading_filtered = filterGyroscopeReading(obj, gyro_reading)       
            %filterGyroscopeReading Applies a low-pass filter to the gyro measurement.
            %   This simulates the smoothing effect often present in real sensor hardware
            %   or preprocessing software.
            %
            %   Inputs:
            %       gyro_reading - 3x1 raw (simulated) gyroscope measurement.
            %   Outputs:
            %       reading_filtered - 3x1 filtered gyroscope measurement.
            
            % Ts_gyro: Sampling time of the gyroscope.
            Ts_gyro = obj.Gyroscope.ts;
            % Discrete first-order low-pass filter implementation.
            alpha = Ts_gyro / (obj.Gyroscope.tau_gyro_filter + Ts_gyro);
            reading_filtered = (1 - alpha) * obj.Gyroscope.last_filtered_omega + alpha * gyro_reading;
            
            % Store the current filtered value for the next iteration's calculation.
            obj.Gyroscope.last_filtered_omega = reading_filtered;
        end
        
        function reading = getAccelerometerReading(obj, R)
            % getAccelerometerReading with hardware constraints[cite: 178].
            true_reading = R' * obj.g_I;
            distorted_reading = obj.Accelerometer.misalignment * obj.Accelerometer.scale_factor * true_reading;
            white_noise = obj.Accelerometer.std .* randn(3, 1);
            reading = distorted_reading + obj.Accelerometer.bias + white_noise;
            
            % Final Step: Apply Hardware Limits.
            reading = obj.applySaturationAndQuantization(reading, -obj.Accelerometer.limit, obj.Accelerometer.limit, obj.Accelerometer.bits);
            reading = reading / (norm(reading) + eps); 
        end
        
        function reading = getMagnetometerReading(obj, R)
            % getMagnetometerReading with hardware constraints[cite: 179].
            true_reading = R' * obj.m_I;
            distorted_reading = obj.Magnetometer.misalignment * obj.Magnetometer.scale_factor * true_reading;
            white_noise = obj.Magnetometer.std .* randn(3, 1);
            reading = distorted_reading + obj.Magnetometer.bias + white_noise;
            
            % Final Step: Apply Hardware Limits.
            reading = obj.applySaturationAndQuantization(reading, -obj.Magnetometer.limit, obj.Magnetometer.limit, obj.Magnetometer.bits);
            reading = reading / (norm(reading) + eps);
        end
    end

    %% Private Methods
    methods (Access = private)
        function y = applySaturationAndQuantization(obj, y, y_min, y_max, n_bits)
            % applySaturationAndQuantization Logic.
            y = min(max(y, y_min), y_max);
            if ~isempty(n_bits) && n_bits > 0
                levels = 2^n_bits - 1;
                delta  = (y_max - y_min) ./ levels;
                delta(delta == 0) = eps;
                y = y_min + delta .* round((y - y_min) ./ delta);
            end
        end
    end
end