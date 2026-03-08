classdef IMU_Sophisticated < handle
    % IMU_Sophisticated Encapsulates a high-fidelity model of an Inertial Measurement Unit.
    %   This class implements advanced error sources including bias instability (random walk), 
    %   scale factor errors, and axis misalignment as described in equations (8) to (12).
    %   It also features saturation and ADC quantization (configurable v_ref and bits) to model hardware constraints.
    
    %% Public Properties
    properties (Access = public)
        Gyroscope       struct  
        Accelerometer   struct  
        Magnetometer    struct  
        
        g_I             double  % Gravity reference vector in the inertial frame
        m_I             double  % Magnetic field reference vector in the inertial frame
    end
    
    %% Public Methods
    methods (Access = public)
        function obj = IMU_Sophisticated(sensors)
            % IMU_Sophisticated Constructor for the advanced IMU model
            
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
                'bits', sensors.gyro.bits, ...          % User configurable bit width
                'v_ref', sensors.gyro.v_ref ...         % User configurable reference voltage
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
                'bits', sensors.acc.bits, ...          
                'v_ref', sensors.acc.v_ref ...         
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
                'bits', sensors.mag.bits, ...          
                'v_ref', sensors.mag.v_ref ...         
            );
            
            obj.g_I = sensors.acc.g_I;
            obj.m_I = sensors.mag.m_I;
        end
        
        function reading = getGyroscopeReading(obj, omega)
            dt = obj.Gyroscope.ts;
            distorted_omega = obj.Gyroscope.misalignment * obj.Gyroscope.scale_factor * omega;
            
            bias_drift = obj.Gyroscope.std_RW .* sqrt(dt) .* randn(3,1);
            obj.Gyroscope.bias = obj.Gyroscope.bias + bias_drift;
            
            white_noise = obj.Gyroscope.std .* randn(3, 1);
            reading = distorted_omega + obj.Gyroscope.bias + white_noise;
            
            % Final Step: Apply ADC Conversion (Saturation + Quantization)
            reading = obj.applyADCConversion(reading, -obj.Gyroscope.limit, obj.Gyroscope.limit, obj.Gyroscope.bits, obj.Gyroscope.v_ref);
        end
        
        function reading_filtered = filterGyroscopeReading(obj, gyro_reading)       
            Ts_gyro = obj.Gyroscope.ts;
            alpha = Ts_gyro / (obj.Gyroscope.tau_gyro_filter + Ts_gyro);
            reading_filtered = (1 - alpha) * obj.Gyroscope.last_filtered_omega + alpha * gyro_reading;
            obj.Gyroscope.last_filtered_omega = reading_filtered;
        end
        
        function reading = getAccelerometerReading(obj, R)
            true_reading = R' * obj.g_I;
            distorted_reading = obj.Accelerometer.misalignment * obj.Accelerometer.scale_factor * true_reading;
            white_noise = obj.Accelerometer.std .* randn(3, 1);
            reading = distorted_reading + obj.Accelerometer.bias + white_noise;
            
            reading = obj.applyADCConversion(reading, -obj.Accelerometer.limit, obj.Accelerometer.limit, obj.Accelerometer.bits, obj.Accelerometer.v_ref);
            reading = reading / (norm(reading) + eps); 
        end
        
        % NOTE: Dynamic m_I argument added and norm() removed to support magnetorquer momentum dumping.
        function reading = getMagnetometerReading(obj, R, current_m_I)
            if nargin < 3
                current_m_I = obj.m_I;
            end
            
            true_reading = R' * current_m_I;
            distorted_reading = obj.Magnetometer.misalignment * obj.Magnetometer.scale_factor * true_reading;
            white_noise = obj.Magnetometer.std .* randn(3, 1);
            reading = distorted_reading + obj.Magnetometer.bias + white_noise;
            
            reading = obj.applyADCConversion(reading, -obj.Magnetometer.limit, obj.Magnetometer.limit, obj.Magnetometer.bits, obj.Magnetometer.v_ref);
            % NO NORMALIZATION HERE. Magnetorquers need the true physical magnitude.
        end
    end
    
    %% Private Methods
    methods (Access = private)
        function y_quantized = applyADCConversion(obj, y, y_min, y_max, n_bits, v_ref)
            % applyADCConversion Simulates the hardware ADC process.
            
            % 1. Sensor Physical Saturation
            y_sat = min(max(y, y_min), y_max);
            
            if ~isempty(n_bits) && n_bits > 0
                % 2. Map physical value to analog voltage range [0, v_ref]
                % V_in = ((y - y_min) / (y_max - y_min)) * V_ref
                range_y = y_max - y_min;
                if range_y == 0; range_y = eps; end % Prevent division by zero
                
                V_in = ((y_sat - y_min) ./ range_y) .* v_ref;
                
                % 3. ADC Digital Output Calculation (Integer digital number)
                levels = 2^n_bits - 1;
                D_out = round((V_in ./ v_ref) .* levels);
                
                % 4. Reconstruct quantized physical value for the simulation
                y_quantized = y_min + (D_out ./ levels) .* range_y;
            else
                % If no bits are specified, just return saturated analog value
                y_quantized = y_sat;
            end
        end
    end
end