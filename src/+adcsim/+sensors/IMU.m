classdef IMU < handle
    %IMU Encapsulates the sensors of an Inertial Measurement Unit.
    %   This class models an IMU, including a gyroscope, accelerometer,
    %   and magnetometer. It simulates realistic sensor readings by adding
    %   noise and bias to the true values.
    %   Date          Author          Notes
    %   08/20/2025    Bespi123        Create simple IMU model
    %% Public Properties
    properties
        Gyroscope       struct  % Contains gyroscope parameters (std, bias, filter state).
        Accelerometer   struct  % Contains accelerometer parameters (std, bias).
        Magnetometer    struct  % Contains magnetometer parameters (std, bias).
        g_I             double  % Gravity reference vector in the inertial frame.
        m_I             double  % Magnetic reference vector in the inertial frame.
    end
    
    methods
        function obj = IMU(sensors)
            %IMU Constructor method for the IMU class.
            %   Initializes the sensor models with parameters provided in a struct.
            %
            %   Inputs:
            %       sensors - A struct containing sensor parameters (.gyro, .acc, .mag).
            %       g_I     - 3x1 gravity reference vector (e.g., [0; 0; 1]).
            %       m_I     - 3x1 magnetic reference vector.
            
            % Initialize the Gyroscope model properties.
            obj.Gyroscope = struct(...
                'std', sensors.gyro.std, ...
                'bias', sensors.gyro.bias, ...
                'tau_gyro_filter', sensors.gyro.tau, ...
                'last_filtered_omega', zeros(3,1),...
                'ts', sensors.gyro.Ts ...
            );
            
            % Initialize the Accelerometer model properties.
            obj.Accelerometer = struct(...
                'std', sensors.acc.std, ...
                'bias', sensors.acc.bias ...
            );
            
            % Initialize the Magnetometer model properties.
            obj.Magnetometer = struct(...
                'std', sensors.mag.std, ...
                'bias', sensors.mag.bias ...
            );

            % Store the inertial frame reference vectors.
            obj.g_I = sensors.acc.g_I;
            obj.m_I = sensors.mag.m_I;
        end
        
        function reading = getGyroscopeReading(obj, omega)
            %getGyroscopeReading Simulates a gyroscope measurement.
            %   Adds Gaussian noise and a constant bias to the true angular velocity.
            %
            %   Inputs:
            %       omega - 3x1 true angular velocity vector (rad/s).
            %   Outputs:
            %       reading - 3x1 simulated gyroscope measurement.

            noise = obj.Gyroscope.std .* randn(3, 1);
            reading = omega + noise + obj.Gyroscope.bias;
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
            %getAccelerometerReading Simulates an accelerometer measurement.
            %   Rotates the inertial gravity vector into the body frame, adds noise
            %   and bias, and then normalizes the result to a unit vector.
            %
            %   Inputs:
            %       R - 3x3 rotation matrix from body to inertial frame.
            %   Outputs:
            %       reading - 3x1 simulated (and normalized) accelerometer measurement.

            noise = obj.Accelerometer.std .* randn(3, 1);
            
            % Rotate gravity vector from inertial to body frame (R') and add imperfections.
            reading = R' * obj.g_I + noise + obj.Accelerometer.bias;
            
            % Normalize to represent a pure direction vector.
            reading = reading / norm(reading);
        end
        
        function reading = getMagnetometerReading(obj, R)
            %getMagnetometerReading Simulates a magnetometer measurement.
            %   Rotates the inertial magnetic vector into the body frame, adds noise
            %   and bias, and then normalizes the result.
            %
            %   Inputs:
            %       R - 3x3 rotation matrix from body to inertial frame.
            %   Outputs:
            %       reading - 3x1 simulated (and normalized) magnetometer measurement.

            noise = obj.Magnetometer.std .* randn(3, 1);
            
            % Rotate magnetic vector from inertial to body frame (R') and add imperfections.
            reading = R' * obj.m_I + noise + obj.Magnetometer.bias;
            
            % Normalize to represent a pure direction vector.
            reading = reading / norm(reading);
        end
    end
end