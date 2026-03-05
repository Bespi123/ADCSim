classdef ReactionWheelAssembly
    % REACTIONWHEELASSEMBLY Models a complete assembly of N reaction wheels.
    % * INCLUDES: Voltage saturation, Anti-windup, Speed and Current saturation.
    % * DIMENSIONS FIXED: States are [2 x N], Control signals are [N x 1].
    
    properties (SetAccess = public)
        States          % State matrix [2 x N_rw], where row 1 is speed, row 2 is current
        NumWheels       % The number of reaction wheels (N)
        MotorParams     % Struct with SCALAR motor parameters 
        voltage_vector  % Voltage input for motors [N x 1]
    end
    
    properties (Access = private)
        pid_Kp, pid_Ki, pid_Kd, pid_N
        pid_IntegratorState     % [N x 1] Column
        pid_DifferentiatorState % [N x 1] Column
        pid_LastError           % [N x 1] Column
    end
    
    methods
        function obj = ReactionWheelAssembly(rw)
            obj.NumWheels   = rw.number;
            obj.MotorParams = rw.motor;
            
            % States [2 x N], Voltages [N x 1]
            obj.States = zeros(2, obj.NumWheels);
            obj.voltage_vector = zeros(obj.NumWheels, 1);
            
            pid_params = rw.motor.pid;
            obj.pid_Kp = pid_params.kp;
            obj.pid_Ki = pid_params.ki;
            obj.pid_Kd = pid_params.kd;
            obj.pid_N  = pid_params.Nu;
            
            % Initialize PID states explicitly as COLUMN vectors [N x 1]
            obj.pid_IntegratorState     = zeros(obj.NumWheels, 1);
            obj.pid_DifferentiatorState = zeros(obj.NumWheels, 1);
            obj.pid_LastError           = zeros(obj.NumWheels, 1);
        end
        
        function [obj, torque_real_vector] = update(obj, speed_cmd_vector, dt)
            % 1. Get speed as a column vector [N x 1]
            current_speed_vector = obj.States(1, :)'; 
            
            % 2. Calculate error as a strict column vector [N x 1]
            speed_error_vector = speed_cmd_vector(:) - current_speed_vector;
            
            % 3. PID logic (all internals are now N x 1)
            [obj, obj.voltage_vector] = obj.calculatePID(speed_error_vector, dt);
            
            % 4. Integrate dynamics (States remain 2 x N)
            currentStateMatrix = obj.States; 
            
            g1 = dt * obj.brushlessModelMatrix(currentStateMatrix, obj.voltage_vector);
            g2 = dt * obj.brushlessModelMatrix(currentStateMatrix + 0.5 * g1, obj.voltage_vector);
            g3 = dt * obj.brushlessModelMatrix(currentStateMatrix + 0.5 * g2, obj.voltage_vector);
            g4 = dt * obj.brushlessModelMatrix(currentStateMatrix + g3, obj.voltage_vector);
            
            newStateMatrix = currentStateMatrix + (1/6) * (g1 + 2*g2 + 2*g3 + g4);
            
            % --- STATE SATURATION ---
            w_max = obj.MotorParams.max_speed; 
            newStateMatrix(1, :) = max(min(newStateMatrix(1, :), w_max), -w_max);
            
            i_max = obj.MotorParams.max_current; 
            newStateMatrix(2, :) = max(min(newStateMatrix(2, :), i_max), -i_max);
            
            % 5. Real torque [N x 1]
            angular_acceleration_row = (newStateMatrix(1, :) - obj.States(1, :)) / dt;
            torque_real_vector = (obj.MotorParams.Jrw * angular_acceleration_row)'; % Transpose to [N x 1]
            
            obj.States = newStateMatrix;
        end
    end
    
    methods (Access = private)
        function [obj, u_vector] = calculatePID(obj, error_vector, dt)
            % error_vector is [N x 1]
            
            error_dot_vector = (error_vector - obj.pid_LastError) / dt;
            
            % Runge-Kutta filter for derivative (all [N x 1])
            k1 = obj.pid_N * (error_dot_vector - obj.pid_DifferentiatorState);
            k2 = obj.pid_N * (error_dot_vector - (obj.pid_DifferentiatorState + k1*dt/2));
            k3 = obj.pid_N * (error_dot_vector - (obj.pid_DifferentiatorState + k2*dt/2));
            k4 = obj.pid_N * (error_dot_vector - (obj.pid_DifferentiatorState + k3*dt));
            obj.pid_DifferentiatorState = obj.pid_DifferentiatorState + (k1 + 2*k2 + 2*k3 + k4)*dt/6;
            
            % Tentative integral [N x 1]
            tentative_integral = obj.pid_IntegratorState + error_vector * dt;
            
            % Unbounded control [N x 1]
            u_unbounded = obj.pid_Kp * error_vector + ...
                          obj.pid_Ki * tentative_integral + ...
                          obj.pid_Kd * obj.pid_DifferentiatorState;
            
            % VOLTAGE SATURATION
            v_max = obj.MotorParams.max_voltage;
            u_vector = max(min(u_unbounded, v_max), -v_max);
            
            % ANTI-WINDUP (Safe element-wise comparison for [N x 1])
            is_saturated = (u_unbounded ~= u_vector);
            is_same_sign = (sign(error_vector) == sign(u_unbounded));
            
            windup_flag = is_saturated & is_same_sign; % [N x 1] boolean
            
            % Conditional integration
            obj.pid_IntegratorState = obj.pid_IntegratorState + error_vector .* (~windup_flag) * dt;
            
            obj.pid_LastError = error_vector;
        end
        
        function dx = brushlessModelMatrix(obj, x, u)
            % x is [2 x N], u is [N x 1]
            w = x(1, :); % [1 x N]
            i = x(2, :); % [1 x N]
            
            % Force 'u' to be a row vector [1 x N] for the physics equation
            u_row = u(:)'; 
            
            p = obj.MotorParams;
            
            w_dot = (1/p.Jrw) * (p.kt * i - p.b * w - p.c * sign(w));
            i_dot = (1/p.L)   * (u_row - p.R * i - p.ke * w); 
            
            dx = [w_dot; i_dot]; % Returns [2 x N]
        end
    end
end