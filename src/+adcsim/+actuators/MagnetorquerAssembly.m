classdef MagnetorquerAssembly
    % MAGNETORQUERASSEMBLY Models a complete assembly of N magnetorquers.
    % * This class manages the state, control, and dynamics for a set of
    % * magnetic torquers (coils/rods) in a fully vectorized manner.
    % * INCLUDES: Voltage saturation, Dipole saturation, and cross-product torque generation.
    
    properties (SetAccess = public)
        NumTorquers     % The number of magnetorquers (usually 3 for X, Y, Z)
        Params          % Struct with SCALAR parameters (max_dipole, resistance, max_voltage, etc.)
        Configuration   % Matrix [3 x N] mapping torquer axes to body axes
        DipoleStates    % Current generated dipole moment [N x 1] (A*m^2)
        VoltageVector   % Current voltage applied to each torquer [N x 1] (V)
    end
    
    methods
        function obj = MagnetorquerAssembly(mtq_params)
            % Constructor: Creates the assembly for N magnetorquers.
            %
            % Inputs:
            %   mtq_params - Struct containing:
            %       .number        (N)
            %       .W             [3 x N] Configuration matrix (e.g., eye(3))
            %       .max_dipole    Maximum magnetic dipole moment [A*m^2]
            %       .max_voltage   Maximum voltage from the power bus [V]
            %       .resistance    Internal resistance of the coil [Ohms]
            %       .area          Cross-sectional area [m^2] (Optional for detailed model)
            %       .turns         Number of wire turns (Optional for detailed model)
            
            obj.NumTorquers = mtq_params.number;
            obj.Configuration = mtq_params.W;
            obj.Params = mtq_params;
            
            % Initialize states as column vectors [N x 1]
            obj.DipoleStates = zeros(obj.NumTorquers, 1);
            obj.VoltageVector = zeros(obj.NumTorquers, 1);
        end
        
        function [obj, torque_body, current_vector] = update(obj, dipole_cmd_vector, b_field_body)
            % update: Simulates the magnetorquer dynamics and generates torque.
            %
            % Inputs:
            %   dipole_cmd_vector - Commanded magnetic dipole for each coil [N x 1] (A*m^2)
            %   b_field_body      - Local magnetic field measured/estimated in Body frame [3 x 1] (Teslas)
            %
            % Outputs:
            %   obj            - Updated object
            %   torque_body    - Real mechanical torque generated on the satellite chasis [3 x 1] (Nm)
            %   current_vector - Electrical current consumed by each coil [N x 1] (Amps)
            
            % 1. Command Saturation (Dipole Limit)
            % The commanded dipole cannot exceed the physical limits of the coil core.
            m_max = obj.Params.max_dipole;
            dipole_cmd_vector = max(min(dipole_cmd_vector(:), m_max), -m_max);
            
            % 2. Electrical Model (Voltage and Current Mapping)
            % In a simple DC model (ignoring fast inductance transients):
            % Dipole = N * I * A  -->  I = Dipole / (N*A)
            % If N*A is not explicitly provided, we can map it directly if we know 
            % the max dipole corresponds to max current/voltage.
            % Let's assume a direct proportional conversion factor: K_m = m_max / I_max
            % I = V / R
            
            % Calculate required current to achieve the commanded dipole
            % (Assuming linear core for simplicity, no hysteresis in this basic model)
            if isfield(obj.Params, 'turns') && isfield(obj.Params, 'area')
                current_required = dipole_cmd_vector / (obj.Params.turns * obj.Params.area);
            else
                % Simplified mapping if geometric parameters are missing
                % We assume a conversion constant: m = K_coil * I
                K_coil = obj.Params.max_dipole / (obj.Params.max_voltage / obj.Params.resistance);
                current_required = dipole_cmd_vector / K_coil;
            end
            
            % Calculate required voltage
            voltage_required = current_required * obj.Params.resistance;
            
            % 3. Voltage Saturation (Power Bus Limit)
            v_max = obj.Params.max_voltage;
            obj.VoltageVector = max(min(voltage_required, v_max), -v_max);
            
            % 4. Realized Current and Dipole
            current_vector = obj.VoltageVector / obj.Params.resistance;
            
            if isfield(obj.Params, 'turns') && isfield(obj.Params, 'area')
                obj.DipoleStates = obj.Params.turns * obj.Params.area * current_vector;
            else
                obj.DipoleStates = K_coil * current_vector;
            end
            
            % 5. Total Dipole Vector in Body Frame [3 x 1]
            % Map the individual coil dipoles to the satellite's axes using the configuration matrix W
            total_dipole_body = obj.Configuration * obj.DipoleStates;
            
            % 6. Physics: Torque Generation (Cross Product Law)
            % Torque = m_body x B_body
            torque_body = cross(total_dipole_body, b_field_body(:));
        end
    end
end