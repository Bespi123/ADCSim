classdef torque_disturbances
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        satellite %satellite object
        
    end

    methods
        function obj = gravity_torque_model(state, sat, orbit, eart)
        % DIST_GRAVITYGRAD Calcula el torque por gradiente de gravedad.
        % Útil para satélites con inercias asimétricas como un CubeSat 12U.
        
            % 1. Extraer cuaternión (asumiendo formato [q_w, q_x, q_y, q_z])
            q = state(7:10);
            q = q / norm(q); % Normalización de seguridad
        
            % 2. Parámetros orbitales
            % Usar mu desde la estructura earth si está disponible, sino el estándar
            if isfield(earth, 'mu')
                mu = earth.mu;
            else
                mu = 3.986004418e14; % [m^3/s^2]
            end
            
            rc = earth.Radius + orbit.altitude;
            omega_o_sq = mu / rc^3; % Velocidad orbital al cuadrado (w_o^2)
        
            % 3. Dirección del Nadir en el marco del cuerpo (Body Frame)
            % Representa la dirección hacia el centro de la Tierra
            c3 = [2*(q(2)*q(4) - q(3)*q(1));
                  2*(q(3)*q(4) + q(2)*q(1));
                  1 - 2*(q(2)^2 + q(3)^2)];
        
            % 4. Cálculo del Torque de Gradiente de Gravedad
            % Tgg = 3 * w_o^2 * (n x (I * n))
            Tgg = 3 * omega_o_sq * cross(c3, sat.Is * c3);
        end

        function outputArg = drag_torque_model(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end

        function outputArg = solar_torque_model(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end