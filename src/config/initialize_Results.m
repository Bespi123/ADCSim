function res = initialize_Results(params)
% INITIALIZE_RESULTS Creates a clean structure to store simulation outputs.
%
%   INPUT:  params - The configuration structure from simulation_Parameters.
%   OUTPUT: res    - A structure with pre-allocated matrices (zeros).

    N = params.sim.nSteps; % Get the total number of steps

    % --- Time Data ---
    res.adcSim.t = zeros(1, N); % Vector time [maneuver]
    res.adcSim.nSteps = 0;      % Vector time steps

    % --- Attitude Results (Quaternions and Angular Rates) ---
    res.attitude.q_true = zeros(N, 4); % Actual orientation
    res.attitude.q_est  = zeros(N, 4); % AHRS estimated orientation
    res.dynamics.omega  = zeros(N, 3); % Angular velocity

    % --- Environmental Torques (Body Frame) ---
    res.disturbances.torque = [];

    % --- Actuator States (Reaction Wheels) ---
    res.actuators.rw_speed   = zeros(N, params.rw.number);
    res.actuators.rw_current = zeros(N, params.rw.number);
    
    % --- Sensor Raw Data (For debugging AHRS) ---
    res.sensors.gyro_raw = zeros(N, 3);
    res.sensors.mag_raw  = zeros(N, 3);
end