function res = initialize_Results(params)
% INITIALIZE_RESULTS Creates a pre-allocated structure to store simulation outputs.
%
%   Pre-allocating memory with 'zeros' instead of empty brackets [] 
%   significantly improves performance by preventing MATLAB from 
%   re-allocating memory during the numerical integration (RK4).
%
%   INPUT:  params - Configuration structure from simulation_Parameters.
%   OUTPUT: res    - Structure with pre-allocated matrices.

    % --- Time & Dimension Calculation ---
    % Generating the time vector to determine the exact number of iterations (n).
    t0 = params.sim.t0;
    tf = params.sim.tf;
    dt = params.sim.step;
    t_vector = t0:dt:tf;
    n = numel(t_vector); 

    % --- Global Metadata ---
    res.adcSim.t = t_vector;
    res.adcSim.nSteps = n;
    
    % --- Orbital Data ---
    % Placeholder for Orekit high-fidelity orbital propagation results.
    res.orbital_propagation = [];

    % --- Simulation Results (Core States) ---
    % x: [q0; q1; q2; q3; w_x; w_y; w_z] -> 7-state vector (Quaternion + Body Rates)
    res.adcSim.x          = zeros(7, n); 
    res.adcSim.x_est      = zeros(7, n); % AHRS/EKF state estimation
    res.adcSim.u          = zeros(3, n); % Control command (Ideal torque request)
    res.adcSim.error_flag = 0;           % Simulation stability flag (0 = Success)
    
    % --- Actuator States (Reaction Wheels) ---
    num_rw = params.rw.number; 
    res.adcSim.torque                = zeros(num_rw, n); % Commanded torque per wheel
    res.adcSim.actuators.torque_real = zeros(num_rw, n); % Actual applied torque (with friction/saturation)
    res.adcSim.actuators.w_cmd       = zeros(num_rw, n); % Commanded wheel speeds [rad/s]
    
    % x_rw: Typically stores [speed; current] for each wheel -> (num_rw * 2)
    res.adcSim.actuators.x_rw        = zeros(num_rw * 2, n); 
    res.adcSim.actuators.u_rw        = zeros(num_rw, n);     % Motor input voltage [V]
    
    % --- Sensor Data & Estimation ---
    % star_number represents the dimension of star tracker measurements
    star_dim = params.sensors.star.numberOfStars; 
    
    % Raw measurements (e.g., IMU, Magnetometer, Sun Sensors, Star Tracker)
    res.adcSim.sensors.meas           = zeros(9 + star_dim, n);
    % Estimated sensor biases or filtered states
    res.adcSim.sensors.est            = zeros(6 + star_dim, n);
    res.adcSim.sensors.omega_filtered = zeros(3, n); % Filtered body angular rates
    
    % --- Performance Indices ---
    % Settling time (ts) for [Roll; Pitch; Yaw]. Initialized as NaN.
    res.adcSim.index.ts       = NaN(3, 1); 
    res.adcSim.index.eulerInt = zeros(1, n); % Cumulative integral of attitude error
    res.adcSim.index.ascct    = zeros(1, n); % Accumulated Spacecraft Control Torque (Effort)
    res.adcSim.index.o        = zeros(1, n); % Computational cost per step (Execution time)
    res.adcSim.index.RMSE     = 0;           % Root Mean Square Error of attitude estimation

    % --- Environmental Torques (Body Frame) ---
    % Pre-allocated for individual torque disturbance components [Nm]
    res.disturbances.torque.total          = zeros(3, n);
    res.disturbances.torque.gravity_torque = zeros(3, n);
    res.disturbances.torque.drag_torque    = zeros(3, n);
    res.disturbances.torque.solar_torque   = zeros(3, n);
    res.disturbances.torque.t              = t_vector;
   
end