function cost = objective_function_control(gains, simParameters, simResults)
    % OBJECTIVE_FUNCTION_CONTROL Evaluates the fitness of controller gains.
    % 
    % This function computes a multi-objective scalar cost based on:
    % 1. Settling Time (Speed)
    % 2. ASCCT (Energy/Actuator effort)
    % 3. Euler Integral (Precision/Steady-state error)
    
    import core.*

    % --- 1. Parameter Update ---
    if simParameters.controller.selector == 1
        % Feedback Controller: Mapping gains to diagonal matrices
        simParameters.controller.feedback.Peye = diag(gains(1:3)); 
        simParameters.controller.feedback.Keye = diag(gains(4:6)); 
    else
        % Adaptive Boskovic Controller: Mapping gains to adaptive laws
        simParameters.controller.boskController.delta = gains(1); 
        simParameters.controller.boskController.gamma = gains(2); 
        simParameters.controller.boskController.k0    = gains(3); 
    end

    % --- 2. Numerical Simulation ---
    try
        % Run the RK4 simulation engine
        [~, ~, ~, ~, indicators, ~, ~, error_flag] = simulation_rk4(simParameters, simResults, false);
        
        % Check for numerical divergence or instability
        if error_flag == 1
            cost = 1e12; % High penalty for unstable systems
            return;      
        end
    catch
        % Catch unexpected solver crashes
        cost = 1e12; 
        return;
    end

    % --- 3. Weighted Cost Calculation with Normalization ---
    % Metrics should be normalized to be within similar orders of 
    % magnitude, or weights should be adjusted to compensate.
    
    % Weights (Adjust based on mission requirements: e.g., Power-constrained vs Speed-critical)
    w_ts    = 0.8;  % Weight for Settling Time
    w_ascct = 0.5;  % Weight for Control Effort (Energy)
    w_euler = 1.2;  % Weight for Integrated Error (Precision)
    
    % Data Extraction
    % Using mean settling time across all axes [Roll, Pitch, Yaw]
    ts_val = mean(indicators.ts);
    
    % If the system never settles (NaN in ts), apply a secondary penalty
    if any(isnan(indicators.ts))
        ts_val = simParameters.sim.tf * 2; % Penalty: twice the simulation time
    end
    val_ascct   = indicators.ascct(end);
    val_euler   = indicators.eulerInt(end);
    % --- 4. Cost Aggregation ---
    
    % Normalize values
    ref_ts  = simParameters.ga.controller.ref_ts; % stlemment time
    norm_ts = (ts_val    / ref_ts)^2;
    norm_ascct = (val_ascct / val_ascct)^2;
    norm_euler = (val_euler / val_euler)^2;

    % Using a weighted sum of performance indices
    cost_speed     = w_ts    * norm_ts;
    cost_energy    = w_ascct * norm_ascct;
    cost_precision = w_euler * norm_euler;
    
    cost = cost_speed + cost_energy + cost_precision;

    % --- 5. Sanity Check & Penalty ---
    % Ensure cost is a real scalar. If the simulation results in NaN states, 
    % it is treated as a failed individual.
    if isnan(cost) || isinf(cost)
        cost = 1e12;
    end
end