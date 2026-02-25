function [optimal_params, min_cost] = optimize_ahrs_ga(simParameters, simResults)
    import optimization.*
    % OPTIMIZE_AHRS_GA Tunes EKF noise covariance parameters using Genetic Algorithms.
    %
    % The goal is to find the optimal Standard Deviations (sigma) for the 
    % sensor noise matrices (R) and process noise (Q) that minimize estimation error.
    
    % --- 1. Objective Function Setup ---
    % We pass the handle so 'ga' can run the simulation internally
    fitnessFcn = @(params) objective_function_ahrs(params, simParameters, simResults);
    
    % --- 2. Variable Definition ---
    % Check if Star Tracker is active to adjust the search vector length
    useStar = simParameters.sensors.star.enable == 1;
    
    if useStar
        nvars = 12; % [Gyro(3), Acc(3), Mag(3), Star(3)]
    else
        nvars = 9;  % [Gyro(3), Acc(3), Mag(3)]
    end
    
    % --- 3. Search Space (Bounds) ---
    % EKF tuning is sensitive. We search for Standard Deviations (std).
    % Lower Bound: Close to 0 but positive (numerical stability).
    % Upper Bound: 0.1 is usually huge for sensors like Gyros (rad/s), 
    % but allows the GA to reject bad sensors if needed.
    lb = 1e-8 * ones(1, nvars); 
    ub = 0.5 * ones(1, nvars);  
    
    % --- 4. Initial Population (Seeding) ---
    % We use the current manual tuning as the starting point.
    % Using (:)' ensures we get a Row Vector regardless of input shape.
    p = simParameters.ahrs.ekf;
    
    gyro0 = p.gyro_std(:)';
    acc0  = p.acc_std(:)';
    mag0  = p.mag_std(:)';
    
    if useStar
        % Verify exact path in your struct for star std
        % Assuming simParameters.ahrs.ekf.star.std based on previous context
        try
            star0 = p.star_std(:)'; 
        catch
            % Fallback if structure path is different
            star0 = [0.001, 0.001, 0.001]; 
        end
        initial_pop = [gyro0, acc0, mag0, star0];
    else
        initial_pop = [gyro0, acc0, mag0];
    end
    
    % --- 5. GA Solver Options ---
    % Critical for convergence. 
    options = optimoptions('ga', ...
        'PopulationSize', 50, ...        % Diversity: 50 individuals is standard for 9-12 vars.
        'MaxGenerations', 30, ...        % Iterations: Give it time to evolve.
        'UseParallel', true, ...         % Speed: Use all CPU cores.
        'Display', 'iter', ...           % Feedback: Show text log.
        'InitialPopulationMatrix', initial_pop, ...
        'FunctionTolerance', 1e-6, ...   % Precision: Stop if improvement is negligible.
        'PlotFcn', {@gaplotbestf, @gaplotdistance, @gaplotrange, @gaplotscores}); 
        
    % --- 6. Execution ---
    disp('Starting AHRS (EKF) Parameter Optimization...');
    disp(['Variables to tune: ', num2str(nvars)]);
    
    [optimal_params, min_cost] = ga(fitnessFcn, nvars, [], [], [], [], lb, ub, [], options);
    
    % --- 7. Reporting ---
    fprintf('\nOptimization Finished.\n');
    fprintf('Best RMSE Cost found: %.6f\n', min_cost);
    
end