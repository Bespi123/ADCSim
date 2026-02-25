function optimal_gains = optimize_control_gains_ga(simParameters, simResults)
    % OPTIMIZE_CONTROL_GAINS_GA Finds ideal controller gains using a Genetic Algorithm.
    %
    % This function optimizes the control law (Feedback or Adaptive Boskovic)
    % by minimizing a fitness function that typically balances attitude error 
    % (ISE/ITAE) and control effort.
    
    import optimization.*

    % --- 1. Decision Logic for Controller Type ---
    % Selector 1: Linear Feedback Control [P and K matrices]
    % Selector 2: Boskovic Adaptive Control [delta, gamma, k0]
    
    if simParameters.controller.selector == 1
        nvars = 6;
        lb = zeros(1, 6);        % Minimum gain limit
        ub = 20 * ones(1, 6);    % Increased upper bound for higher exploration
        
        % Fetch current gains to seed the initial population
        P = simParameters.controller.feedback.Peye;
        K = simParameters.controller.feedback.Keye;
        initial_population = [diag(P)', diag(K)']; 
    else
        nvars = 3; 
        lb = [0, 0, 0]; 
        ub = [20, 5, 20]; 
        
        delta = simParameters.controller.boskController.delta;
        gamma = simParameters.controller.boskController.gamma;
        k0    = simParameters.controller.boskController.k0;
        initial_population = [delta, gamma, k0];
    end

    % --- 2. Genetic Algorithm Solver Options ---
    % Increase PopulationSize and MaxGenerations for PhD-level convergence.
    % Parallel computing is enabled to speed up the evaluation of individuals.
    
    options = optimoptions('ga', ...
        'PopulationSize', simParameters.ga.controller.popultaionSize, ...     % Increased for better diversity
        'MaxGenerations', simParameters.ga.controller.generations, ...        % More generations for fine-tuning
        'Display', 'iter', ...           
        'PlotFcn', {@gaplotbestf, @gaplotdistance, @gaplotrange, @gaplotscores}, ... 
        'UseParallel', true, ...         % Evaluates individuals in parallel (requires Parallel Toolbox)
        'InitialPopulationMatrix', initial_population, ...
        'FunctionTolerance', simParameters.ga.controller.tolerance, ...   % Stopping precision
        'Vectorized', 'off');           % RK4 simulation is not naturally vectorized

    % --- 3. Objective Function Wrapper ---
    % We pass the simulation environment handles to the fitness function.
    fitnessFcn = @(gains) objective_function_control(gains, simParameters, simResults);

    % --- 4. Optimization Execution ---
    fprintf('Starting Genetic Algorithm Optimization...\n');
    
    % Mute warnings regarding near-singularities (common in divergent control trials)
    warning('off', 'MATLAB:nearlySingularMatrix');
    
    [optimal_gains, min_cost] = ga(fitnessFcn, nvars, [], [], [], [], lb, ub, [], options);
    
    warning('on', 'MATLAB:nearlySingularMatrix');

    % --- 5. Reporting ---
    fprintf('Optimization Finished.\n');
    fprintf('Best Fitness Value (Min Cost): %.6f\n', min_cost);
    disp('Optimal Parameter Set:');
    disp(optimal_gains);
end