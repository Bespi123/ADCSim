function costs_vector = objective_function_ahrs(params, simParameters, simResults)
    import core.*
    
    % ekfObjectiveFunction - Evaluates the cost of a set of EKF parameters using a
    %                        Monte Carlo approach to handle stochastic noise.
    %
    % This function is intended to be called by optimize_ekf_ga.m
    % to determine the quality of a given set of EKF noise covariance
    % parameters. It runs multiple simulations to get a stable, averaged cost.
    %
    % Inputs:
    %   params        - Row vector containing the EKF parameters to be optimized.
    %   disturbances  - Structure with external disturbances for the simulation.
    %   simParameters - Structure with the base simulation parameters.
    %   time          - Structure with simulation time parameters.
    %
    % Output:
    %   costs_vector  - A row vector containing the final, averaged Root Mean
    %                   Square Error (RMSE) for each state. A lower value is better.

    % --- Monte Carlo Simulation Setup ---
    num_simulations = 4; % Number of simulations to run for averaging the cost.

    % Initialize a matrix to store the cost from each simulation run.
    % The number of columns should match the number of states being evaluated.
    % NOTE: This is currently set to 1. If 'indicators.RMSE' returns a
    % vector, you must adjust 'num_states' to match its length.
    num_states = 1;
    total_costs_matrix = zeros(num_simulations, num_states);

    % --- 1. Extract Parameters from Optimizer ---
    % The 'params' vector from the optimizer is used to update the EKF's
    % noise standard deviations for this evaluation.
    simParameters.ahrs.ekf.gyro_std = params(1:3)';
    simParameters.ahrs.ekf.acc_std  = params(4:6)';
    simParameters.ahrs.ekf.mag_std  = params(7:9)';
    if simParameters.sensors.star.enable == 1
        simParameters.ahrs.ekf.star_std  = params(10:12)';
    end
    simParameters.ahrs.enable = 1;
    simParameters.ahrs.enable_feedback = 0;     % Ensure EKF is not used as feedback during open-loop test.
    simParameters.ahrs.ekf.equalModel = 0; % Ensure custom EKF parameters are used.
    simParameters.ahrs.flag = 'EKF';

    % --- 2. Monte Carlo Loop ---
    % Run the simulation multiple times to average out the effects of random noise.
    for i = 1:num_simulations
        try
           
            % Run the main simulation function. It's assumed that this function
            % uses different random noise on each call.
            [~, ~, ~, ~, indicators, ~, ~, error_flag] = simulation_rk4(simParameters, simResults, false);
            
            if error_flag
                % If the simulation reports an error (e.g., divergence),
                % assign a very high penalty cost.
                total_costs_matrix(i, :) = 1e10;
            else
                % If the simulation is successful, get the cost.
                current_run_costs = indicators.RMSE;
                
                % Check for and penalize invalid NaN (Not-a-Number) results.
                current_run_costs(isnan(current_run_costs)) = 1e10;
                
                % Store the cost for this run.
                total_costs_matrix(i, :) = current_run_costs;
            end
        catch
            % If the simulation crashes for any unexpected reason, assign the penalty.
            total_costs_matrix(i, :) = 1e10;
        end
    end
    
    % --- 3. Calculate Final Averaged Cost ---
    % Average the costs from all simulation runs to get a stable metric.
    % The mean is taken along the first dimension (down the rows).
    avg_costs_per_state = mean(total_costs_matrix, 1);
    
    % The final cost is a row vector, scaled by 100.
    costs_vector = avg_costs_per_state * 100;
    
    % --- 4. Display Progress (Optional) ---
    % Print the averaged cost for this set of parameters to the console.
    fprintf('Average Cost per State (RMSE): %s\n', mat2str(costs_vector, 4));
end
