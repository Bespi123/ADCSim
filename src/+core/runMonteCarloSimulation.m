function mc = runMonteCarloSimulation(baseParameters, baseResults, options)
%RUNMONTECARLOSIMULATION Runs repeated ADCS simulations with randomized inputs.
%
%   mc = core.runMonteCarloSimulation(baseParameters, baseResults, options)
%
%   Inputs:
%       baseParameters - Parameter struct produced by simulation_Parameters.
%       baseResults    - Results struct produced by initialize_Results.
%       options        - Optional struct:
%           .numRuns              Number of simulations. Default: 100
%           .seed                 Random seed. Default: []
%           .showProgress          Print progress to command window. Default: true
%           .collectTrajectories   Store x, x_est, u for every run. Default: false
%           .perturbationFcn       Function handle for custom perturbations:
%                                 params = f(params, runIndex, options)
%           .perturbations         Struct with Gaussian 1-sigma perturbations:
%               .initialAttitudeStdRad  Small-angle attitude std [rad]. Default: 0
%               .initialRateStdRadSec   Body-rate std [rad/s]. Default: zeros(3,1)
%               .inertiaRelStd          Inertia relative std. Default: 0
%               .gyroBiasStdRadSec      Gyro bias std [rad/s]. Default: zeros(3,1)
%               .accBiasStd             Accelerometer bias std. Default: zeros(3,1)
%               .magBiasStd             Magnetometer bias std. Default: zeros(3,1)
%               .rwMotorRelStd          RW motor relative std. Default: 0
%               .disturbanceTorqueStdNm Additive disturbance std [Nm]. Default: 0
%
%   Output:
%       mc - Struct containing per-run metrics, optional trajectories, and
%            summary statistics.

    import core.simulation_rk4

    if nargin < 3 || isempty(options)
        options = struct();
    end
    options = applyDefaults(options);

    if ~isempty(options.seed)
        rng(options.seed);
    end

    numRuns = options.numRuns;
    metrics = preallocateMetrics(numRuns);
    runTemplate = struct( ...
        'seed', [], ...
        'error_flag', NaN, ...
        'message', '', ...
        'finalState', NaN(7,1), ...
        'finalEstimate', NaN(7,1), ...
        'indicators', []);
    runs = repmat(runTemplate, numRuns, 1);

    if options.collectTrajectories
        n = baseResults.adcSim.nSteps;
        trajectories.x = NaN(7, n, numRuns);
        trajectories.x_est = NaN(7, n, numRuns);
        trajectories.u = NaN(3, n, numRuns);
    else
        trajectories = struct();
    end

    for runIndex = 1:numRuns
        runSeed = randi(2^31 - 1);
        rng(runSeed);

        simParameters = perturbParameters(baseParameters, options, runIndex);
        simResults = perturbResults(baseResults, options);

        try
            [x, u, ~, x_est, indicators, ~, ~, error_flag] = simulation_rk4(simParameters, simResults, false);
            runs(runIndex).error_flag = error_flag;
            runs(runIndex).indicators = indicators;
            runs(runIndex).finalState = x(:, end);
            runs(runIndex).finalEstimate = padFinalEstimate(x_est);

            if error_flag == 0 && isstruct(indicators)
                metrics.success(runIndex) = true;
                metrics.RMSE(runIndex) = indicators.RMSE;
                metrics.finalEulerInt(runIndex) = indicators.eulerInt(end);
                metrics.finalAscct(runIndex) = indicators.ascct(end);
                metrics.ts(runIndex, :) = indicators.ts(:)';
            end

            if options.collectTrajectories
                trajectories.x(:, :, runIndex) = x;
                trajectories.u(:, :, runIndex) = u;
                trajectories.x_est(:, :, runIndex) = padHistory(x_est, size(trajectories.x_est, 2));
            end
        catch ME
            runs(runIndex).error_flag = 1;
            runs(runIndex).message = ME.message;
        end

        runs(runIndex).seed = runSeed;

        if options.showProgress
            fprintf('Monte Carlo run %d/%d complete. Success: %d\n', ...
                runIndex, numRuns, metrics.success(runIndex));
        end
    end

    mc = struct();
    mc.options = options;
    mc.runs = runs;
    mc.metrics = metrics;
    mc.summary = summarizeMetrics(metrics);
    mc.trajectories = trajectories;
end

function options = applyDefaults(options)
    options = setDefault(options, 'numRuns', 100);
    options = setDefault(options, 'seed', []);
    options = setDefault(options, 'showProgress', true);
    options = setDefault(options, 'collectTrajectories', false);
    options = setDefault(options, 'perturbationFcn', []);

    if ~isfield(options, 'perturbations') || isempty(options.perturbations)
        options.perturbations = struct();
    end

    p = options.perturbations;
    p = setDefault(p, 'initialAttitudeStdRad', 0);
    p = setDefault(p, 'initialRateStdRadSec', zeros(3,1));
    p = setDefault(p, 'inertiaRelStd', 0);
    p = setDefault(p, 'gyroBiasStdRadSec', zeros(3,1));
    p = setDefault(p, 'accBiasStd', zeros(3,1));
    p = setDefault(p, 'magBiasStd', zeros(3,1));
    p = setDefault(p, 'rwMotorRelStd', 0);
    p = setDefault(p, 'disturbanceTorqueStdNm', 0);
    options.perturbations = p;
end

function s = setDefault(s, fieldName, value)
    if ~isfield(s, fieldName) || isempty(s.(fieldName))
        s.(fieldName) = value;
    end
end

function metrics = preallocateMetrics(numRuns)
    metrics.success = false(numRuns, 1);
    metrics.RMSE = NaN(numRuns, 1);
    metrics.finalEulerInt = NaN(numRuns, 1);
    metrics.finalAscct = NaN(numRuns, 1);
    metrics.ts = NaN(numRuns, 3);
end

function simParameters = perturbParameters(baseParameters, options, runIndex)
    simParameters = baseParameters;
    p = options.perturbations;

    simParameters.initialValues.q0 = perturbQuaternion( ...
        simParameters.initialValues.q0, p.initialAttitudeStdRad);
    simParameters.initialValues.Wo = simParameters.initialValues.Wo + ...
        columnVector(p.initialRateStdRadSec) .* randn(3,1);

    if p.inertiaRelStd > 0
        rel = p.inertiaRelStd * randn(size(simParameters.spacecraft.I));
        rel = (rel + rel') / 2;
        simParameters.spacecraft.I = simParameters.spacecraft.I .* (1 + rel);
        simParameters.spacecraft.I = (simParameters.spacecraft.I + simParameters.spacecraft.I') / 2;
    end

    simParameters.sensors.gyro.bias = simParameters.sensors.gyro.bias + ...
        columnVector(p.gyroBiasStdRadSec) .* randn(3,1);
    simParameters.sensors.acc.bias = simParameters.sensors.acc.bias + ...
        columnVector(p.accBiasStd) .* randn(3,1);
    simParameters.sensors.mag.bias = simParameters.sensors.mag.bias + ...
        columnVector(p.magBiasStd) .* randn(3,1);

    if isfield(simParameters, 'sensors_soph')
        simParameters.sensors_soph.gyro.initial_bias = simParameters.sensors_soph.gyro.initial_bias + ...
            columnVector(p.gyroBiasStdRadSec) .* randn(3,1);
        simParameters.sensors_soph.acc.bias = simParameters.sensors_soph.acc.bias + ...
            columnVector(p.accBiasStd) .* randn(3,1);
        simParameters.sensors_soph.mag.bias = simParameters.sensors_soph.mag.bias + ...
            columnVector(p.magBiasStd) .* randn(3,1);
    end

    if p.rwMotorRelStd > 0
        motorFields = {'kt', 'Jrw', 'b', 'L', 'R', 'ke'};
        for k = 1:numel(motorFields)
            fieldName = motorFields{k};
            if isfield(simParameters.rw.motor, fieldName)
                scale = max(0, 1 + p.rwMotorRelStd * randn());
                simParameters.rw.motor.(fieldName) = simParameters.rw.motor.(fieldName) * scale;
            end
        end
    end

    if ~isempty(options.perturbationFcn)
        simParameters = options.perturbationFcn(simParameters, runIndex, options);
    end
end

function simResults = perturbResults(baseResults, options)
    simResults = baseResults;
    sigma = options.perturbations.disturbanceTorqueStdNm;
    if sigma <= 0
        return;
    end

    if isfield(simResults, 'disturbances') && isfield(simResults.disturbances, 'torque') && ...
            isfield(simResults.disturbances.torque, 'total')
        simResults.disturbances.torque.total = simResults.disturbances.torque.total + ...
            sigma * randn(size(simResults.disturbances.torque.total));
    end

    if isfield(simResults, 'disturbances') && isfield(simResults.disturbances, 'mask') && ...
            isfield(simResults.disturbances.mask, 'torque_dense')
        simResults.disturbances.mask.torque_dense = simResults.disturbances.mask.torque_dense + ...
            sigma * randn(size(simResults.disturbances.mask.torque_dense));
    end
end

function q = perturbQuaternion(q0, angleStd)
    q0 = columnVector(q0);
    if angleStd <= 0
        q = q0 / norm(q0);
        return;
    end

    dtheta = angleStd * randn(3,1);
    angle = norm(dtheta);
    if angle < eps
        dq = [1; 0; 0; 0];
    else
        axis = dtheta / angle;
        dq = [cos(angle/2); axis * sin(angle/2)];
    end

    q = quatProduct(dq, q0);
    q = q / norm(q);
    if q(1) < 0
        q = -q;
    end
end

function q = quatProduct(q1, q2)
    q = [q1(1)*q2(1) - q1(2:4)'*q2(2:4);
         q1(1)*q2(2:4) + q2(1)*q1(2:4) + cross(q1(2:4), q2(2:4))];
end

function v = columnVector(v)
    v = v(:);
    if isscalar(v)
        v = repmat(v, 3, 1);
    end
end

function xFinal = padFinalEstimate(x_est)
    xFinal = NaN(7,1);
    cols = find(all(~isnan(x_est), 1), 1, 'last');
    if ~isempty(cols)
        xFinal = x_est(:, cols);
    end
end

function padded = padHistory(x_est, n)
    padded = NaN(7, n);
    cols = min(size(x_est, 2), n);
    padded(:, 1:cols) = x_est(:, 1:cols);
    if cols > 0 && cols < n
        padded(:, cols+1:end) = repmat(padded(:, cols), 1, n-cols);
    end
end

function summary = summarizeMetrics(metrics)
    good = metrics.success;
    summary = struct();
    summary.numRuns = numel(good);
    summary.numSuccess = nnz(good);
    summary.successRate = nnz(good) / max(1, numel(good));
    summary.RMSE = describe(metrics.RMSE(good));
    summary.finalEulerInt = describe(metrics.finalEulerInt(good));
    summary.finalAscct = describe(metrics.finalAscct(good));
    summary.ts = describeRows(metrics.ts(good, :));
end

function stats = describe(values)
    values = values(~isnan(values));
    if isempty(values)
        stats = struct('mean', NaN, 'std', NaN, 'min', NaN, 'max', NaN, 'p05', NaN, 'p50', NaN, 'p95', NaN);
        return;
    end

    stats = struct();
    stats.mean = mean(values);
    stats.std = std(values);
    stats.min = min(values);
    stats.max = max(values);
    stats.p05 = percentile(values, 5);
    stats.p50 = percentile(values, 50);
    stats.p95 = percentile(values, 95);
end

function stats = describeRows(values)
    stats = struct('mean', NaN(1,3), 'std', NaN(1,3), 'min', NaN(1,3), ...
        'max', NaN(1,3), 'p05', NaN(1,3), 'p50', NaN(1,3), 'p95', NaN(1,3));
    for k = 1:3
        colStats = describe(values(:, k));
        stats.mean(k) = colStats.mean;
        stats.std(k) = colStats.std;
        stats.min(k) = colStats.min;
        stats.max(k) = colStats.max;
        stats.p05(k) = colStats.p05;
        stats.p50(k) = colStats.p50;
        stats.p95(k) = colStats.p95;
    end
end

function q = percentile(values, pct)
    values = sort(values(:));
    if isempty(values)
        q = NaN;
        return;
    end
    idx = 1 + (numel(values) - 1) * pct / 100;
    lo = floor(idx);
    hi = ceil(idx);
    if lo == hi
        q = values(lo);
    else
        q = values(lo) + (idx - lo) * (values(hi) - values(lo));
    end
end
