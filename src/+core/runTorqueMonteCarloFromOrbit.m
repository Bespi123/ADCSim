function torqueMc = runTorqueMonteCarloFromOrbit(orbitMc, ~, baseResults, options)
%RUNTORQUEMONTECARLOFROMORBIT Converts orbit Monte Carlo histories to torque cases.
%
%   torqueMc = core.runTorqueMonteCarloFromOrbit(orbitMc, params, results, options)
%
%   Each successful Orekit history is converted into gravity-gradient, drag,
%   and SRP torque histories, then resampled onto results.adcSim.t. The output
%   includes ready-to-use simParameters/simResults pairs for simulation_rk4.

    if nargin < 4 || isempty(options)
        options = struct();
    end
    options = applyDefaults(options);

    if ~isempty(options.seed)
        rng(options.seed);
    end

    sourceRuns = orbitMc.runs;
    numRuns = numel(sourceRuns);
    runTemplate = struct( ...
        'seed', [], ...
        'success', false, ...
        'message', '', ...
        'orbitRunIndex', NaN, ...
        'simParameters', [], ...
        'simResults', [], ...
        'torque', [], ...
        'magneticField', []);
    runs = repmat(runTemplate, numRuns, 1);

    for runIndex = 1:numRuns
        runSeed = randi(2^31 - 1);
        rng(runSeed);

        runs(runIndex).seed = runSeed;
        runs(runIndex).orbitRunIndex = runIndex;

        if ~sourceRuns(runIndex).success
            runs(runIndex).message = 'Skipped because the source orbit run failed.';
            continue;
        end

        try
            simParameters = sourceRuns(runIndex).parameters;
            simParameters.disturbances.model_propagator.enable = 0;

            simResults = baseResults;
            torqueData = buildTorqueData(sourceRuns(runIndex).history, simParameters, options);
            torqueData = resampleTorqueData(torqueData, simResults.adcSim.t);

            simResults.orbital_propagation = sourceRuns(runIndex).history;
            simResults.disturbances.torque.gravity_torque = torqueData.gravity;
            simResults.disturbances.torque.drag_torque = torqueData.drag;
            simResults.disturbances.torque.solar_torque = torqueData.srp;
            simResults.disturbances.torque.total = torqueData.total;
            simResults.disturbances.torque.magnetic_field = torqueData.magneticField;
            simResults.disturbances.mask.t_dense = simResults.adcSim.t;
            simResults.disturbances.mask.torque_dense = torqueData.total;
            simResults.disturbances.mask.magnetic_field = torqueData.magneticField;

            runs(runIndex).success = true;
            runs(runIndex).simParameters = simParameters;
            runs(runIndex).simResults = simResults;
            runs(runIndex).torque = torqueData;
            runs(runIndex).magneticField = torqueData.magneticField;
        catch ME
            runs(runIndex).message = ME.message;
        end

        if options.showProgress
            fprintf('Torque Monte Carlo case %d/%d complete. Success: %d\n', ...
                runIndex, numRuns, runs(runIndex).success);
        end
    end

    torqueMc = struct();
    torqueMc.options = options;
    torqueMc.runs = runs;
    torqueMc.summary = summarizeRuns(runs);
end

function options = applyDefaults(options)
    options = setDefault(options, 'seed', []);
    options = setDefault(options, 'showProgress', true);
    options = setDefault(options, 'perturbationFcn', []);

    if ~isfield(options, 'perturbations') || isempty(options.perturbations)
        options.perturbations = struct();
    end

    p = options.perturbations;
    p = setDefault(p, 'inertiaRelStd', 0);
    p = setDefault(p, 'massRelStd', 0);
    p = setDefault(p, 'dragLeverArmStdM', zeros(3,1));
    p = setDefault(p, 'srpLeverArmStdM', zeros(3,1));
    p = setDefault(p, 'torqueWhiteNoiseStdNm', zeros(3,1));
    p = setDefault(p, 'magneticFieldRelStd', 0);
    options.perturbations = p;
end

function s = setDefault(s, fieldName, value)
    if ~isfield(s, fieldName) || isempty(s.(fieldName))
        s.(fieldName) = value;
    end
end

function torqueData = buildTorqueData(history, simParameters, options)
    simParameters = perturbTorqueParameters(simParameters, options);

    disturbanceModel = adcsim.environment.torque_disturbances(simParameters.spacecraft, history);
    gravity = disturbanceModel.gravity_torque_model(history(:, 26:29), simParameters.initialOrbit);
    drag = disturbanceModel.drag_torque_model();
    srp = disturbanceModel.solar_torque_model();
    magneticField = history(:, 23:25);

    total = gravity + drag + srp;
    noiseStd = columnVector(options.perturbations.torqueWhiteNoiseStdNm)';
    if any(noiseStd > 0)
        total = total + randn(size(total)) .* noiseStd;
    end

    if options.perturbations.magneticFieldRelStd > 0
        magneticField = magneticField .* (1 + options.perturbations.magneticFieldRelStd * randn(size(magneticField)));
    end

    if ~isempty(options.perturbationFcn)
        [total, gravity, drag, srp, magneticField] = options.perturbationFcn( ...
            total, gravity, drag, srp, magneticField, simParameters, options);
    end

    torqueData = struct();
    torqueData.time = history(:, 1)';
    torqueData.gravity = gravity';
    torqueData.drag = drag';
    torqueData.srp = srp';
    torqueData.total = total';
    torqueData.magneticField = magneticField';
end

function simParameters = perturbTorqueParameters(simParameters, options)
    p = options.perturbations;

    if p.inertiaRelStd > 0
        rel = p.inertiaRelStd * randn(size(simParameters.spacecraft.I));
        rel = (rel + rel') / 2;
        simParameters.spacecraft.I = simParameters.spacecraft.I .* (1 + rel);
        simParameters.spacecraft.I = (simParameters.spacecraft.I + simParameters.spacecraft.I') / 2;
    end

    simParameters.spacecraft.mass = positiveScale(simParameters.spacecraft.mass, p.massRelStd);
    simParameters.spacecraft.r_cp_cg = simParameters.spacecraft.r_cp_cg + columnVector(p.dragLeverArmStdM)' .* randn(1,3);
    simParameters.spacecraft.r_sun_cg = simParameters.spacecraft.r_sun_cg + columnVector(p.srpLeverArmStdM)' .* randn(1,3);
end

function torqueData = resampleTorqueData(torqueData, targetTime)
    targetTime = targetTime(:)';
    sourceTime = torqueData.time(:);

    torqueData.time = targetTime;
    torqueData.gravity = interpRows(sourceTime, torqueData.gravity', targetTime)';
    torqueData.drag = interpRows(sourceTime, torqueData.drag', targetTime)';
    torqueData.srp = interpRows(sourceTime, torqueData.srp', targetTime)';
    torqueData.total = interpRows(sourceTime, torqueData.total', targetTime)';
    torqueData.magneticField = interpRows(sourceTime, torqueData.magneticField', targetTime)';
end

function values = interpRows(sourceTime, sourceValues, targetTime)
    [sourceTime, uniqueIdx] = unique(sourceTime, 'stable');
    sourceValues = sourceValues(uniqueIdx, :);
    values = interp1(sourceTime, sourceValues, targetTime(:), 'linear', 'extrap');
end

function value = positiveScale(value, relStd)
    if relStd > 0
        value = value * max(0, 1 + relStd * randn());
    end
end

function v = columnVector(v)
    v = v(:);
    if isscalar(v)
        v = repmat(v, 3, 1);
    end
end

function summary = summarizeRuns(runs)
    success = [runs.success]';
    summary = struct();
    summary.numRuns = numel(runs);
    summary.numSuccess = nnz(success);
    summary.successRate = nnz(success) / max(1, numel(runs));
end
