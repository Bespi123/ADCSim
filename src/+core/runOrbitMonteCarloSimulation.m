function orbitMc = runOrbitMonteCarloSimulation(baseParameters, options)
%RUNORBITMONTECARLOSIMULATION Runs Orekit orbit propagation with randomized inputs.
%
%   orbitMc = core.runOrbitMonteCarloSimulation(baseParameters, options)
%
%   The output histories can be passed to core.runTorqueMonteCarloFromOrbit
%   to generate disturbance torque scenarios for attitude simulations.

    import core.runOrekitSimulation

    if nargin < 2 || isempty(options)
        options = struct();
    end
    options = applyDefaults(options);

    if ~isempty(options.seed)
        rng(options.seed);
    end

    numRuns = options.numRuns;
    runTemplate = struct( ...
        'seed', [], ...
        'success', false, ...
        'message', '', ...
        'parameters', [], ...
        'history', []);
    runs = repmat(runTemplate, numRuns, 1);

    for runIndex = 1:numRuns
        runSeed = randi(2^31 - 1);
        rng(runSeed);

        params = perturbOrbitParameters(baseParameters, options, runIndex);
        runs(runIndex).seed = runSeed;
        runs(runIndex).parameters = params;

        try
            runs(runIndex).history = runOrekitSimulation(params);
            runs(runIndex).success = true;
        catch ME
            runs(runIndex).message = ME.message;
        end

        if options.showProgress
            fprintf('Orbit Monte Carlo run %d/%d complete. Success: %d\n', ...
                runIndex, numRuns, runs(runIndex).success);
        end
    end

    orbitMc = struct();
    orbitMc.options = options;
    orbitMc.runs = runs;
    orbitMc.summary = summarizeRuns(runs);
end

function options = applyDefaults(options)
    options = setDefault(options, 'numRuns', 100);
    options = setDefault(options, 'seed', []);
    options = setDefault(options, 'showProgress', true);
    options = setDefault(options, 'perturbationFcn', []);

    if ~isfield(options, 'perturbations') || isempty(options.perturbations)
        options.perturbations = struct();
    end

    p = options.perturbations;
    p = setDefault(p, 'initialDateStdSec', 0);
    p = setDefault(p, 'semiMajorAxisStdM', 0);
    p = setDefault(p, 'eccentricityStd', 0);
    p = setDefault(p, 'inclinationStdRad', 0);
    p = setDefault(p, 'raanStdRad', 0);
    p = setDefault(p, 'argPerigeeStdRad', 0);
    p = setDefault(p, 'trueAnomalyStdRad', 0);
    p = setDefault(p, 'massRelStd', 0);
    p = setDefault(p, 'dragAreaRelStd', 0);
    p = setDefault(p, 'dragCoeffRelStd', 0);
    p = setDefault(p, 'srpAreaRelStd', 0);
    p = setDefault(p, 'srpCoeffRelStd', 0);
    options.perturbations = p;
end

function s = setDefault(s, fieldName, value)
    if ~isfield(s, fieldName) || isempty(s.(fieldName))
        s.(fieldName) = value;
    end
end

function params = perturbOrbitParameters(baseParameters, options, runIndex)
    params = baseParameters;
    p = options.perturbations;

    if p.initialDateStdSec > 0
        params.initialDate = params.initialDate + seconds(p.initialDateStdSec * randn());
    end

    params.initialOrbit.a = max(1, params.initialOrbit.a + p.semiMajorAxisStdM * randn());
    params.initialOrbit.e = min(max(params.initialOrbit.e + p.eccentricityStd * randn(), 0), 0.99);
    params.initialOrbit.i = params.initialOrbit.i + p.inclinationStdRad * randn();
    params.initialOrbit.raan = wrapToTwoPiLocal(params.initialOrbit.raan + p.raanStdRad * randn());
    params.initialOrbit.pa = wrapToTwoPiLocal(params.initialOrbit.pa + p.argPerigeeStdRad * randn());
    params.initialOrbit.ta = wrapToTwoPiLocal(params.initialOrbit.ta + p.trueAnomalyStdRad * randn());

    params.spacecraft.mass = positiveScale(params.spacecraft.mass, p.massRelStd);
    params.spacecraft.dragArea = positiveScale(params.spacecraft.dragArea, p.dragAreaRelStd);
    params.spacecraft.dragCoeff = positiveScale(params.spacecraft.dragCoeff, p.dragCoeffRelStd);
    params.spacecraft.srpArea = positiveScale(params.spacecraft.srpArea, p.srpAreaRelStd);
    params.spacecraft.srpCoeff = positiveScale(params.spacecraft.srpCoeff, p.srpCoeffRelStd);

    if ~isempty(options.perturbationFcn)
        params = options.perturbationFcn(params, runIndex, options);
    end
end

function value = positiveScale(value, relStd)
    if relStd > 0
        value = value * max(0, 1 + relStd * randn());
    end
end

function angle = wrapToTwoPiLocal(angle)
    angle = mod(angle, 2*pi);
end

function summary = summarizeRuns(runs)
    success = [runs.success]';
    summary = struct();
    summary.numRuns = numel(runs);
    summary.numSuccess = nnz(success);
    summary.successRate = nnz(success) / max(1, numel(runs));
end
