%% Orbit-to-torque Monte Carlo example for ADCSim
% This runs Orekit multiple times, converts each orbit history into torque
% disturbance histories, then runs one attitude simulation with the first case.

clear;
clc;
close all;

project_root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(project_root, 'src')));
addpath(genpath(fullfile(project_root, 'lib')));

% Orekit must be configured before running orbit propagation.
run(fullfile(project_root, 'src', 'config', 'import_orekit.m'));

params = simulation_Parameters([]);
results = initialize_Results(params);

orbitOpts = struct();
orbitOpts.numRuns = 5;
orbitOpts.seed = 42;
orbitOpts.showProgress = true;
orbitOpts.perturbations.semiMajorAxisStdM = 100;
orbitOpts.perturbations.inclinationStdRad = deg2rad(0.01);
orbitOpts.perturbations.trueAnomalyStdRad = deg2rad(0.1);
orbitOpts.perturbations.massRelStd = 0.02;
orbitOpts.perturbations.dragAreaRelStd = 0.05;
orbitOpts.perturbations.dragCoeffRelStd = 0.05;
orbitOpts.perturbations.srpAreaRelStd = 0.05;
orbitOpts.perturbations.srpCoeffRelStd = 0.05;

orbitMc = core.runOrbitMonteCarloSimulation(params, orbitOpts);
disp('Orbit Monte Carlo summary:');
disp(orbitMc.summary);

torqueOpts = struct();
torqueOpts.seed = 7;
torqueOpts.showProgress = true;
torqueOpts.perturbations.inertiaRelStd = 0.02;
torqueOpts.perturbations.dragLeverArmStdM = [0.002; 0.002; 0.002];
torqueOpts.perturbations.srpLeverArmStdM = [0.002; 0.002; 0.002];
torqueOpts.perturbations.torqueWhiteNoiseStdNm = [1e-7; 1e-7; 1e-7];
torqueOpts.perturbations.magneticFieldRelStd = 0.01;

torqueMc = core.runTorqueMonteCarloFromOrbit(orbitMc, params, results, torqueOpts);
disp('Torque Monte Carlo summary:');
disp(torqueMc.summary);

firstCase = find([torqueMc.runs.success], 1, 'first');
if ~isempty(firstCase)
    simParams = torqueMc.runs(firstCase).simParameters;
    simResults = torqueMc.runs(firstCase).simResults;
    [x, u, T_winf, x_est, indicators, sensors, actuators, error_flag] = ...
        core.simulation_rk4(simParams, simResults, false);
    fprintf('Attitude simulation error flag: %d\n', error_flag);
    disp(indicators);
end
