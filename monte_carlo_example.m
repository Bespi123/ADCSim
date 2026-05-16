%% Monte Carlo example for ADCSim
% Run this script from the project root after configuring the MATLAB path.

clear;
clc;
close all;

project_root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(project_root, 'src')));
addpath(genpath(fullfile(project_root, 'lib')));

params = simulation_Parameters([]);
results = initialize_Results(params);

opts = struct();
opts.numRuns = 50;
opts.seed = 1234;
opts.showProgress = true;
opts.collectTrajectories = false;

opts.perturbations.initialAttitudeStdRad = deg2rad(1.0);
opts.perturbations.initialRateStdRadSec = deg2rad([0.05; 0.05; 0.05]);
opts.perturbations.inertiaRelStd = 0.03;
opts.perturbations.gyroBiasStdRadSec = deg2rad([0.005; 0.005; 0.005]);
opts.perturbations.accBiasStd = [0.002; 0.002; 0.002];
opts.perturbations.magBiasStd = [0.002; 0.002; 0.002];
opts.perturbations.rwMotorRelStd = 0.02;
opts.perturbations.disturbanceTorqueStdNm = 1e-5;

mc = core.runMonteCarloSimulation(params, results, opts);

disp('Monte Carlo summary:');
disp(mc.summary);

figure('Name', 'Monte Carlo RMSE');
histogram(mc.metrics.RMSE(mc.metrics.success));
xlabel('RMSE');
ylabel('Runs');
grid on;

