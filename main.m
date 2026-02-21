%% 
% =========================================================================
%               MAIN SCRIPT TO LAUNCH THE ADCS SIMULATION APP
% =========================================================================
% This script prepares the MATLAB environment and then launches the main
% GUIDE application.
%
% TO RUN THE PROGRAM: Simply execute this file.
% Authot: bespi123
% -------------------------------------------------------------------------

%% 1. Environment Cleanup
% Start from a clean slate to avoid conflicts with previous sessions.
clear;          % Clear all variables from the workspace
clc;            % Clear the command window
close all;      % Close all open figure windows

disp('Initializing ADCS Simulation Application...');

%% 2. Path Setup
% Add all necessary folders and their subfolders to the MATLAB path.
% This allows MATLAB to find all your functions, classes, and libraries.

disp('Adding project folders to the path...');

% Get the path of the current folder where main.m is located.
project_root = fileparts(mfilename('fullpath'));

% Add the 'src' folder (your source code) and all its subfolders.
addpath(genpath(fullfile(project_root, 'src')));

% Add the 'lib' folder (external libraries) and all its subfolders.
addpath(genpath(fullfile(project_root, 'lib')));

disp('Environment is ready.');
fprintf('\n');

%% 3. Launch the GUIDE Application
% Call the main GUIDE file by its name to open the graphical interface.
% !!! IMPORTANT: Replace 'ADCS_program' with the actual name of your GUIDE file.
ADCS_program;

disp('GUI launched successfully.');