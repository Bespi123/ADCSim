% compile_java_handler.m
% -------------------------------------------------------------------------
% GOAL: Compile the Java step handler (StepStorageHandler.java)
%       and verify that the Orekit environment is accessible.
% EXECUTE: Only once, or whenever the .java file is modified.
% Tested by Bespi123
% -------------------------------------------------------------------------
clear; clc; clear java;

%% 1. Environment and Library Setup
disp('--- COMPILATION SCRIPT ---');
disp('Setting up Orekit environment...');
lib_folder = 'lib/orekit';
java_folder = 'lib/custom_java';

% --- Add all .jar libraries to the Java path ---
% Main Orekit library
javaaddpath(fullfile(lib_folder, 'orekit-13.1.2.jar'));
% Hipparchus library modules (version 4.0.2)
javaaddpath(fullfile(lib_folder, 'hipparchus-clustering-4.0.2.jar'));
javaaddpath(fullfile(lib_folder, 'hipparchus-core-4.0.2.jar'));
javaaddpath(fullfile(lib_folder, 'hipparchus-fft-4.0.2.jar'));
javaaddpath(fullfile(lib_folder, 'hipparchus-filtering-4.0.2.jar'));
javaaddpath(fullfile(lib_folder, 'hipparchus-fitting-4.0.2.jar'));
javaaddpath(fullfile(lib_folder, 'hipparchus-geometry-4.0.2.jar'));
javaaddpath(fullfile(lib_folder, 'hipparchus-ode-4.0.2.jar'));
javaaddpath(fullfile(lib_folder, 'hipparchus-optim-4.0.2.jar'));
javaaddpath(fullfile(lib_folder, 'hipparchus-samples-4.0.2.jar'));
javaaddpath(fullfile(lib_folder, 'hipparchus-stat-4.0.2.jar'));

disp('Libraries added to the Java path.');
fprintf('\n');

%% 2. Compile the Custom Step Handler
disp('Compiling StepStorageHandler.java...');

% Build the dynamic classpath that MATLAB is currently using
classpath_str = strjoin(javaclasspath('-dynamic'), pathsep);

% Create and execute the compilation command
compile_command = ['javac --release 8 -cp "', classpath_str, '" ', fullfile(java_folder,'StepStorageHandler.java')];
[status, cmdout] = system(compile_command);

% Check if the compilation was successful
if status == 0
    disp('SUCCESS! The StepStorageHandler.class file has been created correctly.');
    disp('You can now run the simulation script.');
else
    disp('--- JAVA COMPILATION ERROR ---');
    disp('Check the compiler output to find the error:');
    disp(cmdout);
    error('Could not compile StepStorageHandler.java. The simulation script will fail.');
end