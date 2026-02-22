%% 1. Environment Setup
disp('Setting up Orekit environment...');

% Define the paths to the library folders (adjust if necessary)
orekit_lib_folder  = 'lib/orekit';
orekit_data_folder = 'lib/orekit/orekit-data';
java_lib_folder    = 'lib/custom_java';

% --- Add all required .jar libraries to the MATLAB Java path ---
if exist('org.orekit.propagation.numerical.NumericalPropagator', 'class') ~= 8
    % Main Orekit library
    javaaddpath(fullfile(orekit_lib_folder, 'orekit-13.1.2.jar'));
    % Hipparchus library modules (version 4.0.2)
    javaaddpath(fullfile(orekit_lib_folder, 'hipparchus-clustering-4.0.2.jar'));
    javaaddpath(fullfile(orekit_lib_folder, 'hipparchus-core-4.0.2.jar'));
    javaaddpath(fullfile(orekit_lib_folder, 'hipparchus-fft-4.0.2.jar'));
    javaaddpath(fullfile(orekit_lib_folder, 'hipparchus-filtering-4.0.2.jar'));
    javaaddpath(fullfile(orekit_lib_folder, 'hipparchus-fitting-4.0.2.jar'));
    javaaddpath(fullfile(orekit_lib_folder, 'hipparchus-geometry-4.0.2.jar'));
    javaaddpath(fullfile(orekit_lib_folder, 'hipparchus-ode-4.0.2.jar'));
    javaaddpath(fullfile(orekit_lib_folder, 'hipparchus-optim-4.0.2.jar'));
    javaaddpath(fullfile(orekit_lib_folder, 'hipparchus-samples-4.0.2.jar'));
    javaaddpath(fullfile(orekit_lib_folder, 'hipparchus-stat-4.0.2.jar'));
    
    % --- Add Custom Java Code to the Java Path ---
    % This tells MATLAB where to find our pre-compiled 'StepStorageHandler.class' file.
    javaaddpath(java_lib_folder);
    
    % --- VERIFICATION STEP ---
    % Check if a core Orekit class and our custom handler are visible to MATLAB.
    if exist('org.orekit.propagation.numerical.NumericalPropagator', 'class') ~= 8
        % If the Orekit class is not found, show an error dialog and stop.
        errordlg('Orekit library not found. Check the path in orekit_lib_folder.', 'Setup Error');
        error('Orekit library not found.'); % Also stop script execution with a command-line error.
    end
    if exist('StepStorageHandler', 'class') ~= 8
        % If the custom handler is not found, show an error dialog and stop.
        errordlg('Custom handler StepStorageHandler.class not found. Recompile the .java file.', 'Setup Error');
        error('Custom handler not found.');
    end
    
    % --- SUCCESS MESSAGE ---
    % Display a pop-up window to confirm that all libraries were loaded correctly.
    msgbox('Libraries and custom handler loaded successfully.', 'Setup Complete', 'help');
else
    msgbox('Orekit is already in the Java Path.', 'Setup Complete', 'help');
end

% --- Configure the Orekit Data Provider ---
% This tells Orekit where to find essential data like Earth gravity fields, time scales, etc.
import org.orekit.data.*;
import java.io.File;
manager = DataContext.getDefault().getDataProvidersManager();
manager.clearProviders(); % Good practice to prevent duplicate paths.
manager.addProvider(DirectoryCrawler(File(orekit_data_folder)));

disp('Environment ready.');
fprintf('\n');