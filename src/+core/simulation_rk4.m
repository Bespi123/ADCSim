function [x, u, T_winf_nosat, x_est, indicators, sensors, actuators, error_flag] = simulation_rk4(~,disturbances,simParameters,time)
%SIMULATION_RK4 Simulates CubeSat attitude dynamics with multi-rate sampling.
%
% This function simulates the attitude dynamics of a CubeSat using a Runge-Kutta 4
% integrator. It models a realistic digital control system where sensors (gyro,
% attitude sensors) and the control loop operate at different, independent
% sampling rates.
%
% Author: bespi123
% Creation Date: 2023-10-27
% Last Modified: 2025-10-13
%
% Inputs:
%   ~             - Placeholder for the UI application object, unused here.
%   disturbances  - (3xn) matrix of external disturbance torques over time.
%   simParameters - Struct containing all simulation parameters (initial values,
%                   controller gains, sensor specs, etc.).
%   time          - Struct defining simulation time parameters (.ti, .tf, .step).
%
% Outputs:
%   x             - (7xn) True CubeSat state history [q; w].
%   u             - (3xn) Commanded control torque history (updated at Ts_control).
%   T_winf_nosat  - (RW_num x n) Ideal torque command for each reaction wheel.
%   x_est         - (7xn) EKF estimated CubeSat state history [q_est; gyro_bias].
%   indicators    - Struct containing performance metrics (settling time, RMSE, etc.).
%   sensors       - Struct containing sensor models.
%   actuators     - Struct containing actuator models.
%   error_flag    - Flag indicating simulation stability (1 if NaN was detected).
%
% ----------------------------------------------------------------------------------

%%% Import adsim utilities
import adcsim.utils.*
%%% Import environment to use python libraries
%pyenv("Version", "D:\git\Satellite_ADCS_GUIDE\.venv\Scripts\python.exe"); 

%% 1. Parameter Recovery and Initialization
%%% Time parameters assuming uniform time steps
t = time.t;         % Full time vector.
n = time.nSteps;         % Total number of time points.
dt = t(2) - t(1);   % Simulation (integration) time step.

%%% Desired states (Setpoint)
% Replicated to have a desired value at each time step.
qd = repmat(simParameters.setPoint.qd',n, 1)';    % Desired quaternion.
wd = repmat(simParameters.setPoint.Wd',n, 1)';    % Desired angular velocity.

%%% Disturbances and Inertia
Td = disturbances; % External disturbance torque.
I = simParameters.spacecraft.I; % CubeSat's inertia matrix.

%%% Error flag
error_flag = 0; % Initialized to 0 (no error).

% Controller gains and initial values
% Gains are selected based on the controller chosen in the parameters.
if simParameters.controller.selector == 1
    %%% Feedback Controller gains
    P = simParameters.controller.feedback.Peye; 
    K = simParameters.controller.feedback.Keye;
else
    %%% Boskovic Controller (Adaptive) gains
    delta = simParameters.controller.boskController.delta;
    gamma = simParameters.controller.boskController.gamma;
    k     = simParameters.controller.boskController.k0;
    Umax  = simParameters.controller.boskController.Umax;
    %%% Previous values for integrating the gain 'k'.
    k_ant = k; k_dot_ant = 0;
end

% Retrieve sampling times from the parameters structure.
Ts_gyro = simParameters.sensors.gyro.Ts;
Ts_attitude = simParameters.sensors.attitude.Ts;
Ts_control = simParameters.control.Ts; % Sampling time for the control loop

% Calculate in how many simulation steps (dt) each sampling occurs.
steps_gyro = round(Ts_gyro / dt);
steps_attitude = round(Ts_attitude / dt);
steps_control = round(Ts_control / dt); % Steps for the control loop

% Ensure the number of steps is at least 1 to avoid division by zero or errors.
if steps_gyro < 1, steps_gyro = 1; end
if steps_attitude < 1, steps_attitude = 1; end
if steps_control < 1, steps_control = 1; end

%%% Star Sensor
if simParameters.sensors.star.enable == 1
    number_of_stars = simParameters.sensors.star.numberOfStars;
else
    number_of_stars = 0;
end

%%% Actuator parameters (Reaction Wheels)
number_of_rw = simParameters.rw.number; % Number of reaction wheels.
W  = simParameters.rw.W;                % Configuration matrix of the wheels.
motor = simParameters.rw.motor;         % DC motor model parameters.

%% 2. Simulation containers
% Pre-allocate memory with NaN for all variables to be logged.
% This improves performance and makes debugging easier (unwritten values remain NaN).
x = NaN(7,n); u = NaN(3,n); dq = NaN(4,n-1); T_winf_nosat = NaN(number_of_rw,n);
g_B = NaN(3,n-1); m_B = NaN(3,n-1); stars_B = NaN(3*number_of_stars,n-1);
omega_meas = NaN(3,n-1); omega_meas_filtered = NaN(3,n-1);
x_est = NaN(7,n-1); y_est = NaN(6+3*number_of_stars,n-3); 
o = NaN(1,n-1); 
x_rw = NaN(2*number_of_rw,n);
u_rw = NaN(number_of_rw,n);
w_cmd = NaN(number_of_rw,n);
torque_real = NaN(number_of_rw,n);

%% 3. Initial conditions
% Set initial true state, control, and estimated state
x(:,1)=[simParameters.initialValues.q0; simParameters.initialValues.Wo];
u(:,1) = zeros(3,1);  x_est(:,1) = [1, zeros(1,6)]'; % Initial estimate: identity quaternion, zero velocity/bias.
T_winf_nosat(:,1) = zeros(number_of_rw, 1);
x_rw(:,1) = repmat([motor.init.w_rw, motor.init.current],1,number_of_rw)';
u_rw(:,1) = zeros(number_of_rw,1);
w_cmd_ant = zeros(number_of_rw,1);

% Initialize containers for the last known measurement/calculation (for Zero-Order Hold).
last_u = zeros(number_of_rw,1);
last_T_winf_nosat = zeros(number_of_rw, 1);

% Initialize simulation objects from custom classes.
mySatellite = adcsim.satellite.Satellite(simParameters.initialValues, I);
reactionWheels = adcsim.actuators.ReactionWheelAssembly(simParameters.rw);
if (simParameters.sensors_selector.imu.flag)
    myIMU =  adcsim.sensors.IMU(simParameters.sensors);
else
    myIMU =  adcsim.sensors.IMU_Sophisticated(simParameters.sensors_soph);
end
myStarSensor = adcsim.sensors.StarTracker(simParameters.sensors.star);

% Initialize AHRS (Attitude and Heading Reference System) objects.
myEKF = adcsim.AHRS_algorithms.AttitudeEKF(simParameters.ahrs.ekf, dt, myIMU, myStarSensor);
myMadwick = adcsim.AHRS_algorithms.MadgwickAHRS(simParameters.ahrs.madwick, dt, myIMU, myStarSensor);
myUKF = adcsim.AHRS_algorithms.AttitudeUKF(simParameters.ahrs.ukf, dt, myIMU, myStarSensor);
myMahony = adcsim.AHRS_algorithms.Mahony(simParameters.ahrs.mahony, dt, myIMU, myStarSensor);
myQuest= adcsim.AHRS_algorithms.QUEST(simParameters.ahrs.quest, dt, myIMU, myStarSensor);
myRequest= adcsim.AHRS_algorithms.REQUEST(simParameters.ahrs.request, dt, myIMU, myStarSensor);
myCKF = adcsim.AHRS_algorithms.AttitudeCKF(simParameters.ahrs.ckf, dt, myIMU, myStarSensor);

% Get the initial rotation matrix from the true state.
R_init = quat2rot(x(1:4, 1));

% Simulate an initial sensor reading to initialize the attitude estimator.
g_B(:,1) = myIMU.getAccelerometerReading(R_init); 
m_B(:,1) = myIMU.getMagnetometerReading(R_init); 

% Initialize the EKF using the TRIAD method with the first measurements.
myEKF.initializeWithTriad(g_B(:,1), m_B(:,1));

if simParameters.sensors.star.enable == 1
    stars_B(:,1) = myStarSensor.getReading(R_init);
end

%% 4. Pre-Loop Calculations for Optimization
% The derivative of the desired angular velocity is pre-calculated to avoid doing it inside the loop.
Wd_dot = diff(wd')'./diff(t);
% Initialize the progress bar for user feedback.
hWaitbar = waitbar(0, 'Progress: 0%','Name', 'Attitude simulation progress..');

%% 5. Main Simulation Loop
for i = 1:n-1
    %% 5.1. Sensor models (add noise, bias and SAMPLING)
    % Gyroscope model with sampling.
    % The gyroscope is read only if the current step is a multiple of `steps_gyro`.
    if mod(i-1, steps_gyro) == 0 || i==1
        current_omega_meas  = myIMU.getGyroscopeReading(x(5:7, i));
        filtered_omega_meas = myIMU.filterGyroscopeReading(current_omega_meas);
    end
    % The measurement (raw or held from the previous step) is saved to the history.
    omega_meas(:, i) = current_omega_meas;
    omega_meas_filtered(:,i) = filtered_omega_meas;
    
    %% 5.2. AHRS - Attitude and Heading Reference System
    if simParameters.ahrs.enable
        if strcmp(simParameters.ahrs.flag, 'EKF')
            % --- EKF Prediction Step (always runs at the high frequency of the simulation) ---
            % The prediction step runs at each time step 'dt' using the
            % most recent (potentially held) gyroscope measurement.
            myEKF.predict(filtered_omega_meas);
            
            %%% --- EKF Correction Step (runs only when new attitude measurements are available) ---
            if mod(i-1, steps_attitude) == 0
                % Generate new attitude sensor measurements from the true state.
                R = quat2rot(x(1:4, i));
                g_B(:,i) = myIMU.getAccelerometerReading(R);
                m_B(:,i) = myIMU.getMagnetometerReading(R);
                if simParameters.sensors.star.enable == 1 
                    stars_B(:,i)  = myStarSensor.getReading(R);
                    y = [g_B(:,i); m_B(:,i); stars_B(:,i)]; % Measurement vector
                else
                    y = [g_B(:,i); m_B(:,i)];
                end
                % Perform the EKF correction step with the new measurements.
                myEKF.correct(y); 
            end
            % Save the updated state estimate.
            x_est(:, i + 1) = myEKF.x_est';
    
        elseif strcmp(simParameters.ahrs.flag, 'MADGWICK')
            %% Madgwick Algorithm
            myMadwick = myMadwick.Predict(filtered_omega_meas);
            % Correction step runs only when new attitude measurements are available.
            if mod(i-1, steps_attitude) == 0
                R = quat2rot(x(1:4, i));
                g_B(:,i) = myIMU.getAccelerometerReading(R);
                m_B(:,i) = myIMU.getMagnetometerReading(R);
                if simParameters.sensors.star.enable == 1 
                    stars_B(:,i)  = myStarSensor.getReading(R);
                    y = [g_B(:,i); m_B(:,i); stars_B(:,i)];
                else
                    y = [g_B(:,i); m_B(:,i)];
                end
                myMadwick = myMadwick.Update(filtered_omega_meas, y);
            end
            % Save the updated state estimate.
            x_est(1:4, i + 1) = myMadwick.x_est;
            x_est(5:7, i + 1) = myMadwick.gyro_bias;
        
        elseif strcmp(simParameters.ahrs.flag, 'MAHONY')
            myMahony = myMahony.Predict(filtered_omega_meas);
            % Correction step runs only when new attitude measurements are available.
            if mod(i-1, steps_attitude) == 0
                R = quat2rot(x(1:4, i));
                g_B(:,i) = myIMU.getAccelerometerReading(R);
                m_B(:,i) = myIMU.getMagnetometerReading(R);
                if simParameters.sensors.star.enable == 1 
                    stars_B(:,i)  = myStarSensor.getReading(R);
                    y = [g_B(:,i); m_B(:,i); stars_B(:,i)];
                else
                    y = [g_B(:,i); m_B(:,i)];
                end
                myMahony = myMahony.Update(filtered_omega_meas, y);
            end
            % Save the updated state estimate.
           x_est(1:4, i + 1) = myMahony.x_est;
           x_est(5:7, i + 1) = myMahony.gyro_bias; % Assuming UKF doesn't estimate gyro bias in this implementation.

        elseif strcmp(simParameters.ahrs.flag, 'QUEST')
             if mod(i-1, steps_attitude) == 0
                R = quat2rot(x(1:4, i));
                g_B(:,i) = myIMU.getAccelerometerReading(R);
                m_B(:,i) = myIMU.getMagnetometerReading(R);
                if simParameters.sensors.star.enable == 1 
                    stars_B(:,i)  = myStarSensor.getReading(R);
                    y = [g_B(:,i); m_B(:,i); stars_B(:,i)];
                else
                    y = [g_B(:,i); m_B(:,i)];
                end
                myQuest = myQuest.Update(filtered_omega_meas, y);
            end
            % Save the updated state estimate.
           x_est(1:4, i + 1) = myQuest.x_est;
           x_est(5:7, i + 1) = myQuest.gyro_bias; % Assuming UKF doesn't estimate gyro bias in this implementation.        

        elseif strcmp(simParameters.ahrs.flag, 'REQUEST')
            myRequest = myRequest.Predict(filtered_omega_meas);
            if mod(i-1, steps_attitude) == 0
                R = quat2rot(x(1:4, i));
                g_B(:,i) = myIMU.getAccelerometerReading(R);
                m_B(:,i) = myIMU.getMagnetometerReading(R);
                if simParameters.sensors.star.enable == 1 
                    stars_B(:,i)  = myStarSensor.getReading(R);
                    y = [g_B(:,i); m_B(:,i); stars_B(:,i)];
                else
                    y = [g_B(:,i); m_B(:,i)];
                end
                myRequest = myRequest.Update(filtered_omega_meas, y);
            end
            % Save the updated state estimate.
            x_est(1:4, i + 1) = myRequest.x_est;
            x_est(5:7, i + 1) = myRequest.gyro_bias; 

        elseif strcmp(simParameters.ahrs.flag, 'CKF')
            myCKF.Predict(filtered_omega_meas);
            if mod(i-1, steps_attitude) == 0
                R = quat2rot(x(1:4, i));
                g_B(:,i) = myIMU.getAccelerometerReading(R);
                m_B(:,i) = myIMU.getMagnetometerReading(R);
                if simParameters.sensors.star.enable == 1 
                    stars_B(:,i)  = myStarSensor.getReading(R);
                    y = [g_B(:,i); m_B(:,i); stars_B(:,i)];
                else
                    y = [g_B(:,i); m_B(:,i)];
                end
                myCKF.Correct(y);
            end
            % Save the updated state estimate.
            x_est(1:4, i + 1) = myCKF.x_est;
            x_est(5:7, i + 1) = myCKF.gyro_bias; 
            
        else
            %% Unscented Kalman Filter (UKF)
            myUKF.Predict(filtered_omega_meas);
            % Correction step runs only when new attitude measurements are available.
            if mod(i-1, steps_attitude) == 0
                R = quat2rot(x(1:4, i));
                g_B(:,i) = myIMU.getAccelerometerReading(R);
                m_B(:,i) = myIMU.getMagnetometerReading(R);
                if simParameters.sensors.star.enable == 1 
                    stars_B(:,i)  = myStarSensor.getReading(R);
                    y = [g_B(:,i); m_B(:,i); stars_B(:,i)];
                else
                    y = [g_B(:,i); m_B(:,i)];
                end
                myUKF.Correct(y);
            end
            % Save the updated state estimate.
           x_est(1:4, i + 1) = myUKF.x_est;
           x_est(5:7, i + 1) = zeros(3,1); % Assuming UKF doesn't estimate gyro bias in this implementation.
        end
    end
    %% 5.3. Control Law (runs at its own sampling rate Ts_control)
    % A new control command is calculated only when it's the correct time step.
    if mod(i-1, steps_control) == 0
        % --- Calculate new control command ---
        % Use either the estimated state or the ideal true state as feedback.
        if simParameters.ahrs.enable_feedback && simParameters.ahrs.enable
            % Use sensor model output (estimated state) as feedback signal.
            % Gyro bias is subtracted from the filtered measurement.
            feed_est = [x_est(1:4,i); omega_meas_filtered(:,i)-x_est(5:7,i)];
        else
            % Use ideal (true) signals as feedback for baseline comparison.
            feed_est = x(:, i);
        end
        
        % Calculate the error quaternion.
        dq(:, i) = Error_quaternio(qd(:, i), feed_est(1:4));
        
        % Measure the computation time of the control law.
        tic
        if simParameters.controller.selector == 1
            % Standard Feedback Controller.
            u_new = adcsim.controllers.ControlFeedback_rw(I, feed_est, dq(:, i), wd(:, i), Wd_dot(:, i), P, K); 
        else
            % Boskovic Adaptive Controller.
            k_dot = adcsim.controllers.Gain_estimator_bosk(feed_est(5:7), wd(:, i), dq(:, i), delta, gamma, k, Umax);
            k = k_ant + Ts_control / 6 * (k_dot_ant + 2 * (k_dot_ant + k_dot) + k_dot); % Integration step for adaptive gain
            u_new = adcsim.controllers.Boskovic_control(feed_est(5:7), wd(:, i), dq(:, i), delta, k, Umax);
        end
        o(i) = toc; % Log computation time.
        
        % --- Torque Allocator ---
        % Distributes the desired 3-axis torque `u_new` among the available reaction wheels.
        if number_of_rw == 3
            % For 3 orthogonal RWs, the solution is a direct matrix inversion.
            T_winf_new = W\u_new;
        else
            % For more than 3 RWs (redundant system), use optimization.
            T_w_L2norm = allocator_L2norm(W, u_new); % Pseudoinverse (minimum energy)
            T_winf_new = allocator_LinfNorm(T_w_L2norm); % Another allocation strategy
        end
        
        % --- Update last known values for Zero-Order Hold (ZOH) ---
        last_u = u_new;
        last_T_winf_nosat = T_winf_new;
        % --- Calculate commanded angular rate for each reaction wheel ---
        last_w_rw_cmd = w_cmd_ant + (last_T_winf_nosat / motor.Jrw) * Ts_control;
    else
        % --- Hold previous control command (ZOH) ---
        % If it's not a control step, do not compute a new command. The previous command is held.
        o(i) = NaN; % No control computation cost for this step.
    end
    
    % Log the control signal for every simulation step.
    % This will result in a piecewise constant signal due to the ZOH.
    u(:, i + 1) = last_u;
    T_winf_nosat(:, i + 1) = last_T_winf_nosat;
    w_cmd(:, i) = last_w_rw_cmd;

    %% Actuator model
    % Updates the reaction wheels' state based on the command and returns the actual torque produced.
    [reactionWheels, torque_real_vector] = reactionWheels.update(w_cmd(:, i), dt);
    x_rw(:, i + 1) = reactionWheels.States(:);     % Log RW states (speed, current)
    u_rw(:, i + 1) = reactionWheels.voltage_vector(:); % Log applied voltage
    torque_real(:,i) = torque_real_vector';      % Log actual produced torque
    
    % Calculate the total torque from the actuators projected onto the satellite body axes.
    T_u  = W*torque_real(:,i);
    
    %% Satellite Dynamics
    % Integrate the satellite's dynamics forward by one time step `dt`.
    mySatellite = mySatellite.updateState(T_u, Td(:, i), dt);
    x(:, i + 1) = mySatellite.State; % Log the new true state.
    
    % Check if the simulation has become unstable (produced a NaN).
    if isnan(x(:, i + 1))
        error_flag = 1;
        break;
    end

    %% Progress bar update
    if mod(i, round(n / 10)) == 0
        progress = i / (n - 1);
        if isa(hWaitbar, 'handle') && isvalid(hWaitbar)
            waitbar(progress, hWaitbar, sprintf('Progress: %.1f%%', progress * 100));
        end
    end
    % Update the "previous" commanded wheel speed for the next iteration's calculation.
    w_cmd_ant = w_cmd(:,i);
end

% --- Package sensor data for output ---
% Fill missing values using the 'previous' method to create continuous signals for plotting.
sensors.meas = [omega_meas; fillmissing([g_B; m_B; stars_B], 'previous',2)];
dq_array = fillmissing(dq, 'previous',2);
sensors.est = fillmissing(y_est, 'previous',2);
sensors.omega_filtered = omega_meas_filtered;

% --- Package actuator data for output ---
actuators.w_cmd = w_cmd;
actuators.x_rw  = x_rw;
actuators.torque_real = torque_real;
actuators.u_rw  = u_rw;
close(hWaitbar);

%% 6. Update performance indices
% Calculate performance indicators if no error occurred during the simulation.
if error_flag == 0
    % Convert error quaternion to angle-axis representation to get the pointing error angle.
    ang_and_axis = quat2axang(dq_array'); 
    eulerang = ang_and_axis(:,4); % The 4th element is the angle of rotation.
    
    % Integral of the angle error (IAE).
    indicators.eulerInt = cumtrapz(t(1:end-1),eulerang); 
    
    % Accumulated squared control torque (approximates control effort).
    ascct_dot=vecnorm(u,2,1).^2;
    indicators.ascct = cumtrapz(t,ascct_dot); 
    
    % Settling time calculation.
    tol = 5/100*(180/pi*quat2eul(simParameters.setPoint.qd')); % 5% tolerance
    euler_angles_deg = 180/pi*quat2eul(dq_array');
    indicators.ts(1) = calculateSettlementTime(euler_angles_deg(:,1), t, tol(1));
    indicators.ts(2) = calculateSettlementTime(euler_angles_deg(:,2), t, tol(2));
    indicators.ts(3) = calculateSettlementTime(euler_angles_deg(:,3), t, tol(3));

    % Fill missing computation times.
    indicators.o = fillmissing(o, 'previous');
    
    % Root Mean Square Error (RMSE) of the attitude estimation.
    if simParameters.ahrs.enable
        ang_and_axis_real = quat2axang(x(1:4,:)');
        ang_and_axis_est = quat2axang(x_est(1:4,:)');
        error_vec = ang_and_axis_real(:,4) - ang_and_axis_est(:,4);
        indicators.RMSE = sqrt(mean(error_vec.^2));         
    else 
        indicators.RMSE = 0.00;
    end
else
    % If an error occurred, show a dialog box and return NaNs for indicators.
    errordlg('Unstable System', 'ERROR');
    indicators = NaN;
end
end