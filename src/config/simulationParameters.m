function parameters = simulationParameters(~)
    import org.orekit.utils.Constants;

    % This function defines all configurable parameters for the CubeSat ADCS simulation.
    % It returns a single structure 'parameters' that contains all nested sub-structures
    % for different components of the simulation (e.g., initial values, sensors, control).

    % --- CubeSat Physical Parameters ---
    % Inertia tensor of the CubeSat [kg*m^2]. This represents the distribution of mass
    % and how the satellite resists rotational motion. Off-diagonal terms represent
    % products of inertia.
    parameters.initialValues.I = 1E-3 * [3.084 0.082 -0.054;
                                         0.082 3.132 0.016;
                                        -0.054 0.016 3.540];
    % Initial attitude quaternion [w, x, y, z]' (scalar-first format). 
    % Represents the initial orientation of the CubeSat relative to the inertial frame.
    parameters.initialValues.q0 = [1; 0; 0; 0]; % Starts with no rotation.
    
    % Initial angular velocity of the CubeSat body [rad/s].
    parameters.initialValues.Wo = [0; 0; 0]; % Starts stationary.
    
    % --- System Timing and Sampling Parameters ---
    % These parameters simulate a digital system where sensors and control loops
    % run at discrete time intervals, not continuously.
    
    % Sample time for the gyroscope sensor [s]. Gyros are often sampled at a higher
    % rate because they are used for high-frequency state propagation in the AHRS.
    parameters.sensors.gyro.Ts      = 0.005; % 200 Hz
    
    % Sample time for other attitude sensors (Accelerometer, Magnetometer, Star Tracker) [s].
    % These sensors are typically slower and used for correcting the gyro drift.
    parameters.sensors.attitude.Ts  = 0.01;  % 100 Hz
    
    % Sample time for the control loop execution [s]. This defines how often the
    % control law is calculated and a new command is sent to the actuators.
    parameters.control.Ts           = 0.005; % 200 Hz
    
    % --- Sensor Model Parameters ---
    % These sections model the imperfections of real-world sensors.
    % Flag to use sophisticated model (0) or simple model (1)
    parameters.sensors_selector.imu.flag = 1;

    %%% Accelerometer sensor
    % Gravitational acceleration vector in the inertial frame (e.g., pointing towards Earth).
    % This is the 'true' vector that the accelerometer measures to help determine 'down'.
    parameters.sensors.acc.g_I  = [0,0,1]';  
    % Accelerometer constant bias [m/s^2]. A persistent offset in the measurement.
    parameters.sensors.acc.bias = [0.013127156; 0.010186677; 0.000138897];
    % Accelerometer noise standard deviation [m/s^2]. Models random, fluctuating errors.
    parameters.sensors.acc.std  = [0.011817097; 0.012306208; 0.015949601];
    
    %%% Magnetometer sensor
    % Earth's magnetic field vector in the inertial frame (e.g., pointing towards magnetic north).
    % This is the 'true' vector that the magnetometer measures to get a heading reference.
    parameters.sensors.mag.m_I  = [0.7071 ,0 ,0.7071]';  
    % Magnetometer constant bias [uT or normalized units].
    parameters.sensors.mag.bias = [0; 0; 0];
    % Magnetometer noise standard deviation [uT or normalized units].
    parameters.sensors.mag.std  = [0.011817097; 0.012306208; 0.015949601];
    
    %%% Gyroscope sensor
    % Gyroscope constant bias [rad/s]. Represents a slow drift in the angular rate measurement.
    parameters.sensors.gyro.bias = [-0.0000541; -0.0001553; -0.00006458];
    % Gyroscope noise standard deviation (Angular Random Walk) [rad/s]. Models high-frequency noise.
    parameters.sensors.gyro.std  = [0.009; 0.008; 0.007];
    % Gyroscope low-pass filter time constant [s]. Models the sensor's response time.
    parameters.sensors.gyro.tau = 0.1;
    
    %%% Star Tracker sensor (high-precision attitude sensor)
    % Star Tracker constant bias [rad].
    parameters.sensors.star.bias = 1E-6 * [1; 1; 1];
    % Star Tracker noise standard deviation [rad].
    parameters.sensors.star.std  = 1E-3 * [1; 1; 1];
    % Flag to enable (1) or disable (0) the star tracker in the simulation.
    parameters.sensors.star.enable = 1;
    % Number of stars (or direction vectors) to be simulated by the star tracker.
    parameters.sensors.star.numberOfStars = 10;

    % --- Sophisticated Sensor Model Parameters ---
    % These parameters are for the IMU_Sophisticated class, which includes
    % more complex error sources like bias instability, scale factor, and misalignment.

    %%% Gyroscope (Sophisticated Model)
    % Initial bias at the start of the simulation [rad/s].
    parameters.sensors_soph.gyro.initial_bias = [-0.0000541; -0.0001553; -0.00006458];
    % Gyroscope low-pass filter time constant [s]. Models the sensor's response time.
    parameters.sensors_soph.gyro.tau = parameters.sensors.gyro.tau;
    % Standard deviation of the high-frequency white noise [rad/s].
    parameters.sensors_soph.gyro.std = [0.009; 0.008; 0.007];
    % Bias Instability (Random Walk) standard deviation [rad/s/sqrt(Hz) or rad/sqrt(s)].
    % This value dictates how fast the bias drifts. (e.g., from a 10 deg/hr spec).
    parameters.sensors_soph.gyro.std_RW = 4.8e-5 * [1; 1; 1];
    % Scale factor errors (e.g., 0.1% error) [dimensionless].
    parameters.sensors_soph.gyro.scale_factor_errors = 0.001 * [0.5; -0.8; 0.3];
    % Axis misalignment errors (e.g., 0.2 degrees) [radians].
    % [m12, m13, m21, m23, m31, m32]
    parameters.sensors_soph.gyro.misalignment_errors = 0.0035 * [1; -0.5; 0.8; 1.1; -0.9; 0.6];
    
    %%% Accelerometer (Sophisticated Model)
    parameters.sensors_soph.acc.g_I = parameters.sensors.acc.g_I;
    % Constant bias [m/s^2].
    parameters.sensors_soph.acc.bias = [0.013127156; 0.010186677; 0.000138897];
    % Standard deviation of the high-frequency white noise [m/s^2].
    parameters.sensors_soph.acc.std  = [0.011817097; 0.012306208; 0.015949601];
    % Scale factor errors (e.g., 0.5% error) [dimensionless].
    parameters.sensors_soph.acc.scale_factor_errors = 0.005 * [-0.7; 0.4; 0.9];
    % Axis misalignment errors (e.g., 0.5 degrees) [radians].
    parameters.sensors_soph.acc.misalignment_errors = 0.0087 * [1; 1; 1; 1; 1; 1];

    %%% Magnetometer (Sophisticated Model)
    parameters.sensors_soph.mag.m_I = parameters.sensors.mag.m_I;
    % Constant bias [uT or normalized units].
    parameters.sensors_soph.mag.bias = [0; 0; 0];
    % Standard deviation of the high-frequency white noise [uT or normalized units].
    parameters.sensors_soph.mag.std  = [0.011817097; 0.012306208; 0.015949601];
    % Scale factor errors (magnetometers are often less accurate, e.g., 2% error) [dimensionless].
    parameters.sensors_soph.mag.scale_factor_errors = 0.02 * [1.2; -0.8; 0.5];
    % Axis misalignment errors (e.g., 1.0 degree) [radians].
    parameters.sensors_soph.mag.misalignment_errors = 0.0175 * [1; 1; 1; 1; 1; 1];
    
    % --- AHRS (Attitude and Heading Reference System) Selection ---
    % The AHRS is the algorithm that fuses sensor data to estimate the attitude.
    % Selector for the AHRS algorithm: 'MAHONY', 'MADGWICK', 'EKF', or 'UKF'.
    %parameters.ahrs.flag = 'EKF';
    parameters.ahrs.flag = 'MAHONY';
    % Flag: 1 enables AHRS to the controller (realistic).
    %       0 don't use AHRS algorithm.
    parameters.ahrs.enable = true;
    % Flag: 1 enables feedback from the AHRS to the controller (realistic).
    %       0 uses the ideal, true state for feedback (for baseline comparison).
    parameters.ahrs.enable_feedback = true;
    
    % MAHONY filter parameters
    parameters.ahrs.mahony.K_P = 0.5;   % Proportional gain
    parameters.ahrs.mahony.K_I = 0.001; % Integral gain
    parameters.ahrs.mahony.K_1 = 1;     % Acc gain
    parameters.ahrs.mahony.K_2 = 1;   % Gyro gain
    parameters.ahrs.mahony.K_3 = 1;   % Star tracker gain
    parameters.ahrs.mahony.dt = 0.01; % Filter sampling time [s]

    % MADGWICK filter parameters
    parameters.ahrs.madwick.beta = 1; % Filter gain.
    parameters.ahrs.madwick.bias_gain = 0.01; % Gain for gyro bias estimation.
    
    % EKF (Extended Kalman Filter) parameters
    % The EKF needs a model of the sensor noise (R matrix).
    % Flag: 1 sets the EKF's noise parameters equal to the sensor model's noise (above).
    %       0 allows you to specify different noise values, e.g., if you want to
    %       intentionally mismatch the model to test robustness.
    parameters.ahrs.ekf.equalModel = 0;
    % EKF accelerometer noise standard deviation (used in R matrix).
    parameters.ahrs.ekf.acc_std    = [0.011817097; 0.012306208; 0.015949601];
    % EKF magnetometer noise standard deviation (used in R matrix).
    parameters.ahrs.ekf.mag_std    = [0.011817097; 0.012306208; 0.015949601];
    % EKF gyroscope noise standard deviation (used in Q matrix - process noise).
    parameters.ahrs.ekf.gyro_std   = [0.009; 0.008; 0.007];
    % EKF star tracker noise standard deviation (used in R matrix).
    parameters.ahrs.ekf.star_std   = 1E-3 * [1; 1; 1];
    
    % UKF (Unscented Kalman Filter) parameters
    % These parameters control the generation and weighting of sigma points.
    parameters.ahrs.ukf.alpha = 1e-3; % Scaling parameter, spread of sigma points.
    parameters.ahrs.ukf.beta  = 2.0;  % Used to incorporate prior knowledge of the state distribution.
    parameters.ahrs.ukf.kappa = 0.0;  % Secondary scaling parameter.
    parameters.ahrs.ukf.gyro_std = [0.009; 0.008; 0.007];
    parameters.ahrs.ukf.acc_std  = [0.011817097; 0.012306208; 0.015949601];
    parameters.ahrs.ukf.mag_std  = [0.011817097; 0.012306208; 0.015949601];
    parameters.ahrs.ukf.star_std  = 1E-3 * [1; 1; 1];
    
    % --- Controller Setpoint Parameters ---
    % The 'setpoint' is the target state that the controller tries to achieve.
    % Desired target attitude quaternion [w, x, y, z]'.
    parameters.setPoint.qd = [0.9808; 0.1691; 0.0933; 0.0277];
    % Desired target angular velocity [rad/s].
    parameters.setPoint.Wd = [0; 0; 0]; % e.g., hold a steady attitude.
    
    % --- Controller Selection ---
    % Selector: 1 for a standard PD Feedback Controller, 2 for Boskovic's Adaptive Controller.
    parameters.controller.selector = 1;
    
    % --- Feedback Controller (PD) Parameters ---
    % Derivative gain matrix (acts on angular velocity error).
    parameters.controller.feedback.Keye = 0.2 * eye(3);
    % Proportional gain matrix (acts on quaternion error).
    parameters.controller.feedback.Peye = 0.3 * eye(3);
    
    % --- Boskovic Adaptive Controller Parameters ---
    parameters.controller.boskController.delta = 0.5;
    parameters.controller.boskController.gamma = 0.001; % Adaptation rate.
    parameters.controller.boskController.Umax  = 0.10;  % Maximum control torque magnitude [Nm].
    parameters.controller.boskController.k0    = 1;     % Initial gain value.
    
    % --- Disturbance Torque Parameters ---
    % These model external torques acting on the satellite (e.g., solar pressure, drag).
    % Constant disturbance torque components [Nm].
    parameters.disturbances.constant.Tdx = 0.001;
    parameters.disturbances.constant.Tdy = 0.001;
    parameters.disturbances.constant.Tdz = 0.001;
    % Flag to enable (2) or disable (1) constant disturbances.
    parameters.disturbances.constant.enable = 1;
    
    % Parameters for variable (sinusoidal) disturbances [Nm].
    parameters.disturbances.variable.A   = 0.01;   % Amplitude
    parameters.disturbances.variable.w1  = 1.00;   % Frequency for x-axis
    parameters.disturbances.variable.w2  = 2.00;   % Frequency for y-axis
    parameters.disturbances.variable.w3  = 3.00;   % Frequency for z-axis
    % Flag to enable (2) or disable (1) variable disturbances.
    parameters.disturbances.variable.enable = 1;
    
    % --- Reaction Wheels (Actuators) Configuration ---
    % Number of reaction wheels in the assembly.
    parameters.rw.number = 3;
    % Maximum momentum storage per wheel [Nms]. A physical limitation.
    parameters.rw.maxH   = 1;
    % Mounting angles for wheels (beta1 and beta2 are typically used for cant angles).
    parameters.rw.beta1  = 0;
    parameters.rw.beta2  = 0;
    % Distribution matrix 'W'. Maps the 3-axis desired body torque to the
    % individual torque commands for each wheel based on their mounting orientation.
    % For a standard 3-orthogonal-wheel setup, this is the identity matrix.
    parameters.rw.W      = eye(3);
    
    % --- Reaction Wheel DC Motor Model ---
    % These parameters define the electrical and mechanical properties of the motors.
    parameters.rw.motor.kt   = 0.0254;   % Torque constant [Nm/A].
    parameters.rw.motor.Jrw  = 3.475E-5; % Rotor inertia [kg*m^2].
    parameters.rw.motor.b    = 5.58E-4;  % Viscous friction coefficient [Nms].
    parameters.rw.motor.c    = 0.0000;   % Coulomb friction (static friction) [Nm].
    parameters.rw.motor.L    = 0.000403; % Armature inductance [H].
    parameters.rw.motor.R    = 1.16;     % Armature resistance [Ohms].
    parameters.rw.motor.ke   = 0.0254;   % Back-EMF (voltage) constant [Vs/rad].
    
    % --- Reaction Wheel DC motor initial conditions ---
    parameters.rw.motor.init.w_rw = 0;      % Initial wheel speed [rad/s].
    parameters.rw.motor.init.current = 0; % Initial motor current [A].
    
    % --- Reaction Wheel PID Speed Controller Gains ---
    % This internal PID controller's job is to ensure the wheel's actual speed
    % matches the commanded speed from the ADCS control law.
    parameters.rw.motor.pid.kp = 0.504519552012394;  % Proportional gain
    parameters.rw.motor.pid.ki = 28.7345189644697;   % Integral gain
    parameters.rw.motor.pid.kd = -0.000554999924740984; % Derivative gain
    parameters.rw.motor.pid.Nu = 289.873462888377;  % Derivative filter coefficient.

    %% Orekit parameters
    % --- Initial Date ---
    parameters.initialDate = datetime(2025, 10, 10, 23, 30, 00);
    
    % --- Initial Orbit (Keplerian Elements) ---
    parameters.initialOrbit.a    = Constants.WGS84_EARTH_EQUATORIAL_RADIUS + 420e3; % Semi-major axis (m)
    parameters.initialOrbit.e    = 0.0001;          % Eccentricity
    parameters.initialOrbit.i    = 51.6 * pi / 180; % Inclination (rad)
    parameters.initialOrbit.raan = 0.0;             % Right Ascension of the Ascending Node (rad)
    parameters.initialOrbit.pa   = 0.0;             % Argument of Perigee (rad)
    parameters.initialOrbit.ta   = 0.0;             % True Anomaly (rad)
    
    % --- Spacecraft Physical Properties ---
    parameters.spacecraft.mass       = 10.0;  % kg
    parameters.spacecraft.dragArea   = 1.0;   % m^2
    parameters.spacecraft.dragCoeff  = 2.2;   % unitless
    parameters.spacecraft.srpArea    = 1.0;   % m^2
    parameters.spacecraft.srpCoeff   = 1.2;   % unitless (Cr)
    
    % --- Simulation Duration & Output ---
    parameters.duration.num_orbits = 3;
    parameters.output.stepSize     = 60.0; % Time step for data logging (seconds)
    
    % --- Integrator Settings ---
    parameters.integrator.minStep       = 0.001;  % s
    parameters.integrator.maxStep       = 1000;   % s
    parameters.integrator.posTolerance  = 1.0;    % m

end