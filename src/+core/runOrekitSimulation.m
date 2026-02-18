function history = runOrekitSimulation(simParams)
%RUNOREKITSIMULATION High-fidelity orbit propagation using Orekit (MATLAB ↔ Java).
%
%   history = runOrekitSimulation(simParams)
%
%   OVERVIEW
%   --------
%   This function performs a high-fidelity numerical orbit propagation using
%   Orekit (Java) from within MATLAB. It sets up a NumericalPropagator with an
%   adaptive integrator and several perturbation force models, and logs the
%   propagation history at a fixed output step size using a custom Java
%   step handler (StepStorageHandler).
%
%   FORCE MODELS INCLUDED
%   ---------------------
%   1) Earth gravity field: normalized spherical harmonics (21×21)
%   2) Atmospheric drag: Harris–Priester model
%   3) Third-body gravity: Sun and Moon
%   4) Solar radiation pressure (SRP): isotropic model
%
%   GEOMAGNETIC FIELD
%   -----------------
%   The geomagnetic field is evaluated using the IGRF model (Orekit) and is
%   stored at each output step by StepStorageHandler as:
%       BN, BE, BD (North-East-Down components) in Tesla.
%
%   INPUT
%   -----
%   simParams : struct with fields:
%
%     initialDate : datetime
%         Simulation start date (UTC).
%
%     initialOrbit : struct (Keplerian elements)
%         .a    [m]    Semi-major axis
%         .e    [-]    Eccentricity
%         .i    [rad]  Inclination
%         .raan [rad]  Right ascension of ascending node
%         .pa   [rad]  Argument of perigee
%         .ta   [rad]  True anomaly
%
%     spacecraft : struct
%         .mass      [kg]  Spacecraft mass
%         .dragArea  [m^2] Drag reference area
%         .dragCoeff [-]   Drag coefficient
%         .srpArea   [m^2] SRP reference area
%         .srpCoeff  [-]   SRP reflectivity coefficient
%
%     duration : struct
%         .num_orbits [-] Number of Keplerian orbital periods to propagate
%
%     integrator : struct (Dormand–Prince 8(5,3))
%         .minStep      [s] Minimum integration step
%         .maxStep      [s] Maximum integration step
%         .posTolerance [m] Position tolerance (used for abs/rel tolerances)
%
%     output : struct
%         .stepSize [s] Fixed logging step size used by StepStorageHandler
%
%   OUTPUT
%   ------
%   history : double matrix [N × 22]
%       Each row corresponds to one fixed output time step.
%
%       Columns:
%         1      : time since initialDate [s]
%         2–4    : position (x, y, z) in inertial frame (GCRF) [m]
%         5–7    : Earth gravity acceleration [m/s^2]
%         8–10   : atmospheric drag acceleration [m/s^2]
%         11–13  : Sun third-body acceleration [m/s^2]
%         14–16  : Moon third-body acceleration [m/s^2]
%         17–19  : SRP acceleration [m/s^2]
%         20–22  : geomagnetic field (BN, BE, BD) in local NED [Tesla]
%
%   REQUIREMENTS
%   ------------
%   - Orekit must be correctly configured in MATLAB Java classpath.
%   - orekit-data must be available (gravity field, EOP, IGRF coefficients, etc.).
%   - StepStorageHandler.class must be compiled and accessible.
%
%   See also: StepStorageHandler, org.orekit.propagation.numerical.NumericalPropagator

    %% Import Orekit Java Classes
    % --- Propagation & Integration ---
    import org.orekit.propagation.numerical.NumericalPropagator;
    import org.orekit.propagation.SpacecraftState;
    import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
    % --- Time & Coordinates ---
    import org.orekit.frames.FramesFactory;
    import org.orekit.time.AbsoluteDate;
    import org.orekit.time.TimeScalesFactory;
    import org.orekit.orbits.PositionAngleType;
    import org.orekit.orbits.KeplerianOrbit;
    % --- Geomagnetic field ---
    import org.orekit.models.earth.GeoMagneticFieldFactory;
    import org.orekit.models.earth.GeoMagneticField;
    import org.orekit.models.earth.GeoMagneticElements;
    % --- Force Models & Bodies ---
    import org.orekit.forces.gravity.HolmesFeatherstoneAttractionModel;
    import org.orekit.forces.gravity.ThirdBodyAttraction;
    import org.orekit.forces.radiation.SolarRadiationPressure;
    import org.orekit.forces.radiation.IsotropicRadiationSingleCoefficient;
    import org.orekit.forces.drag.DragForce;
    import org.orekit.forces.drag.IsotropicDrag;
    import org.orekit.utils.Constants;
    import org.orekit.bodies.CelestialBodyFactory;
    import java.util.ArrayList;

    %% STEP 1: Configure numerical integrator
    disp('Configuring numerical propagator...');
    integrator = DormandPrince853Integrator(simParams.integrator.minStep, ...
                                           simParams.integrator.maxStep, ...
                                           simParams.integrator.posTolerance, ...
                                           simParams.integrator.posTolerance);

    %% STEP 2: Define initial state (epoch + Keplerian orbit)
    utc = TimeScalesFactory.getUTC();

    % Convert MATLAB datetime to Orekit AbsoluteDate (UTC)
    dt = simParams.initialDate;
    initialDate = AbsoluteDate(dt.Year, dt.Month, dt.Day, dt.Hour, dt.Minute, dt.Second, utc);

    % Build geomagnetic model (IGRF) at initial epoch (decimal year)
    decYear0  = GeoMagneticField.getDecimalYear(dt.Day, dt.Month, dt.Year);
    magModel0 = GeoMagneticFieldFactory.getIGRF(decYear0);

    inertialFrame = FramesFactory.getGCRF();
    p = simParams.initialOrbit;

    initialOrbit = KeplerianOrbit(p.a, p.e, p.i, p.pa, p.raan, p.ta, ...
                                  PositionAngleType.TRUE, inertialFrame, ...
                                  initialDate, Constants.WGS84_EARTH_MU);

    initialState = SpacecraftState(initialOrbit, simParams.spacecraft.mass);

    %% STEP 3: Create propagator
    propagator = NumericalPropagator(integrator);

    %% STEP 4: Add detailed force models
    force_models_list = ArrayList();
    sc = simParams.spacecraft;

    % Force 1: Earth gravity (21×21 normalized spherical harmonics)
    gravityProvider = org.orekit.forces.gravity.potential.GravityFieldFactory.getNormalizedProvider(21, 21);
    gravityForce = HolmesFeatherstoneAttractionModel(inertialFrame, gravityProvider);
    propagator.addForceModel(gravityForce);
    force_models_list.add(gravityForce);

    % Force 2: Atmospheric drag (Harris–Priester)
    sun = CelestialBodyFactory.getSun();
    earth = org.orekit.bodies.OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS, ...
        Constants.WGS84_EARTH_FLATTENING, ...
        FramesFactory.getITRF(org.orekit.utils.IERSConventions.IERS_2010, true));

    atmosphere = org.orekit.models.earth.atmosphere.HarrisPriester(sun, earth);
    dragForce = DragForce(atmosphere, IsotropicDrag(sc.dragArea, sc.dragCoeff));
    propagator.addForceModel(dragForce);
    force_models_list.add(dragForce);

    % Force 3: Third-body gravity (Sun and Moon)
    moon = CelestialBodyFactory.getMoon();
    sunAttraction = ThirdBodyAttraction(sun);
    moonAttraction = ThirdBodyAttraction(moon);

    propagator.addForceModel(sunAttraction);
    force_models_list.add(sunAttraction);

    propagator.addForceModel(moonAttraction);
    force_models_list.add(moonAttraction);

    % Force 4: Solar radiation pressure (isotropic)
    isotropicRadiation = IsotropicRadiationSingleCoefficient(sc.srpArea, sc.srpCoeff);
    srpForce = SolarRadiationPressure(sun, earth, isotropicRadiation);

    propagator.addForceModel(srpForce);
    force_models_list.add(srpForce);

    disp('Propagator created with high-fidelity force models.');
    fprintf('\n');

    %% STEP 5: Configure handler and run propagation
    fixedStepHandler = StepStorageHandler(force_models_list, initialDate, earth, magModel0);
    propagator.setStepHandler(simParams.output.stepSize, fixedStepHandler);

    propagation_duration = initialOrbit.getKeplerianPeriod() * simParams.duration.num_orbits;
    disp(['Propagating orbit for ', num2str(simParams.duration.num_orbits), ' orbits...']);

    propagator.setInitialState(initialState);
    finalState = propagator.propagate(initialDate.shiftedBy(propagation_duration)); %#ok<NASGU>
    disp('Propagation complete.');
    fprintf('\n');

    %% STEP 6: Retrieve and return results
    disp('Retrieving simulation history...');
    history_java = fixedStepHandler.getHistory();

    % Convert Java double[][] to MATLAB double matrix
    [rows, cols] = size(history_java);
    history = zeros(rows, cols);
    for r = 1:rows
        for c = 1:cols
            history(r,c) = history_java(r,c);
        end
    end
    disp('Data successfully retrieved.');
end