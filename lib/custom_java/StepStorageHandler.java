import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.sampling.OrekitFixedStepHandler;
import org.orekit.forces.ForceModel;
import org.orekit.time.AbsoluteDate;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import java.util.ArrayList;
import java.util.List;

// --- Frames + Tierra + Geomag ---
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.Transform;
import org.orekit.utils.IERSConventions;
import org.orekit.bodies.OneAxisEllipsoid;
import org.orekit.bodies.GeodeticPoint;

import org.orekit.models.earth.GeoMagneticField;
import org.orekit.models.earth.GeoMagneticElements;

/**
 * A custom Orekit step handler designed to store detailed simulation data for analysis in MATLAB.
 *
 * Stores at each fixed step:
 *  - time since initialDate [s]
 *  - position (x,y,z) in the state's frame
 *  - acceleration (ax,ay,az) for each ForceModel (up to 5 forces -> 15 values)
 *  - geomagnetic field components (BN, BE, BD) in local NED [Tesla]
 *
 * Total columns:
 *  1 (time) + 3 (pos) + 5*3 (forces) + 3 (B) = 22
 */
public class StepStorageHandler implements OrekitFixedStepHandler {

    // --- CONFIG: expected forces block size ---
    // Indices:
    // 0    = time
    // 1..3 = position
    // 4..18 = force accelerations (max 5 forces * 3)
    // 19..21 = BN, BE, BD
    private static final int FORCES_BLOCK_END_EXCLUSIVE = 19; // last valid force index is 18
    private static final int N_COLS = 22;

    // --- CLASS FIELDS ---
    private final List<ForceModel> forceModels;
    private final AbsoluteDate initialDate;
    private final ArrayList<double[]> history = new ArrayList<>();

    // --- Geomagnetic dependencies ---
    private final OneAxisEllipsoid earth;
    private final Frame itrf;
    private final GeoMagneticField magModel;

    /**
     * Constructs the step handler.
     *
     * @param forceModels The list of ForceModel objects whose accelerations will be recorded.
     * @param initialDate The absolute start date of the simulation.
     * @param earth       OneAxisEllipsoid Earth model (used to convert ITRF position to lat/lon/alt).
     * @param magModel    GeoMagneticField model instance (e.g., WMM or IGRF) for a given epoch.
     */
    public StepStorageHandler(final List<ForceModel> forceModels,
                              final AbsoluteDate initialDate,
                              final OneAxisEllipsoid earth,
                              final GeoMagneticField magModel) {
        this.forceModels = forceModels;
        this.initialDate = initialDate;
        this.earth = earth;
        this.magModel = magModel;
        this.itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, true);
    }

    @Override
    public void handleStep(final SpacecraftState currentState) {

        // Time since start [s]
        final double time = currentState.getDate().durationFrom(this.initialDate);

        // Position in the state's frame (likely inertial)
        final Vector3D pos = currentState.getPosition();

        // Allocate row
        final double[] stepData = new double[N_COLS];

        // Fill time + position
        stepData[0] = time;
        stepData[1] = pos.getX();
        stepData[2] = pos.getY();
        stepData[3] = pos.getZ();

        // --- Acceleration per force model ---
        // Forces occupy indices [4..18]. We do NOT write beyond 18 to protect B columns.
        for (int i = 0; i < this.forceModels.size(); i++) {

            final int baseIndex = 4 + (i * 3);

            // Only write if it fits inside the dedicated forces block.
            if (baseIndex + 2 < FORCES_BLOCK_END_EXCLUSIVE) {
                final ForceModel force = this.forceModels.get(i);
                final Vector3D a = force.acceleration(currentState, force.getParameters());

                stepData[baseIndex]     = a.getX();
                stepData[baseIndex + 1] = a.getY();
                stepData[baseIndex + 2] = a.getZ();
            } else {
                // If you have more than 5 forces, extra forces are ignored here.
                // (Alternatively, you could resize N_COLS and store them.)
                break;
            }
        }

        // --- Geomagnetic field BN, BE, BD (Tesla) ---
        try {
            final AbsoluteDate date = currentState.getDate();

            // Transform position to ITRF
            final Transform t = currentState.getFrame().getTransformTo(itrf, date);
            final Vector3D rITRF = t.transformPosition(currentState.getPVCoordinates().getPosition());

            // Convert to geodetic lat/lon/alt
            final GeodeticPoint gp = earth.transform(rITRF, itrf, date);
            final double lat = gp.getLatitude();   // rad
            final double lon = gp.getLongitude();  // rad
            final double h   = gp.getAltitude();   // m

            // Compute geomagnetic field (NED components, Tesla)
            final GeoMagneticElements e = magModel.calculateField(lat, lon, h);
            final Vector3D Bned = e.getFieldVector(); // (BN, BE, BD)
            
            // Build NED unit vectors expressed in ITRF
            double sLat = Math.sin(lat), cLat = Math.cos(lat);
            double sLon = Math.sin(lon), cLon = Math.cos(lon);

            Vector3D eE = new Vector3D(-sLon,  cLon, 0.0);             // East
            Vector3D eN = new Vector3D(-sLat*cLon, -sLat*sLon, cLat);  // North
            Vector3D eD = new Vector3D(-cLat*cLon, -cLat*sLon, -sLat); // Down (= -Up)
            
            // Convert NED to ITRF
            double BN = Bned.getX();
            double BE = Bned.getY();
            double BD = Bned.getZ();
            
            Vector3D Bitrf = eN.scalarMultiply(BN)
                        .add(eE.scalarMultiply(BE))
                        .add(eD.scalarMultiply(BD));
            
            // Convert ITRF -> GCRF (vector transform, NOT position)
            Frame gcrf = FramesFactory.getGCRF();
            Transform itrfToGcrf = itrf.getTransformTo(gcrf, date);
            Vector3D Bgcrf = itrfToGcrf.transformVector(Bitrf);
            
            // Store inertial B components (Tesla) in history columns 20–22
            stepData[19] = Bgcrf.getX(); // B_X_GCRF
            stepData[20] = Bgcrf.getY(); // B_Y_GCRF
            stepData[21] = Bgcrf.getZ(); // B_Z_GCRF

            // If you prefer nanoTesla:
            // stepData[19] = Bned.getX() * 1e9;
            // stepData[20] = Bned.getY() * 1e9;
            // stepData[21] = Bned.getZ() * 1e9;

        } catch (Exception ex) {
            stepData[19] = Double.NaN;
            stepData[20] = Double.NaN;
            stepData[21] = Double.NaN;
        }

        // Store row
        history.add(stepData);
    }

    /**
     * Retrieves the complete simulation history as a 2D array.
     *
     * @return A 2D double array where each row is one time step.
     */
    public double[][] getHistory() {
        final double[][] historyArray = new double[history.size()][N_COLS];
        for (int i = 0; i < history.size(); i++) {
            historyArray[i] = history.get(i);
        }
        return historyArray;
    }
}