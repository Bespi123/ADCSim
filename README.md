# ADCSim: Software for Attitude Determination and Control System Design and Simulation

## Overview
This file presents ADCSim, an open-source MATLAB software package for the design, simulation, and validation of satellite Attitude Determination and Control Systems (ADCS). The software integrates detailed models for satellite dynamics, high-fidelity sensors (gyroscopes, accelerometers, magnetometers, star trackers), and actuators (reaction wheels). A key feature is its integration with the Orekit space dynamics library for realistic orbital environment modeling.

ADCSim features an intuitive Graphical User Interface (GUI) for Model-in-the-Loop (MIL) simulations, and implements multiple estimation algorithms (EKF, UKF, Madgwick) and non-linear controllers. Its advanced functionalities include a robust algorithm to compute and visualize the 3D momentum envelope for reaction wheel arrays and the use of genetic algorithms to auto-tune gains for both the EKF and attitude controllers. By providing a comprehensive suite of tools, ADCSim offers a powerful, cost-effective platform for developing and testing satellite attitude control systems.

![ADCSim main window](src/figs/main_Guide.jpg)
Figure 1. ADCSim main window.

![ADCSim propagation results](src/figs/propagation_results.jpg)
Figure 2. Orbit propagation results obtained using `Orekit` simulator engine.

---

## Installation and Setup

### 1. Prerequisites
Before you begin, ensure you have the following software installed:
* **MATLAB**
* **Git** (for cloning the repository)
* **Java Development Kit (JDK)** (This is required to compile the `StepStorageHandler`).

### 2. Clone the Repository
Open your terminal and clone the repository:

```sh
git clone https://github.com/Bespi123/Satellite_ADCS_GUIDE
```
### 3. Project Structure
ADCSim follows the following structure:
```
Satellite_ADCS_GUIDE/
├── lib/
│   ├── orekit/
│   └── custom_java/
├── src/
│   ├── +adcsim/
│   ├── +core/
│   ├── +optimization/
│   ├── config/
│   ├── figs/
│   └── gui/
├── main.m
├── LICENCE.txt
├── .gitignore
└── README.md

```
### 4. Download Dependencies
For the software to work, you must place the required files in the correct folders.

1.  **Download Orekit:** Download `orekit-13.1.2.jar` and place it inside the `lib/orekit` folder.
    * **Link 1:** [https://www.orekit.org/static/downloads.html](https://www.orekit.org/static/downloads.html).
    * **Link 2:** [https://gitlab.orekit.org/orekit/orekit/-/releases](https://gitlab.orekit.org/orekit/orekit/-/releases)

3.  **Download Hipparchus:** Download all the Hipparchus-4.0.2 binary `.jar` files (e.g., `hipparchus-core-X.X.jar`, `hipparchus-filtering-X.X.jar`, etc.) and place them inside the `lib/orekit` folder.
    * **Link:** [https://hipparchus.org/downloads.html](https://hipparchus.org/downloads.html)

4.  **Download Orekit Data:** Download the `orekit-data.zip` file from the tutorials link.
    * **Link:** [https://www.orekit.org/static/downloads.html](https://www.orekit.org/static/downloads.html).
    * **Action:** Unzip the file and place the resulting `orekit-data` folder **inside the `lib/orekit/orekit-data` folder**.

4.  **Download IGRF-14 geomagnetic model:** Download the `IGRF14_Windows.COF` file from the link.
    * **Link:** [https://www.ncei.noaa.gov/products/international-geomagnetic-reference-field](https://www.ncei.noaa.gov/products/international-geomagnetic-reference-field)
    * **Action:** Unzip the file, rename it as `IGRF.COF`, and place the resulting file **inside the `lib/orekit/orekit-data` folder**.

6. **Compile the Java Handler:** Before running the main program, you must compile the custom Java handler.
    * Open MATLAB
    * Navigate to the `lib/custom_java` directory.
    * Run the `compile_java_handler.m` compile script in the MATLAB Command Window.
    * Verify that the `StepStorageHandler.class` file was created successfully.

Your final `/lib` folder structure must match this diagram:

```
Satellite_ADCS_GUIDE/
├── lib/
│   ├── orekit/
│   │   ├── orekit-13.1.2.jar
│   │   ├── hipparchus-core-4.0.2.jar
│   │   ├── ... (all other jars)
│   │   └── orekit-data/
│   │       ├── (data files and folders)
│   └── custom_java/
│       ├── compile_java_handler.m
│       ├── StepStorageHandler.class
│       └── StepStorageHandler.java
│
...
└── main.m
```

### Run the Application
Run the `main.m` file.

## Credits
Brayan Espinoza (bespi123)
