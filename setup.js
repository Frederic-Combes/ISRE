function Settings() {
    var parent = this;

    this.dimension      = undefined;    // Space dimension, 1, 2 or 3

    this.duration       = undefined;    // Length of the simulation (in seconds)
    this.timesteps      = undefined;    // Number of timesteps
    this.saveInterval   = 0.02;         // Time between two saved points (in seconds)

    this.energySamples  = 1000;         // Number of energy samples for the simulation
    this.energyBins     = [];           // Used to save partial spins
}

// Operations (allowing the realization of pi-pulse and other transformations of the spins)
function Operation() {
    this.isFixedTime = false;
    this.time        = undefined;
    this.operator    = [[1,0,0],[0,1,0],[0,0,1]];
}

// Control parameters (allowing a simulation to be run for various parameters)
function Parameter(name, values) {
    this.name   = name;
    this.values = values instanceof Array ? values : [values];
}

function Function(name, code) {
    this.name = name;
    this.code = code;
}

// Simulation parameters
/*function Simulation() {
    this.frequency  = "0";
    this.exchange   = "0";
    this.damping    = "0";

    this.atomLoss = {
        fromF1 : "0",
        fromF2 : "0"
    };

    this.initialSpin = ["1", "0", "0"];

    this.operations = [];
    this.parameters = [];
    this.functions = [];
    this.settings = new Settings();

    this.name = undefined;
}*/

// Application adds to the script environment :
//  * simulate(settings, simulation, parameters)    // starts the simulation using the given settings, simulation and control parameters
//  * loadScript(filename)                          // loads the given script /!\ the current directory is the directory where the ISRE command is executed

// Function expressions are declared as string (ex : simulation.frequency = "sin(energy) * exp(-time)")
// The application reserves the following names (they must not be used as function names or as parameters names) :
//  * timestep,         // Contains the current timestep of the simulation
//  * time,             // Contains the current time of the simulation
//  * energyIndex,      // (Caching purpose only)
//  * energy,           // Contains the current energy value
//  * duration,         // Contains the duration (end time) of the simulation (matches settings.duration if all operations are fixed-time)
//  * frequency,        // Contains the value of the frequency (for a given time and energy)
//  * exchange,         // Contains the value of the exchange frequency (for a given time and energy)
//  * damping,          // Contains the value of the damping (for a given time and energy)
//  * atomLossFromF1,   // Contains the value of the atom loss from the state F=1,mf=0 to F=1,mf!=0
//  * atomLossFromF2,   // Contains the value of the atom loss from the state F=2,mf=0 to F=2,mf!=0
//  * initialSpinX      // Initial spin in the X direction
//  * initialSpinY      // Initial spin in the Y direction
//  * initialSpinZ      // Initial spin in the Z direction
// Available built-in functions are :
//  * the usual mathematical functions : cos, sin, log, exp, ...
//  * DefaultDensityEnergyDependence    // The default density energy dependence (exp(-var/2) * BesselI(0, func(dimension) * var)

// Example:
//      lifetime            = new Parameter("lifetime", [1, 2, 4, Math.Inf])
//      baseFrequency       = new Parameter("c_param_0", 1)
//      DLSInhomogeneity    = new Parameter("DLS", [1,2,3])
//      MFInhomogeneity     = new Parameter("MF", [4,5])
//      densityInitial      = new Parameter("n_0", [1,2,3,4,5]);
//
//      functionDensityTimePart     = new Function("n_t", "exp(-time / lifetime)");
//      functionDensityEnergyPart   = new Function("n_E", "DefaultDensityEnergyDependence(energy)");
//      functionDensity             = new Function("density", "n_0 * n_t * n_E");
//
//      sim = new Simulation()
//      sim.frequency   = "c_param_0 + DLS * energy + MF * density"
//      sim.exchange    = "1"
//      ...
//      sim.settings = globalSettings
//
//      sim.run([lifetime, baseFrequency, DLSInhomogeneity, MFInhomogeneity, densityInitial])

