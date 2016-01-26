simulation.dimension    = 1;                                                // Space dimension, 1,2 or 3                                        // Must be 1,2 or 3                                                         // Defaults to 1

simulation.time.length  = 2;                                                // Length of the simulation in sec                                  // Must be a single value                                                   // Defaults to 2 sec
simulation.time.steps   = 2000;                                             // Number of time steps                                             // Must be a single value                                                   // Defaults to 2000
simulation.time.points  = 100;                                              // Number of time points to save                                    // Must be a single value, must divide simulation.time.steps                // Defaults to 100

simulation.energy.count = 1000;                                             // Number of energy samples                                         // Must be a single value                                                   // Defaults to 1000
simulation.energy.bins  = [[0,Math.log(2)],[Math.log(2),Infinity]];         // Energy bins                                                      // Must be an array of the form [[bin1min,bin2max],[bin2min,bin2max],...]   // Defaults to [] (no energy bining)

simulation.spinEcho.enabled = true;                                         // Spin echo                                                        // Must be true of false                                                    // Defaults to false (no spin echo)
simulation.spinEcho.time    = [0.25,0.75];                                  // Spin echo time as a fraction of the simulation time              // Must be a value or a list of values between 0 and 1                      // Defaults to [0.25,0.75] (Spin echo happening at 1/4 and 3/4 of the simulation time)

simulation.frequency.larmor         = 0;                                    // Larmor frequency in rad / s.                                     // Can be a single value or a list      // Defaults to 0
simulation.frequency.exchange       = 1;                                    // Exchange frequency in rad / s                                    // Can be a single value or a list      // Defaults to 0
simulation.inhomogeneity.larmor     = simulation.makeList(0,5,6);           // Inhomogeneity in larmor frequency inside the trap in rad /sec    // Can be a single value or a list      // Defaults to 0
simulation.inhomogeneity.density    = [0,1,2];                              // Inhomogeneity related to density                                 // Can be a single value or a list      // Defaults to 0
simulation.relaxation               = 0;                                    // Relaxation rate in rad / s                                       // Can be a single value or a list      // Defaults to 0

simulation.density.value            = [1];                                  // Prefactor that multiplies all quantities proportional to the density // Can be a single value or a list  // Defaults to [1]

var tau = 6;
// Use either one of these
simulation.density.lifetime         = Infinity;                             // Lifetime in s. If simulation.density.timeDependence is not specified, the simulation.density.timeDependencen  = function(t) {return Math.exp(-t/simulation.density.lifetime);}
simulation.density.timeDependence   = function(t) {return Math.exp(-t/tau);}// Time dependence of the density. During the simulation, at time t (in seconds), the quantities related to the density are scaled by timeDependence(t)

simulation.name       = undefined                                           // name of the simulation                                         // Must be a string                     // Defaults to undefined
