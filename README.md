# Molecular Simulation of Materials: <br/> Basic Molecular Dynamics and Monte Carlo Simulations

## Project Description
Molecular simulation code written from scratch in C as coursework for [24-623: Molecular Simulation of Materials](https://www.meche.engineering.cmu.edu/education/courses/24-623.html).

Provides functionality for running MD simulations in NVE or NVT, complete with material property logging, user-defined potentials, particle trajectory output, and more. MC simulation is also included.


The units for all quantities are dimensionless "Lennard-Jones" units, each derived in terms of the base units of length, mass, and energy. [HOOMD-Blue's documentation](https://hoomd-blue.readthedocs.io/en/stable/units.html) has some examples of this.

## Features

#### Simulation
- NVE (microcanonical) Molecular Dynamics simulations with Velocity-Verlet integration
- NVT (canonical) Molecular Dynamics simulations using the Nos&#x00E9;-Hoover thermostat
- Markov Chain Monte Carlo simulations using the [Rosenbluth-Metropolis-Hastings Algorithm](https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm)

#### Measurement
- Measure the following properties of the system instantaneously:
  - Temperature and pressure
  - Volume and density
  - Kinetic energy, potential energy, and Hamiltonian
  - Mean squared displacement (MSD)
  - x-, y-, and z-components of momentum and center of mass

#### File I/O
- Read particle locations from a .txt file
- Output a text file log of material properties over time during a simulation
- Output particle locations as a .xyz or .pdb file
- Output particle trajectories as a .xyz animation file


#### Potentials
- Lennard-Jones potential with variable $\ \varepsilon $ and $\ \sigma $
- Harmonic bond (bonding is unfinished)
- Potentials from user-specified functions


#### Other features
- Built-in cutoff radius for potentials
- Periodic boundary conditions (with optionally separate components)
- Randomize particle velocities to match a target value
- Set number density of the system
- Basic live system rendering on the terminal during a simulation


## Usage
The code can be compiled on Windows using `mingw32-make`.
Run `test.bat` to both compile and run the executable.



## Examples

### MD Simulation
This snippet demonstrates how to create an NVT simulation. To measure properties, it is recommended to let the system equilibrate in NVT, then turn off the thermostat to run in NVE, and log system properties there, as seen in this example.
```
MDSystem sys;                            // Create system
sys_init(&sys);                          // Initialize system with default Lennard-Jones potential

io_load_txt(&sys, "liquid256.txt");      // Load particles locations from file "liquid256.txt"
sys_set_boxlen(&sys, 6.8);               // Periodic boundary with box side lengths 6.8

anim_init(&sys, 100, "trajectory.xyz");  // Output locations to "trajectory.xyz" every 100 steps
log_init(&sys, 10, "log.txt");           // Output material properties to "log.txt" every 10 steps

sys_random_velocities(&sys, 1.773);      // Randomize particle speeds to about 1.773 on average 
                                         // (but keep zero overall momentum)

sys_set_dt(&sys, 0.002);                 // Set time step to 0.002 time units
sys_nvt_ensemble(&sys, temp, 0.05);      // Set system ensemble to NVT with tau = 0.05 time units
sys_run(&sys, 50000, 50);                // Run in NVT for 50,000 steps; print progress every 50 steps

sys_zero_trupos(&sys);                   // Reset mean-squared-displacement particle locations
sys_nve_ensemble(&sys);                  // Set system ensemble to NVE
sys_run(&sys, 50000, 50);                // Run in NVE for 50,000 steps; print progress every 50 steps

sys_destroy(&sys);                       // Free the system and close log/trajectory files
```

### MC Simulation
This example shows a basic MC simulation.
```
MDSystem sys;                               // Create system
sys_init(&sys);                             // Initialize system with default Lennard-Jones potential

io_load_txt(&sys,"liquid256.txt");          // Load particles locations from file "liquid256.txt"
sys_set_boxlen(&sys, 6.8);                  // Periodic boundary with box side lengths 6.8
log_init(&sys, 1, "log.txt");               // Output material properties to "log.txt" every 1 step

sys_run_mc(&sys, 0.831716, 0.1, 2.5e6, 0);  // Run a Monte Carlo simulation...
                                            //   Dimensionless temperature:  kT = 0.831716
                                            //   Max particle perturbation:  dr = 0.1
                                            //   2,500,000 steps
                                            //   0 --> Do not print progress
                                            
io_export_xyz(&sys, "final_snapshot.xyz");  // Output final particle locations
sys_destroy(&sys);                          // Free the system and close log file
```