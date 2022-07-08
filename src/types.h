#ifndef TYPES_H
#define TYPES_H

#include <stdio.h>

typedef unsigned int uint;

/* Vector with 3 components
Can refer to components as xyz or as indices
For example: vec.x == vec.V[0] */
typedef union {
	struct{double x,y,z;};
	double V[3];
} Vec3;


/* One particle with position, velocity, current acting force, and mass*/
typedef struct {
	Vec3 pos;       // Position of particle
	Vec3 vel;       // Velocity of particle
	Vec3 force;     // Current force acting on particle
	Vec3 trupos;    // Position w.r.t. original location, ignoring PBC's
	double mass;    // Mass of particle
} Particle;


/* Function pointer that takes in a distance and a struct of parameters as arguments, for potential/force calculations*/
typedef double (*dist_fun) (double, const void*);


/* A potential, contains the functions used for force/potential calculations, 
as well as the parameters needed, like epsilon and sigma in Lennard Jones*/
typedef struct {
	dist_fun U_fun;    // Pointer to potential (U) function
	dist_fun F_fun;    // Pointer to force (F = -dU/dr) function
	
	double r_cut;      // User-specified cutoff radius
	double U_cut;      // Internally-generated U at cutoff
	double F_cut;      // Internally-generated F at cutoff
	
	double* params;    // Pointer to an array of potential parameters.
} Potential;



/* Contains all information about an MD system 
Particle list, number of particles, whether particles have been set up;
The potential to be used for every pair of particles;
Number of time steps run so far;
Size of each time step, in dimensionless time units;
*/
typedef struct {
	long time_steps;     // Number of time steps run so far 
	double dt;           // Time step in dimensionless units
	double runtime;      // Cumulative time it took simulation to run
	
	Vec3 boxlen;         // xyz-components of box side length
	Vec3 xyz_min;        // Lower bounds of box
	Vec3 xyz_max;        // Upper bounds of box
	double virial;       // Current virial component of pressure
	
	int N_particles;            // Number of particles
	Potential* potential;       // Pointer to potential acting between all particle pairs
	Particle* particles;        // Pointer to particle array
	int particles_initialized;  // Whether particles are initialized
	
	int anim_initialized;    // Whether trajectory animation file has been initialized
	int anim_every;          // Add particle snapshot to trajectory every _ frames
	FILE *anim_file;         // Pointer to trajectory animation file
	
	int log_initialized;     // Whether property log file has been initialized
	int log_every;           // Add system properties to log every _ frames
	FILE *log_file;          // Pointer to log file
	
	int nvt_on;          // Whether to use an NVT thermostat 
	double zeta;         // Nose-Hoover thermodynamic friction coefficient
	double tau_damp;     // Nose-Hoover dimensionless damping time constant
	double T_des;        // Desired system temperature for NVT
	
} MDSystem;

#endif