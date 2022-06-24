#ifndef _TYPES_H_
#define _TYPES_H_

#include <stdio.h>

typedef unsigned int uint;


/* Vector with 3 components
Can refer to components as xyz or as indices
For example: vec.x == vec.V[0] */
typedef union{
	struct{double x,y,z;};
	double V[3];
}Vec3;


/* One particle with position, velocity, current acting force, and mass*/
typedef struct{
	Vec3 pos;
	Vec3 vel;
	Vec3 force;
	Vec3 trupos;
	double mass;
}Particle;


/* Function pointer that takes in a distance and an array of parameters as arguments, for potential/force calculations*/
typedef double (*dist_fun) (double, double*);


/* A potential, contains the functions used for force/potential calculations, 
as well as the parameters needed, like epsilon and sigma in Lennard Jones*/
typedef struct{
	dist_fun U_fun;
	dist_fun F_fun;
	
	double r_cut;
	double U_cut;
	double F_cut;
	
	double* params;
}Potential;



/* Contains all information about an MD system 
Particle list, number of particles, whether particles have been set up;
Pair list, number of pairs, whether pairs have been set up;
The potential to be used for every pair of particles;
Number of time steps run so far;
Size of each time step, in dimensionless time units;
*/
typedef struct{
	long time_steps;
	double dt;
	double runtime;
	
	Vec3 boxlen;
	Vec3 xyz_min;
	Vec3 xyz_max;
	double virial;
	
	int N_particles;
	Potential* potential;
	Particle* particles;
	int particles_initialized;
	
	int anim_initialized;
	int anim_every;
	FILE *anim_file;
	
	int log_initialized;
	int log_every;
	FILE *log_file;
	
	int nvt_on;
	double zeta;
	double tau_damp;
	double T_des;
	
}MDSystem;



#endif