#ifndef UTILS_H
#define UTILS_H

#include "types.h"

/* For a system 'sys' with a potential acting between all pairs of particles,
   this function calculates the energy and force between a pair separated by
   r, for r between 0 and 5.
   The results are tabulated in 'filename' */
int check_potential(MDSystem *sys, const char* filename);


/* This function creates an MD system, runs it in the NVT ensemble at a 
   specified temperature, then runs it in NVE, and then deletes the
   system.
   Particle are initialized to the locations in 'input_file'.
   Material property information is exported to 'log_file'.
   A trajectory is exported to 'traj_file'*/
int md_simulation(const char* input_file, const char* log_file, const char* traj_file, 
                  int log_every, int traj_every, double time_step,
                  double boxlen, double temp, double init_speed, double tau,
                  long nvt_steps, long nve_steps, int progress, unsigned int seed);


/* This function creates a system and runs an MC simulation at a specified
   temperature.
   Particle are initialized to the locations in 'input_file'. 
   Material property information is exported to 'log_file'.
   A snapshot of the final system is exported to 'snapshot_file'.
   */
int mc_simulation(const char* input_file, const char* log_file, const char* snapshot_file, 
                  int log_every, double boxlen, double temp, double max_trial_move,
                  long num_steps, int verbose, unsigned int seed);

#endif