#ifndef UTILS_H
#define UTILS_H

#include "types.h"

int check_potential(MDSystem *sys, const char* filename);

int md_simulation(const char* input_file, const char* log_file, const char* traj_file, 
                  int log_every, int traj_every, double time_step,
                  double boxlen, double temp, double init_speed, double tau,
                  long nvt_steps, long nve_steps, int progress, unsigned int seed);

int mc_simulation(const char* input_file, const char* log_file, const char* snapshot_file, 
                  int log_every, double boxlen, double temp, double max_trial_move,
                  long num_steps, int verbose, unsigned int seed);

#endif