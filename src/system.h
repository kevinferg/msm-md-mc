#ifndef SYSTEM_H
#define SYSTEM_H

#define INITIALIZED 1417

#include "types.h"


int sys_zero_forces(MDSystem *sys);
int sys_unit_masses(MDSystem *sys);
int sys_zero_positions(MDSystem* sys);
int sys_zero_velocities(MDSystem* sys);
int sys_zero_trupos(MDSystem* sys);
int sys_zero_all(MDSystem* sys);
int sys_zero_momentum(MDSystem* sys);
int sys_random_velocities(MDSystem* sys,double A);



int sys_set_N_particles(MDSystem* sys,int N);

int sys_init(MDSystem* sys);
int sys_set_dt(MDSystem* sys, double dt);

int sys_nvt_ensemble(MDSystem *sys, double T_des, double tau);
int sys_nve_ensemble(MDSystem *sys);

int sys_run(MDSystem* sys, long numsteps,int progress);
int sys_destroy(MDSystem* sys);

int sys_print_positions(MDSystem *sys);
int sys_set_boxlen(MDSystem *sys,double boxlen);

double sys_get_runtime(MDSystem* sys);
int sys_get_N_particles(MDSystem* sys);

int sys_rescale(MDSystem* sys, Vec3 *new_L);
int sys_set_density(MDSystem* sys, double density);

int sys_run_mc(MDSystem* sys, double kT, double maxstep,unsigned long numsteps, int verbose);

/**************************************************/
#endif