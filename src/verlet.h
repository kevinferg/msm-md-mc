#ifndef VERLET_H
#define VERLET_H

#include "types.h"

int vv_step(MDSystem* sys);

int nve_step(MDSystem* sys);
int nvt_step(MDSystem* sys);

int nve_update_velocities(MDSystem* sys);
int nvt_update_velocities(MDSystem *sys,int step);
int nvt_update_zeta(MDSystem *sys, double T);

int vv_update_positions(MDSystem* sys);
int vv_update_forces(MDSystem* sys);
int vv_update_velocities(MDSystem* sys);

/************************************/
#endif