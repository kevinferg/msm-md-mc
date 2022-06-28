#ifndef VERLET_H
#define VERLET_H

#include "types.h"

/* Apply one velocity-verlet step to the system;
   Uses the currently-enabled thermodynamic ensemble */
int vv_step(MDSystem* sys);

int vv_update_positions(MDSystem* sys);
int vv_update_forces(MDSystem* sys);

#endif