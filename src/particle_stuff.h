#ifndef _PARTICLE_STUFF_H_
#define _PARTICLE_STUFF_H_

#include "types.h"

int pt_set_force(Particle *pt,Vec3 *vec);
int pt_set_pos(Particle *pt,Vec3 *vec);
int pt_set_vel(Particle *pt,Vec3 *vec);
int pt_set_mass(Particle *pt,double val);

int pt_get_force(Particle *pt,Vec3 *vec);
int pt_get_pos(Particle *pt,Vec3 *vec);
int pt_get_vel(Particle *pt,Vec3 *vec);
int pt_get_mass(Particle *pt,double *val);

int pt_add_force(Particle *pt,Vec3 *vec);
int pt_add_pos(Particle *pt,Vec3 *vec);
int pt_add_trupos(Particle *pt,Vec3 *vec);
int pt_add_vel(Particle *pt,Vec3 *vec);
#endif
