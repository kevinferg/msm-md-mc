#include <stdio.h>

#include "types.h"
#include "vector.h"
#include "particle.h"


int pt_set_force(Particle *pt,Vec3 *vec){
	vec_copy(vec, &(pt->force));
	return(0);
}

int pt_set_pos(Particle *pt,Vec3 *vec){
	vec_copy(vec,&(pt->pos));
	return(0);
}

int pt_set_vel(Particle *pt,Vec3 *vec){
	vec_copy(vec, &(pt->vel));
	return(0);
}
int pt_set_mass(Particle *pt,double val){
	pt->mass = val;
	return(0);
}

int pt_get_force(Particle *pt,Vec3 *vec){
	vec_copy(&(pt->force),vec);
	return(0);
}

int pt_get_pos(Particle *pt,Vec3 *vec){
	vec_copy(&(pt->pos),vec);
	return(0);
}

int pt_get_vel(Particle *pt,Vec3 *vec){
	vec_copy(&(pt->vel),vec);
	return(0);
}

int pt_get_mass(Particle *pt,double *val){
	*val = pt->mass;
	return(0);
}

int pt_add_force(Particle *pt,Vec3 *vec){
	vec_add(vec, &(pt->force));
	return(0);
}

int pt_add_pos(Particle *pt,Vec3 *vec){
	vec_add(vec, &(pt->pos));
	return(0);
}

int pt_add_trupos(Particle *pt,Vec3 *vec){
	vec_add(vec, &(pt->trupos));
	return(0);
}


int pt_add_vel(Particle *pt,Vec3 *vec){
	vec_add(vec, &(pt->vel));
	return(0);
}

