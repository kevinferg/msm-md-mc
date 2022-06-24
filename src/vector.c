#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "types.h"
#include "vector.h"

Vec3 ZERO_VEC = {.x=0.,.y=0.,.z=0.};

int vec_copy(Vec3 *src, Vec3 *dest){
	dest->x = src->x;
	dest->y = src->y;
	dest->z = src->z;
	return 0;
}

int vec_add(Vec3 *src, Vec3 *dest){
	(dest->x) += (src->x);
	(dest->y) += (src->y);
	(dest->z) += (src->z);
	return 0;
}

int vec_times_scalar(Vec3 *vec, double scalar){
	(vec->x) *= scalar;
	(vec->y) *= scalar;
	(vec->z) *= scalar;
	return 0;
}

int vec_get_diff(Vec3 *vecdiff, Vec3 *vec1, Vec3 *vec2){
	vecdiff->x = vec1->x - vec2->x;
	vecdiff->y = vec1->y - vec2->y;
	vecdiff->z = vec1->z - vec2->z;
	return 0;
	
}

double vec_get_mag(Vec3 *vec){
	//return sqrt(pow(vec->x,2.)+pow(vec->y,2.)+pow(vec->z,2.));
	// ^ The above takes way longer than the below:
	return sqrt(vec->x * vec->x  +  vec->y * vec->y  +  vec->z * vec->z);
}

int vec_print(Vec3 *vec){
	printf("[%9.4f   %9.4f   %9.4f]\n",vec->x,vec->y,vec->z);
	return 0;
}

double vec_dot(Vec3 *v1, Vec3 *v2){
	return (v1->x)*(v2->x) + (v1->y)*(v2->y) + (v1->z)*(v2->z);
}



