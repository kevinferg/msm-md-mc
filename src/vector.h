#ifndef VECTOR_H
#define VECTOR_H

#include "types.h"

extern Vec3 ZERO_VEC;

int vec_copy(Vec3 *src, Vec3 *dest);
int vec_add(Vec3 *src, Vec3 *dest);
int vec_times_scalar(Vec3 *vec, double scalar);
int vec_get_diff(Vec3 *vecdiff, Vec3 *vec1, Vec3 *vec2);
double vec_get_mag(Vec3 *vec);
int vec_print(Vec3 *vec);
double vec_dot(Vec3 *v1, Vec3 *v2);

#endif