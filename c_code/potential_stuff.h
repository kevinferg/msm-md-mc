#ifndef _POTENTIAL_STUFF_H_
#define _POTENTIAL_STUFF_H_

#include "types.h"


int potential_init(Potential* potential, double r_cut, dist_fun pot_fun, dist_fun frc_fun,double* params);
double F_lj(double r,double* params);
double U_lj(double r,double* params);

double F_harmonic(double r,double* params);
double U_harmonic(double r,double* params);

double F_get(Potential* potential,double r);
double U_get(Potential* potential,double r);


#endif