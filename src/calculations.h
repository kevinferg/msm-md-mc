#ifndef CALCULATIONS_H
#define CALCULATIONS_H

#include "types.h"

// Pair property calculations
int calc_get_rvec(MDSystem *sys, Vec3* rvec, Particle* p1, Particle* p2);
int calc_correct_position(MDSystem *sys, Particle* p);

double calc_single_potential(MDSystem *sys, int id, Particle* p1);
double calc_single_virial(MDSystem *sys, int id);

// System property calculations
double calc_system_potential(MDSystem *sys);
double calc_system_kinetic(MDSystem *sys);

double calc_system_virial(MDSystem *sys);
double calc_system_temperature(MDSystem *sys, double KE);
double calc_system_pressure(MDSystem *sys, double kT, double virial);

int calc_system_momentum(MDSystem *sys, Vec3 *mv);
int calc_system_com(MDSystem *sys, Vec3 *center);

double calc_system_volume(MDSystem* sys);
double calc_system_density(MDSystem* sys);
double calc_system_msd(MDSystem *sys);

#endif 