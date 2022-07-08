#ifndef POTENTIALS_H
#define POTENTIALS_H

#include "types.h"

/* Initialize a Potential

  potential: pointer to Potential
  r_cut:     cutoff radius
  pot_fun:   potential energy function
  frc_fun:   force function
  params:    struct of parameters

  Warning! 
  Make sure 'params' and 'potential' exist within the same scope as the simulation system.
  i.e. If they are created within a function, a potential initialized with them will no
       longer work once the function is returned, unless 'potential' or 'params' were
       static, or 'params' is an array literal. */
int potential_init(Potential* potential, double r_cut, dist_fun pot_fun, dist_fun frc_fun, const void* params);

/* Evaluates the energy for 'potential' at a given distance 'r', applying the cutoff radius */
double U_get(Potential* potential, double r);

/* Evaluates the force for 'potential' at a given distance 'r', applying the cutoff radius */
double F_get(Potential* potential, double r);


/***********************************  Potentials  *********************************************/

/* Lennard-Jones Potential */
typedef struct LJParams {
   double epsilon; // Energy
   double   sigma; // Distance
} LJParams;
double U_lj(double r, const void* params);
double F_lj(double r, const void* params);


/* Morse Potential */
typedef struct MorseParams {
   double D_e; // Well depth
   double   a; // Well width (larger a = wider well)
   double r_e; // Equilibrium distance
} MorseParams;
double U_Morse(double r, const void* params);
double F_Morse(double r, const void* params);


/* Harmonic Bond Potential */
typedef struct HarmonicParams {
   double  k; // Stiffness
   double L0; // Equilibrium distance
} HarmonicParams;
double U_harmonic(double r, const void* params);
double F_harmonic(double r, const void* params);





#endif