#ifndef POTENTIALS_H
#define POTENTIALS_H

#include "types.h"

/* Initialize a Potential

  potential: pointer to Potential
  r_cut:     cutoff radius
  pot_fun:   potential energy function
  frc_fun:   force function
  params:    array of parameters

  Warning! 
  Make sure 'params' and 'potential' exist within the same scope as the simulation system.
  i.e. If they are created within a function, a potential initialized with them will no
       longer work once the function is returned, unless 'potential' or 'params' were
       static, or 'params' is an array literal. */
int potential_init(Potential* potential, double r_cut, dist_fun pot_fun, dist_fun frc_fun, double* params);

/* Evaluates the energy for 'potential' at a given distance 'r', applying the cutoff radius */
double U_get(Potential* potential, double r);

/* Evaluates the force for 'potential' at a given distance 'r', applying the cutoff radius */
double F_get(Potential* potential, double r);


/***********************************  Potentials  *********************************************/

/* Lennard-Jones Potential
   params: 
   [0] = epsilon
   [1] = sigma */
double U_lj(double r, double* params);
double F_lj(double r, double* params);

/* Morse Potential
   params: 
   [0] = D_e, well depth
   [1] = a, well width (larger a = wider well)
   [2] = r_e, equilibrium distance */
double U_Morse(double r, double* params);
double F_Morse(double r, double* params);

/* Harmonic Bond Potential
   params: 
   [0] = k,  stiffness
   [1] = L0, equilibrium distance */
double U_harmonic(double r, double* params);
double F_harmonic(double r, double* params);





#endif