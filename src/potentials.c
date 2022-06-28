#include <stdio.h>
#include <math.h>
#include "types.h"
#include "potentials.h"

int potential_init(Potential* potential, double r_cut, dist_fun pot_fun, dist_fun frc_fun, double* params){
	potential->U_fun = pot_fun;
	potential->F_fun = frc_fun;
	potential->params = params;
	
	potential->r_cut = r_cut;
	potential->U_cut = (*pot_fun) (r_cut,params);
	potential->F_cut = (*frc_fun) (r_cut,params);
	
	return 0;
}

double F_get(Potential* potential,double r){

	return (r>potential->r_cut)? 0 :
		   (*(potential->F_fun)) (r, potential->params) - potential->F_cut;
	
}
double U_get(Potential* potential,double r){

	return (r>potential->r_cut)? 0 :
		(*(potential->U_fun)) (r, potential->params) - potential->U_cut + (r-potential->r_cut)*potential->F_cut;
}




/*******************************************  Potentials  *****************************************************/

/* Lennard-Jones Potential */
double U_lj(double r, double* params){
	double eps = params[0];
	double sig = params[1];
	double s_over_r = sig/r;
	double sr3 = s_over_r * s_over_r * s_over_r;
	double sr6 = sr3 * sr3;

	return 4.*eps*sr6*(sr6-1.);
}
double F_lj(double r, double* params){
	double eps = params[0];
	double sig = params[1];

	double r3 = r * r * r;
	double r6 = r3 * r3;
	double r12 = r6 * r6;

	double s3 = sig * sig * sig;
	double s6 = s3 * s3;

	return -24.*s6*eps*(r6-2.*s6)/r12/r;
}



/* Harmonic Bond Potential */
double U_harmonic(double r,double* params){
	double k = params[0];
	double L = params[1];
	r = r-L;
	double energy = r*r/2.*k;

	return energy;
}
double F_harmonic(double r,double* params){
	double k = params[0];
	double L = params[1];
	double force = -k*(r-L);

	return force;
}