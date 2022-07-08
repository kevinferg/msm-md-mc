#include <stdio.h>
#include <math.h>
#include "types.h"
#include "potentials.h"

int potential_init(Potential* potential, double r_cut, dist_fun pot_fun, dist_fun frc_fun, void* params){
	potential->U_fun = pot_fun;
	potential->F_fun = frc_fun;
	potential->params = params;
	
	potential->r_cut = r_cut;
	potential->U_cut = (*pot_fun) (r_cut, params);
	potential->F_cut = (*frc_fun) (r_cut, params);
	
	return 0;
}

double F_get(Potential* potential, double r) {

	return (r>potential->r_cut)? 0 :
		   (*(potential->F_fun)) (r, potential->params) - potential->F_cut;
	
}
double U_get(Potential* potential, double r) {

	return (r>potential->r_cut)? 0 :
		(*(potential->U_fun)) (r, potential->params) - potential->U_cut + (r-potential->r_cut)*potential->F_cut;
}




/*******************************************  Potentials  *****************************************************/

/* Lennard-Jones Potential */
double U_lj(double r, const void* params) {
	const LJParams* P = (const LJParams*) params;

	double s_over_r = P->sigma / r;
	double sr3 = s_over_r * s_over_r * s_over_r;
	double sr6 = sr3 * sr3;

	return 4. * P->epsilon * sr6 * (sr6-1.);
}
double F_lj(double r, const void* params) {
	const LJParams* P = (const LJParams*) params;

	double r3 = r * r * r;
	double r6 = r3 * r3;
	double r12 = r6 * r6;

	double s3 = P->sigma * P->sigma * P->sigma;
	double s6 = s3 * s3;

	return -24.*s6 * P->epsilon * (r6-2.*s6) / r12 / r;
}

/* Morse Potential */
double U_Morse(double r, const void* params) {
	const MorseParams* P = (const MorseParams*) params;
    double part = 1 - exp(P->a * (P->r_e - r));
    return P->D_e * part * part;
}
double F_Morse(double r, const void* params) {
    const MorseParams* P = (const MorseParams*) params;
    double part = exp(P->a * (P->r_e - r));
    return -2 * P->a * P->D_e * part * (1 - part);
}

/* Harmonic Bond Potential */
double U_harmonic(double r, const void* params) {
	const HarmonicParams* P = (const HarmonicParams*) params;
	r = r - P->L0;
	return r * r/2. * P->k;
}
double F_harmonic(double r, const void* params) {
	const HarmonicParams* P = (const HarmonicParams*) params;
	return -P->k * (r - P->L0);
}