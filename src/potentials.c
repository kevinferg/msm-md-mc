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

/* Oscillating Pair Potential */

int get_two_powers(double r, int p1, int p2, double* A1, double* A2) {
	// Takes input r to 2 different integer powers
	// A1: r^(p1)
	// A2: r^(p2) 
	*A1 = 1.;
	*A2 = 1.;
	while (p1 || p2) {
		if (p1 % 2) {
			*A1 *= r;
			printf("A1 * %f --> %f\n",r,*A1);
		}
		if (p2 % 2) {
			*A2 *= r;
			printf("A2 * %f --> %f\n",r,*A2);
		}
		p2 /= 2;
		p1 /= 2;
		if (!p1 && !p2) {
			break;
		}
		r *= r;
	}
	return 0;
}

double U_opp(double r, const void* params) {
	const OPPParams* P = (const OPPParams*) params;
	double A1, A2;
	get_two_powers(1./r, P->eta1, P->eta2, &A1, &A2);
	return P->C1 * A1 + P->C2 * A2 * cos(P->k * r - P->phi);
}

double F_opp(double r, const void* params) {
	const OPPParams* P = (const OPPParams*) params;
	double A1, A2;
	get_two_powers(1./r, P->eta1 - 1, P->eta2 - 1, &A1, &A2);
	double term1 = P->eta1 * P->C1 * A1;
	double term2 = P->eta2 * P->C2 * A2 * cos(P->k * r - P->phi);
	double term3 = P->k * P->eta2 * P->C2 * A2 * r * sin(P->k * r - P->phi);;
	return -(term1 + term2 + term3);
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