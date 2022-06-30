#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "types.h"
#include "vector.h"
#include "calculations.h"
#include "potentials.h"
#include "system.h"
#include "particle.h"

double wrap_nearest(double number, double reference, double range) {
	double alt = (number<reference)? (number+range): (number-range);
	return (fabs(number-reference) <= fabs(alt-reference))? number: alt;
}

double wrap_within(double number, double min, double max) {
	double range = max - min;
	return number + range*(number < min) - range*(number > max);
}

int calc_get_rvec(MDSystem *sys, Vec3* rvec, Particle* p1, Particle* p2) {
	/* Get vector from p1 to p2, satsifying the PBCs in sys */
	
	register double x = wrap_nearest(p2->pos.x,p1->pos.x,sys->boxlen.x);
	register double y = wrap_nearest(p2->pos.y,p1->pos.y,sys->boxlen.y);
	register double z = wrap_nearest(p2->pos.z,p1->pos.z,sys->boxlen.z);

	Vec3 p2_nearest = {.x=x, .y=y, .z=z};

	vec_get_diff(rvec, &p2_nearest, &(p1->pos));
	
	return 0;
}

int calc_correct_position(MDSystem *sys, Particle* p) {
	register double x = wrap_within(p->pos.x, sys->xyz_min.x, sys->xyz_max.x);
	register double y = wrap_within(p->pos.y, sys->xyz_min.y, sys->xyz_max.y);
	register double z = wrap_within(p->pos.z, sys->xyz_min.z, sys->xyz_max.z);
	Vec3 newpos = {.x=x, .y=y, .z=z};
	pt_set_pos(p, &newpos);
	return 0;
}

/*****************************************************************************/

int calc_system_com(MDSystem *sys,Vec3 *center) {
	int i;
	Vec3 com, pm;
	vec_copy(&ZERO_VEC, &com);
	
	for (i = 0; i < sys->N_particles; i++) {
		vec_copy(&(sys->particles[i].pos),&pm);
		vec_times_scalar(&pm,sys->particles[i].mass);
		vec_add(&pm,&com);
	}
	vec_times_scalar(&com,1./((double) sys->N_particles));
	vec_copy(&com,center);
	return 0;
}


double calc_single_potential(MDSystem *sys, int id, Particle* p1) {
	Particle* p2;
	Vec3 rvec;
	double r;
	double U=0;
	int j;
	
	for (j = 0; j < sys->N_particles; j++) {
		if (j != id){
			p2 = &((sys->particles)[j]);

			calc_get_rvec(sys,&rvec,p1,p2);
			r = vec_get_mag(&rvec);
			U += U_get(sys->potential,r);

		}
	}


	return U;
}

double calc_system_potential(MDSystem *sys){
	int i,j;
	Particle *p1, *p2;
	double r, U, U_tot = 0.;
	Vec3 a_to_b;

	for (i = 0; i < (sys->N_particles-1); i++) {
		for (j=i+1; j<sys->N_particles; j++) {
			
			p1 = &((sys->particles)[i]);
			p2 = &((sys->particles)[j]);

			calc_get_rvec(sys, &a_to_b, p1, p2);
			r = vec_get_mag(&a_to_b);
			U = U_get(sys->potential, r);
			U_tot += U;
		}
	}
	return U_tot;
}

double calc_system_kinetic(MDSystem *sys){
	int i;
	Particle *p;
	double mass, speed, E, E_tot=0.;
	Vec3 vel;
	for (i = 0; i < sys->N_particles; i++) {
		p = &(sys->particles[i]);
		pt_get_vel(p, &vel);
		pt_get_mass(p, &mass);
		speed = vec_get_mag(&vel);
		E = 0.5*speed*speed*mass;
		E_tot += E;
	}
	return E_tot;
}


double calc_single_virial(MDSystem *sys,int id){
	int j;
	Particle *p1, *p2;
	double r, F, virial=0.;
	Vec3 rvec;
	
	p1 = &((sys->particles)[id]);
	for (j = 0; j < sys->N_particles; j++) {
		if (id != j) {
			p2 = &((sys->particles)[j]);
			calc_get_rvec(sys,&rvec,p1,p2);
			r = vec_get_mag(&rvec); 
			F = F_get(sys->potential,r); 
			virial += F*r;
		}
		
	}
	
	return virial;
}





double calc_system_virial(MDSystem *sys){
	int i, j;
	Particle *p1, *p2;
	double r, F, virial=0.;
	Vec3 rvec;
	for (i = 0; i < (sys->N_particles-1); i++) {
		for (j = i+1; j < sys->N_particles; j++) {

			p1 = &((sys->particles)[i]);
			p2 = &((sys->particles)[j]);

			calc_get_rvec(sys, &rvec, p1, p2);// rvec: Vector from 1 to 2
			r = vec_get_mag(&rvec); 
			F = F_get(sys->potential, r); 
			virial += F*r;
		}
	}
	return virial;
}



double calc_system_temperature(MDSystem *sys, double KE) {
	double kT = 2. * KE / (3. * ((double) sys->N_particles - 1.));
	return kT;
}

double calc_system_pressure(MDSystem *sys, double kT, double virial) {
	double P, V, ideal_gas_component, virial_component;
	V = calc_system_volume(sys);
	
	ideal_gas_component = ((double) sys->N_particles) * kT / V;
	virial_component = (1./(3. * V)) * virial;
	
	P = ideal_gas_component + virial_component;
	return P;
}

int calc_system_momentum(MDSystem *sys,Vec3 *mv) {
	int i;
	Vec3 mv_total;
	Vec3 mv_single;
	vec_copy(&ZERO_VEC, &mv_total);
	
	for (i = 0; i < sys->N_particles; i++){
		vec_copy(&(sys->particles[i].vel), &mv_single);
		vec_times_scalar(&mv_single, sys->particles[i].mass);
		vec_add(&mv_single, &mv_total);
	}
	vec_times_scalar(&mv_total, 1./((double) sys->N_particles));
	vec_copy(&mv_total, mv);
	return 0;
}

double calc_system_volume(MDSystem* sys) {
	double V = (sys->boxlen).x * (sys->boxlen).y * (sys->boxlen).z;
	return V;
}


double calc_system_density(MDSystem* sys) {
	volatile double density = ((double) (sys->N_particles)) / calc_system_volume(sys);
	return density;
}

double calc_system_msd(MDSystem *sys) {
	int i;
	double msd = 0.;
	Vec3* pos;
	for (i = 0; i < sys->N_particles; i++) {
		pos = &((sys->particles)[i].trupos);
		msd += pos->x * pos->x + 
			   pos->y * pos->y +
			   pos->z * pos->z;
	}
	msd /= ((double) sys->N_particles);
	return msd;
}


