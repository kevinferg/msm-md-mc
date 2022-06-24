#include <stdio.h>
#include <math.h>

#include "calculations.h"
#include "types.h"
#include "vector.h"
#include "particle.h"
#include "system.h"
#include "potentials.h"
#include "verlet.h"

int vv_step(MDSystem* sys){
	if (sys->nvt_on){
		nvt_step(sys);
	} else {
		nve_step(sys);
	}
	
	return 0;
}

int nve_step(MDSystem* sys){
	nve_update_velocities(sys);
	vv_update_positions(sys);
	vv_update_forces(sys);
	nve_update_velocities(sys);
	return 0;
}


int nve_update_velocities(MDSystem *sys){
	/* Updates 1/2 step as in Velocity Verlet method */
	int i=0;
	double m,dt_over_2m;
	Vec3 update;
	Particle* p;
	for (i=0;i<sys->N_particles;i++){
		p = &(sys->particles[i]);
		pt_get_mass(p,&m);
		dt_over_2m = (sys->dt)/m/2.;
		vec_copy(&(p->force),&update);
		vec_times_scalar(&update,dt_over_2m);
		pt_add_vel(p,&update);
	}
	return 0;
}

int vv_update_positions(MDSystem *sys){
	int i=0;
	Vec3 update;
	Particle* p;
	for (i=0;i<sys->N_particles;i++){
		p = &(sys->particles[i]);
		vec_copy(&(p->vel),&update);
		vec_times_scalar(&update,(sys->dt));
		pt_add_pos(p,&update);
		pt_add_trupos(p,&update);
		calc_correct_position(sys,p);
	}
	return 0;
}

int vv_update_forces(MDSystem* sys){
	// Updates forces to the next time step...
	// Also updates virial component of pressure
	
	long i=0,j=0;
	Vec3 rvec,update;
	Particle *p1, *p2;
	double r,F;
	
	sys_zero_forces(sys);
	sys->virial = 0.;

	for (i=0;i<(sys->N_particles-1);i++){
		for (j=i+1;j<sys->N_particles;j++){

			p1 = &((sys->particles)[i]);
			p2 = &((sys->particles)[j]);

			calc_get_rvec(sys,&rvec,p1,p2);// rvec: Vector from 1 to 2

			vec_copy(&rvec,&update);
			r = vec_get_mag(&update); 
			F = F_get(sys->potential,r); 

			vec_times_scalar(&update,F/r); // update: Force on 2 from 1
			
			(sys->virial) += F*r;

			pt_add_force(p2,&update);
			vec_times_scalar(&update,-1.);

			pt_add_force(p1,&update);
			
		}
	}
	
	return 0;
}


int nvt_step(MDSystem* sys){
	double KE,T;
	nvt_update_velocities(sys,1);
	vv_update_positions(sys);
	vv_update_forces(sys);
	KE = calc_system_kinetic(sys);
	T = calc_system_temperature(sys,KE);
	nvt_update_zeta(sys,T);
	nvt_update_velocities(sys,2);
	return 0;
}



int nvt_update_velocities(MDSystem *sys,int step){
	/* Updates 1/2 step as in Velocity Verlet method, for either part */
	int i=0;
	double m,dt_over_2,correction;
	Vec3 update, vel;
	Particle* p;
	
	if (step==1){
		
		for (i=0;i<sys->N_particles;i++){
			p = &(sys->particles[i]);
			pt_get_mass(p,&m);
			dt_over_2 = (sys->dt)/2.;
			vec_copy(&(p->force),&update);
			vec_times_scalar(&update,m); // update: F(t)/m

			vec_copy(&(p->vel),&vel);
			vec_times_scalar(&vel,-(sys->zeta));
			vec_add(&vel,&update); //update: F(t)/m-zeta(t)*v(t)
			
			vec_times_scalar(&update,dt_over_2); //update: (F(t)/m-zeta(t)*v(t))*dt/2
			
			pt_add_vel(p,&update); //v(t+dt/2)=v(t)+(F(t)/m-zeta(t)*v(t))*dt/2
		}
		
	} else if (step==2){
		
		for (i=0;i<sys->N_particles;i++){
			p = &(sys->particles[i]);
			pt_get_mass(p,&m);
			dt_over_2 = (sys->dt)/2.;
			vec_copy(&(p->force),&update);
			vec_times_scalar(&update,dt_over_2/m); // update: F(t+dt)/m*dt/2
			pt_add_vel(p,&update); //v(t+dt)=v(t+dt/2)+(F(t+dt)/m*dt/2)
			
			correction = 1.+(sys->zeta)*dt_over_2;
			vec_times_scalar(&(p->vel),1./correction); // v(t+dt)=(v(t+dt/2)+(F(t+dt)/m*dt/2))/(1+zeta(t+dt)*dt/2)
		}
	}
	
	return 0;
}

int nvt_update_zeta(MDSystem *sys, double T){
	double tau = sys->tau_damp;
	sys->zeta = 1./(tau*tau)*(T/(sys->T_des)-1.)*(sys->dt);
	
	return 0;
}


