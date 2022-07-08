#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "types.h"
#include "vector.h"
#include "verlet.h"
#include "particle.h"
#include "potentials.h"
#include "system.h"
#include "calculations.h"
#include "logging.h"
#include "io.h"

double rand_double(double lower, double upper) {
	double number;
	
	number = (double) ((double) rand()/ (double) RAND_MAX);
	number = lower + number*(upper - lower);
	
	return number;
}

int sys_zero_forces(MDSystem *sys) {
	int i;
	for (i = 0; i < sys->N_particles; i++){
		pt_set_force(&(sys->particles[i]), &ZERO_VEC);
	}
	return 0;
}

int sys_unit_masses(MDSystem *sys) {
	int i;
	for (i = 0; i < sys->N_particles; i++){
		pt_set_mass(&(sys->particles[i]), 1.0);
	}
	return 0;
}

int sys_zero_positions(MDSystem* sys) {
	int i;
	for (i = 0; i < sys->N_particles; i++){
		pt_set_pos(&(sys->particles[i]), &ZERO_VEC);
	}
	return 0;
}

int sys_zero_trupos(MDSystem* sys) {
	int i;
	for (i = 0; i < sys->N_particles; i++){
		vec_copy(&ZERO_VEC, &((sys->particles[i]).trupos));
	}
	return 0;
}



int sys_zero_velocities(MDSystem* sys) {
	int i;
	for (i = 0; i < sys->N_particles; i++){
		pt_set_vel(&(sys->particles[i]), &ZERO_VEC);
	}
	return 0;
}

int sys_random_velocities(MDSystem* sys, double A) {
	// Randomizes velocities between (approximately) -A and A
	// Also ensures total momentum is zero
	Vec3 vel;
	int i;
	if (A<0) A *= -1.;
	
	for (i=0; i<sys->N_particles; i++){
		vel.x = rand_double(-A,A);
		vel.y = rand_double(-A,A);
		vel.z = rand_double(-A,A);
		pt_set_vel(&(sys->particles[i]),&vel);
	}
	sys_zero_momentum(sys);
	return 0;
}

int sys_zero_momentum(MDSystem* sys) {
	Vec3 avg_mv, vel;
	int i;
	calc_system_momentum(sys, &avg_mv);
	for (i = 0; i < sys->N_particles; i++){
		vec_copy(&avg_mv, &vel);
		vec_times_scalar(&vel, -1./(((sys->particles)[i]).mass));
		vec_add(&vel, &(((sys->particles)[i]).vel));
	}
	return 0;
}

int sys_zero_all(MDSystem* sys) {
	sys_zero_positions(sys);
	sys_zero_trupos(sys);
	sys_zero_velocities(sys);
	sys_zero_forces(sys);
	sys_unit_masses(sys);
	sys->time_steps = 0;
	sys->runtime = 0.;
	return 0;
}

int sys_set_N_particles(MDSystem* sys,int N) {
	if (sys->particles_initialized != INITIALIZED) {
		sys->particles = malloc(sizeof(Particle)*N);
		if (sys->particles == NULL) {
			fprintf(stderr, "Failed to allocate memory for particle list!\n");
			return -1;
		}
		
	} else {
		sys->particles = realloc(sys->particles, sizeof(Particle)*N);
		if (sys->particles == NULL) {
			fprintf(stderr, "Failed to allocate memory for particle list!\n");
			return -1;
		}
	}
	sys->N_particles = N;
	sys->particles_initialized = INITIALIZED;
	return 0;
}

int sys_init(MDSystem* sys) {
	sys_set_N_particles(sys, 10);
	sys_zero_all(sys);
	sys_set_dt(sys, 0.002);
	sys_set_boxlen(sys,1.);
	
	static LJParams lj_params = {.epsilon = 1., .sigma = 1.};
	static Potential lj;
	potential_init(&lj, 2.5, &U_lj, &F_lj, &lj_params);
	sys->potential = &lj;

	sys->log_initialized = 0;
	sys->time_steps = 0;
	sys->anim_initialized = 0;
	
	sys->nvt_on = 0;
	sys->zeta = 0;

	return 0;
}

int sys_nvt_ensemble(MDSystem *sys, double T_des, double tau) {
	sys->nvt_on = 1;
	sys->zeta = 0;
	sys->tau_damp = tau;
	sys->T_des = T_des;
	return 0;
}

int sys_nve_ensemble(MDSystem *sys) {
	sys->nvt_on = 0;
	return 0;
}


int sys_set_dt(MDSystem* sys, double dt) {
	sys->dt = dt;
	return 0;
}

int sys_run(MDSystem* sys, long numsteps, int progress) {
	long i;
	double time;
	long time_left;
	clock_t t;
	t = clock();


	if (sys->time_steps == 0)
		vv_update_forces(sys);
	
	for (i = 0; i < numsteps; i++) {
		
		if (sys->log_initialized == INITIALIZED) { 
			if (!((sys->time_steps) % (sys->log_every))) {
				log_print_line(sys);
			}
		}
		if (sys->anim_initialized == INITIALIZED) { 
			if (!((sys->time_steps) % (sys->anim_every))) {
				anim_export_frame(sys);
			}
		}

		vv_step(sys);
		
		if (progress > 0) {
			if (!((sys->time_steps) % progress)) {
				time = ((double) clock()-t)/CLOCKS_PER_SEC;
				time_left = (long) ((double)(numsteps-i)*time / (double) (i+1));
				printf("\rProgress: %4.1f%%,    Est. time remaining: %3dmin %2dsec",
				       100.*(double) (i+1)/(double) numsteps, 
					   (int) time_left / 60, (int) time_left % 60);
				fflush(stdout);
			}
		}
		
		(sys->time_steps)++;
	}
	
	time = ((double) clock()-t)/CLOCKS_PER_SEC;
	(sys->runtime) += time;
	if (progress > 0) {
		printf("\rDone.    Time elapsed: %3dmin %2dsec                               \n",
	    (int) ((long) time)/60, (int) ((long) time)%60);
		fflush(stdout);
	}
	return 0;
}

int sys_destroy(MDSystem *sys) {
	if (sys->particles_initialized == INITIALIZED) {
		free(sys->particles);
		sys->particles_initialized = 0;
	}

	if (sys->log_initialized == INITIALIZED) {
		fclose(sys->log_file);
		sys->log_initialized = 0;
	}
	
	if (sys->anim_initialized == INITIALIZED) {
		fclose(sys->anim_file);
		sys->anim_initialized = 0;
	}
	
	sys->N_particles = 0;

	return 0;
}

int sys_print_positions(MDSystem *sys) {
	int i;
	for (i = 0; i < sys->N_particles; i++){
		vec_print(&(sys->particles[i].pos));
	}
	fflush(stdout);
	return 0;
}

int sys_set_boxlen(MDSystem *sys,double boxlen) {
	int status = 0;
	if (boxlen <= 0.){
		boxlen = 1.;
		status = -1;
	}

	(sys->boxlen).x = boxlen; 
	(sys->boxlen).y = boxlen;
	(sys->boxlen).z = boxlen; // Now, sys->boxlen is (L,L,L)

	vec_copy(&ZERO_VEC, &(sys->xyz_min)); // min: (0,0,0)
	vec_copy(&ZERO_VEC, &(sys->xyz_max)); // max: (0,0,0)
	vec_add(&(sys->boxlen), &(sys->xyz_max)); // max: (L,L,L)

	return status;
}

double sys_get_runtime(MDSystem* sys) {
	return sys->runtime;
}

int sys_get_N_particles(MDSystem* sys) {
	return sys->N_particles;
}

int sys_set_density(MDSystem* sys, double density) {
	double current = calc_system_density(sys);
	double scale_factor = current / density;
	scale_factor = pow(scale_factor, 1./3.);

	Vec3 new_L;
	vec_copy(&(sys->boxlen),&new_L);
	vec_times_scalar(&new_L, scale_factor);

	sys_rescale(sys, &new_L);
	return 0;
}

int sys_rescale(MDSystem* sys, Vec3 *new_L) {
	int i, j;
	double pos, minpos, frac;
	for (i = 0; i < sys->N_particles; i++) {
		for (j = 0; j < 3; j++){
			pos = (sys->particles[i]).pos.V[j];
			minpos = sys->xyz_min.V[j];
			frac = (pos-minpos) / (sys->boxlen.V[j]);
			(sys->particles[i]).pos.V[j] = minpos + frac * (new_L->V)[j];
		}
	}
	
	vec_copy(new_L, &(sys->boxlen));
	vec_copy(&(sys->xyz_min), &(sys->xyz_max));
	vec_add(new_L, &(sys->xyz_max));
	return 0;
}

int sys_run_mc(MDSystem* sys, double kT, double maxstep, unsigned long numsteps, int verbose) {
	unsigned long i;
	int j, axis;
	double virial, old_virial_part, new_virial_part;
	double pressure, pressure_avg = 0;
	double U, U_current, U_trial, dU, U_avg = 0;
	double r_current, r_trial;
	double prob, randnum;
	unsigned long accepts = 0;
	Vec3 sep_current, sep_trial;
	double time;
	clock_t t;
	t = clock();
	
	printf("\nRunning MC with a max step of %f in x, y, and z\n at kT=%f for %lu steps\n",
	       maxstep, kT, numsteps);
	fflush(stdout);

	int particle_id;
	Particle p, *p1, *p2;
	
	virial = calc_system_virial(sys);
	pressure = calc_system_pressure(sys, kT, virial);
	U = calc_system_potential(sys);
	
	for (i = 0; i < numsteps; i++) {
		// Select random particle to perturb
		particle_id = (unsigned int) rand_double((double) 0, (double) sys->N_particles);
		particle_id = (particle_id < sys->N_particles)? particle_id: sys->N_particles - 1;
		p1 = &((sys->particles)[particle_id]);
		pt_set_pos(&p,&(p1->pos));
		
		// Calculate trial move as perturbed particle location
		for (axis=0; axis<3; axis++){
			(p.pos.V[axis]) += rand_double(-maxstep,maxstep);
		}
		
		if (verbose) {
			printf("\n\n----- Iteration: %lu -----\n", i);
		    printf("Particle to perturb: %d\n", particle_id);
		    printf(" Original location: [%6.3f,%6.3f,%6.3f]\n",
			       p1->pos.x, p1->pos.y, p1->pos.z);
		}
		calc_correct_position(sys, &p); // Wrap across PBCs
		if (verbose) printf("  Trial   location: [%6.3f,%6.3f,%6.3f]\n", p.pos.x, p.pos.y, p.pos.z);
		
		// Calc energies for that particle in these 2 states
		U_current = calc_single_potential(sys,particle_id,p1);
		U_trial = calc_single_potential(sys,particle_id,&p);
		
		// Calc probability of acceptance
		dU = U_trial-U_current;
		prob = exp(-1./kT*dU);
		prob = (prob>1.)?1.:prob;
		if (verbose) printf("dU = %f --> Probability: %f\n", dU, prob);
		
		// Accept or reject and update particle position and total energy
		randnum = rand_double(0., 1.);
		if (randnum<prob) {
			if (verbose) printf("Random number %f --> ACCEPT trial move\n", randnum);
			accepts++;
			old_virial_part = calc_single_virial(sys,particle_id);
			pt_set_pos(p1,&(p.pos));
			new_virial_part = calc_single_virial(sys,particle_id);
			U = U + dU;
			virial = virial - old_virial_part + new_virial_part;
			pressure = calc_system_pressure(sys,kT,virial);
			
		} else {
			if (verbose) printf("Random number %f --> REJECT trial move\n", randnum);
		}
		fprintf(sys->log_file,"%lu %f %f\n", i, U, pressure);
		
		if (verbose) {
			printf("Energy:   U = %f\n",U);
			printf("Pressure: P = %f\n",pressure);
			fflush(stdout);
		}
		
		U_avg += U/((double) numsteps);
		pressure_avg += pressure/((double) numsteps);
	}
	
	
	time = ((double) clock()-t)/CLOCKS_PER_SEC;
	(sys->runtime) += time;
	
	printf("\n----------RESULTS-----------\n");
	printf("%lu accepts of %lu total steps\n", accepts, numsteps);
	printf("Time elapsed: %0.1f seconds\n", time);
	printf("Average energy: %f\n", U_avg);
	printf("Average pressure: %f\n", pressure_avg);
	return 0;
}