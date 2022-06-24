#include <stdio.h>
#include <stdlib.h>

#include "utils.h"
#include "md_all.h"


int check_potential(MDSystem *sys, const char* filename) {
	FILE* f;
	int i;
	double r;
	f=fopen(filename,"w");
    if (f == NULL) {
        fprintf(stderr, "Could not open/create %s\n", filename);
        return -1;
    }
	for (i=0; i<1000; i++) {
		r = (double)((5./1000.) * (double) (i+1));
		fprintf(f,"%f %f %f\n", r, U_get(sys->potential,r), F_get(sys->potential,r));
	}
	fclose(f);
	printf("Potential (and force) values saved to %s)\n", filename);
	printf("...format of each row: (r) (U) (F)\n");
	return 0;
}

int md_simulation(const char* input_file, const char* log_file, const char* traj_file, 
                  int log_every, int traj_every, double time_step,
                  double boxlen, double temp, double init_speed, double tau,
                  long nvt_steps, long nve_steps, int progress, unsigned int seed) {
	MDSystem sys;
	sys_init(&sys);
	
	io_load_txt(&sys, input_file);
	sys_set_boxlen(&sys, boxlen);

	srand(seed);
	sys_random_velocities(&sys, init_speed);

    log_init(&sys, log_every, log_file);
	anim_init(&sys, traj_every, traj_file);

	sys_set_dt(&sys, time_step);
	sys_nvt_ensemble(&sys, temp, tau);
	sys_run(&sys, nvt_steps, progress);
	
	sys_zero_trupos(&sys);
	sys_nve_ensemble(&sys);
	sys_run(&sys, nve_steps, progress);
	
	sys_destroy(&sys);
	return 0;
}

int mc_simulation(const char* input_file, const char* log_file, const char* snapshot_file, 
                  int log_every, double boxlen, double temp, double max_trial_move,
                  long num_steps, int verbose, unsigned int seed) {
	MDSystem sys;
	sys_init(&sys);
	io_load_txt(&sys, input_file);
	sys_set_boxlen(&sys, boxlen);
	log_init(&sys, log_every, log_file);

	srand(seed);
	sys_run_mc(&sys, temp, max_trial_move, num_steps, verbose);
	io_export_xyz(&sys, snapshot_file);
    
	sys_destroy(&sys);
    return 0;
}