#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "types.h"
#include "logging.h"
#include "calculations.h"
#include "system_stuff.h"
#include "vec_stuff.h"

int log_init(MDSystem* sys, int every, const char* filename){
	if (every < 1){
		every = 1;
	}
	
	static FILE* f;
	
	f = fopen(filename, "w");
	if (f==NULL){
		fprintf(stderr,"Could not open the input file: %s\n",filename);
		return -2;
	}
	
	fprintf(f,"time_step  time  K  U  H  kT  p  msd \n");//  mv_x  mv_y  mv_z\n");
	
	sys->log_every = every;
	sys->log_file = f;
	sys->log_initialized = INITIALIZED;
	printf("System will log to [%s] every %d time steps.\n",filename,every);
	return 0;
}

int log_print_line(MDSystem* sys){
	FILE* f = sys->log_file;
	double t,K,U,kT,pressure, msd;
	long ts;
	
	ts = sys->time_steps;
	//Vec3 momentum;
	t = (sys->dt) * (double) ts;
	K = calc_system_kinetic(sys);
	U = calc_system_potential(sys);
	msd = calc_system_msd(sys);
	
	kT = calc_system_temperature(sys,K);
	pressure = calc_system_pressure(sys,kT,sys->virial);

	//calc_system_momentum(sys,&momentum);

	fprintf(sys->log_file,"%d %f %f %f %f %f %f %f\n",//  %f %f %f\n",
			ts,t,
			K,U,K+U,kT,pressure,msd); //,
			//momentum.x,momentum.y,momentum.z);
	return 0;
}