#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "md_all.h"

double check_potential(MDSystem *sys);

int run_at_temp(double temp, int num);

int main(int argc, char** argv){
	MDSystem sys;
	sys_init(&sys);
	io_load_txt(&sys,"liquid256.txt");
	sys_set_boxlen(&sys, 6.8);
	log_init(&sys,1,"mc_test_log3.txt");
	srand(1);
	sys_run_mc(&sys, 0.831716, 0.1, 2.5e6,0);
	io_export_xyz(&sys, "mc_final_frame2.xyz");
	sys_destroy(&sys);
	printf("MC: Success.\n");

	
	run_at_temp(0.831716,100);
	printf("MD: Success.\n");
	return 0;
}

int run_at_temp(double temp, int num){
	MDSystem sys;
	sys_init(&sys);
	
	io_load_txt(&sys,"liquid256.txt");

	sys_set_boxlen(&sys, 6.8);
	
	char log_filename[40];
	char anim_filename[40];
	sprintf(log_filename,"nvt_log_%d.txt",num);
	sprintf(anim_filename,"nvt_anim_%d.xyz",num);

	srand(0);
	sys_random_velocities(&sys,1.773);
	

	
	anim_init(&sys,100,anim_filename);
	log_init(&sys,10,log_filename);
	
	
	sys_set_dt(&sys,0.002);
	sys_nvt_ensemble(&sys,temp,0.05);
	sys_run(&sys,50000,50);
	
	sys_zero_trupos(&sys);
	sys_nve_ensemble(&sys);
	sys_run(&sys,50000,50);
	
	
	sys_destroy(&sys);
	printf("Success!\n");
	return 0;

}



double check_potential(MDSystem *sys){
	FILE* f;
	int i;
	double r;
	f=fopen("check_potential.txt","w");
	for (i=0;i<1000;i++){
		r = (double)((5./1000.) * (double) (i+1));
		fprintf(f,"%f %f %f\n",r,U_get(sys->potential,r),F_get(sys->potential,r));
	}
	fclose(f);
	printf("Potential (and force) values for 0 saved to 'check_potential.txt')\n");
	printf("...format of each row: (r) (U) (F)\n");
	return 0;
}

