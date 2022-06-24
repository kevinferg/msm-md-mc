#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "types.h"
#include "calculations.h"
#include "io.h"
#include "system.h"
#include "particle.h"
#include "vector.h"

static int count_numbers(FILE* f);
static int next_number(FILE* f, double *val);

int io_load_txt(MDSystem* sys,const char* filename) {
	FILE* f;
	f = fopen(filename, "r");
	if (f == NULL) {
		fprintf(stderr,"Could not open the input file: %s\n",filename);
		return -2;
	}
	
	int numbers = 0;
	char buf[NUM_CHARS_TO_BUFFER];
	int status;
	int i,j;
	Vec3 position;
	
	numbers = count_numbers(f);

	if ((numbers<0) || (numbers%3)) {
		fprintf(stderr,"Could not get xyz coordinates from input file: %s\n",filename);
		return -3;
	}
	
	status = sys_set_N_particles(sys,numbers/3);
	if (status<0) return status;
	
	sys_zero_all(sys);

	for (i=0; i<(numbers/3); i++){
		for (j=0; j<3; j++){
			status = next_number(f, &(position.V[j]));
			if (status<0){
				fprintf(stderr,"An error ocurred while reading the file: %s\n",filename);
				return status;
			}
		}
		pt_set_pos(&(sys->particles[i]), &position);
	}
	
	fclose(f);
	return 0;
}


static int count_numbers(FILE* f) {
	int N = 0;
	char current = 0;
	
	do {
		current = fgetc(f);
		if (IS_NUM_CHAR(current)) {
			while(IS_NUM_CHAR(current))
				current = fgetc(f);
			N++;
		}
	} while (current != EOF);
	rewind(f);
	return N;
}

static int next_number(FILE* f, double* val) {
	char current = 0;
	static char buf[NUM_CHARS_TO_BUFFER];
	int index = 0;
	int status;
	do{
		current = fgetc(f);
		if (current == EOF) {
			return -3;
		}
	} while(!IS_NUM_CHAR(current));
	
	
	while(IS_NUM_CHAR(current)) {
		buf[index] = current;
		index++;
		current = fgetc(f);
		if (index >= NUM_CHARS_TO_BUFFER)
			break;
		
	}
	buf[index] = '\0';
	status = sscanf(buf, "%lf", val);
	if (status==EOF)
		return -3;
	return 0;
}


int io_export_pdb(MDSystem* sys,const char* filename) {
	FILE *f;
	
	f = fopen(filename, "w");
	if (f==NULL) {
		fprintf(stderr,"Could not open the output file: %s\n",filename);
		return -2;
	}
	int i;
	
	for(i=0;i<sys->N_particles;i++) {
		fprintf(f,"ATOM%7d  C   PRO A%4d    %8.3f%8.3f%8.3f  1.00 10.00           C\n",
		i, 1,
		sys->particles[i].pos.x,
		sys->particles[i].pos.y,
		sys->particles[i].pos.z);
	}

	fclose(f);
	return 0;
}


int io_export_xyz(MDSystem* sys, const char* filename) {
	FILE *f;
	
	f = fopen(filename, "w");
	if (f==NULL) {
		fprintf(stderr,"Could not open the output file: %s\n",filename);
		return -2;
	}
	int i;
	fprintf(f,"%d\n",sys->N_particles);
	for(i=0;i<sys->N_particles;i++) {
		fprintf(f,"\nC %f %f %f",
		sys->particles[i].pos.x,
		sys->particles[i].pos.y,
		sys->particles[i].pos.z);
	}

	fclose(f);
	return 0;
}

int anim_export_frame(MDSystem* sys) {
	FILE* f = sys->anim_file;
	int i;

	fprintf(f,"%d\n%d\n",sys->N_particles,(int) sys->time_steps/sys->anim_every);

	for(i=0;i<sys->N_particles;i++) {
		fprintf(f,"%d %f %f %f\n",i,
		sys->particles[i].pos.x,
		sys->particles[i].pos.y,
		sys->particles[i].pos.z);
	}
	
	fflush(f);
	return 0;
}

int anim_init(MDSystem* sys, int every, const char* filename) {
	if (every < 1) {
		every = 1;
	}
	
	FILE* f;
	
	f = fopen(filename, "w");
	if (f==NULL) {
		fprintf(stderr, "Could not open the input file: %s\n", filename);
		return -2;
	}

	sys->anim_every = every;
	sys->anim_file = f;
	sys->anim_initialized = INITIALIZED;
	return 0;
}

