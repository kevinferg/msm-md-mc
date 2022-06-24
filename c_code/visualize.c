#include <stdio.h>

#if defined(WIN32) || defined(_WIN32)
#include <windows.h>
#endif

#include "types.h"
#include "visualize.h"


static int cursor_position(unsigned char r, unsigned char c);
static uint snap_to_grid(double x, double y, uint* row, uint* col, double* limits);

#if defined(WIN32) || defined(_WIN32)
static DWORD initialMode;
#endif


int vis_init(void) {
	#if defined(WIN32) || defined(_WIN32)
		DWORD outMode = 0;
		HANDLE hOut = GetStdHandle(STD_OUTPUT_HANDLE);
		GetConsoleMode(hOut, &outMode);
		initialMode = outMode;

		outMode |= ENABLE_VIRTUAL_TERMINAL_PROCESSING;

		SetConsoleMode(hOut, outMode);
	#endif

	fputs("\x1b[?25l",stdout);
	fputs("\x1b[?12l",stdout);
}


int vis_end(void) {
	fputs("\x1b[?12h",stdout);
	fputs("\x1b[?25h",stdout);
	
	#if defined(WIN32) || defined(_WIN32)
		HANDLE hOut = GetStdHandle(STD_OUTPUT_HANDLE);
		SetConsoleMode(hOut,initialMode);
	#endif
}

int vis_render_system(MDSystem *sys, unsigned int dim1, unsigned int dim2) {
	if ((dim1>2) || (dim2>2) || (dim1==dim2)){
		dim1 = 0;
		dim2 = 1;
	}
	
	double limits[4]={(sys->xyz_min.V)[dim1],(sys->xyz_max.V)[dim1],
					(sys->xyz_min.V[dim2]),(sys->xyz_max.V)[dim2]};
	
	static int first=1;
	uint r,c,i;
	int now;
	static char console_buffer[NUM_CHARS];
	if (first==1){
		for (i=0;i<NUM_CHARS;i++){
			console_buffer[i] = ' ';
		}
		
		for (i=0;i<CONSOLE_HEIGHT;i++) {
			console_buffer[(i+1)*(CONSOLE_WIDTH+1)] = '\n';
		}
		console_buffer[NUM_CHARS-1] = '\0';
	}


	for (i=0;i<sys->N_particles;i++) {
		snap_to_grid(sys->particles[i].pos.V[dim1],sys->particles[i].pos.V[dim2],&r,&c,limits);
		console_buffer[r*(CONSOLE_WIDTH+1)+c] = 'O';
	}
	
	cursor_position(0,0);
	fputs(console_buffer,stdout);
	fflush(stdout);
	
	for (i=0;i<sys->N_particles;i++) {
		snap_to_grid(sys->particles[i].pos.V[dim1],sys->particles[i].pos.V[dim2],&r,&c,limits);
		console_buffer[r*(CONSOLE_WIDTH+1)+c] = ' ';
	}
	
	first = 0;
	return 0;
} 

static uint snap_to_grid(double x, double y, uint* row, uint* col, double* limits) {
	double x_min = limits[0];
	double x_max = limits[1];
	double y_min = limits[2];
	double y_max = limits[3];
	
	if (x<x_min) x = x_min+1e-6;

	if (y<y_min) y = y_min+1e-6;
	
	if (x>x_max) x = x_max-1e-6;
	
	if (y>y_max) y = y_max-1e-6;

	*row = (uint) (1.+(double) CONSOLE_HEIGHT*((y_max-y)/(y_max-y_min)));
	*col = (uint) (1.+(double) CONSOLE_WIDTH*((x-x_min)/(x_max-x_min)));
	return 0;
}

static int cursor_position(unsigned char r, unsigned char c) {
	static char str[STR_SIZE];
	sprintf(str,"\x1b[%u;%uH",r,c);
	fputs(str,stdout);
	fflush(stdout);
	return 0;
}
