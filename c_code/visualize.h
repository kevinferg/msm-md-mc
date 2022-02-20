#ifndef _VISUALIZATION_H_
#define _VISUALIZATION_H_

#define CONSOLE_WIDTH 98
#define CONSOLE_HEIGHT 49
#define NUM_CHARS ((CONSOLE_WIDTH+1)*CONSOLE_HEIGHT+1)
#define ENABLE_VIRTUAL_TERMINAL_PROCESSING 0x0004
#define STR_SIZE 20

#include "types.h"

int vis_init(void);
int vis_end(void);

int vis_render_system(MDSystem *sys,unsigned int dim1,unsigned int dim2);




#endif