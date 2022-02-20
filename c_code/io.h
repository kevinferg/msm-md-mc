#ifndef _IO_H_
#define _IO_H_

#define NUM_CHARS_TO_BUFFER 32
#define IS_NUM_CHAR(c) ((current=='-') || (current == '.') || (current>='0' && current <='9'))

#include "types.h"

int io_load_txt(MDSystem* sys,const char* filename);


int io_export_pdb(MDSystem* sys,const char* filename);
int io_export_xyz(MDSystem* sys,const char* filename);
int anim_export_frame(MDSystem* sys);
int anim_init(MDSystem* sys, int every, const char* filename);

#endif