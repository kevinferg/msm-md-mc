#ifndef _LOGGING_H_
#define _LOGGING_H_

#include "types.h"

int log_init(MDSystem* sys, int every, const char* filename);
int log_print_line(MDSystem* sys);


#endif