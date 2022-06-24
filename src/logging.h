#ifndef LOGGING_H
#define LOGGING_H

#include "types.h"

int log_init(MDSystem* sys, int every, const char* filename);
int log_print_line(MDSystem* sys);


#endif