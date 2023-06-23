#ifndef RDF_H
#define RDF_H

#include "types.h"

int rdf_compute_from_system(MDSystem* sys, int N_bins, double* r_vals, double* g_vals);
int rdf_export_from_system(MDSystem* sys, int N_bins, char* filename);

#endif