#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "types.h"
#include "rdf.h"
#include "calculations.h"

#define PI 3.14159265358979

double min3(double A, double B, double C) {
    double D = A < B ? A : B;
    double E = D < C ? D : C;
    return E;
}

double periodic_distance(double A, double B, double L) {
    double abs_d = fabs(A - B);
    return abs_d < L/2. ? abs_d : L - abs_d;
}

int rdf_count_pairs(MDSystem* sys) {
    return (sys->N_particles)*(sys->N_particles - 1)/2;
}

double get_dr(double r_min, double r_max, int N_bins) {
    return  (r_max - r_min) / ((double) N_bins);
}

double rdf_pair_distance(MDSystem* sys, int idxA, int idxB) {
    double dx = periodic_distance(sys->particles[idxA].pos.x, sys->particles[idxB].pos.x, sys->boxlen.x);
    double dy = periodic_distance(sys->particles[idxA].pos.y, sys->particles[idxB].pos.y, sys->boxlen.y);
    double dz = periodic_distance(sys->particles[idxA].pos.z, sys->particles[idxB].pos.z, sys->boxlen.z);
    return sqrt(dx*dx + dy*dy + dz*dz);
}
int rdf_pair_distances(MDSystem* sys, double* distances, int N_distances) {
    int i, j, index = 0;
	for (i=0; i<(sys->N_particles-1); i++) {
		for (j = i + 1; j<sys->N_particles; j++) {
            if (index >= N_distances) {
                return -1;
            }
            distances[index++] = rdf_pair_distance(sys, i, j);
        }
    }
    return 0;
}

int rdf_sort_distances(double* distances, int N_distances, int* bins, int N_bins, double r_min, double r_max) {
    int i, bin_index;
    double dr = get_dr(r_min, r_max, N_bins);

    // Empty all bins
    for (i = 0; i < N_bins; i++) {
        bins[i] = 0;
    }
    
    // Fill bins as needed
    for (i = 0; i < N_distances; i++) {
        bin_index = (int) ((distances[i]-r_min)/dr);
        if (bin_index < 0 || bin_index >= N_bins) {
            continue; // Out-of-bounds index for this r value
        } 
        bins[bin_index]++;
    }
    return 0;
}

int rdf_get_bin_centers(double* r_vals, int N_bins, double r_min, double r_max) {
    int i;
    double dr = get_dr(r_min, r_max, N_bins);
    double r = r_min + dr/2.;
    for (i = 0; i < N_bins; i++) {
        r_vals[i] = r;
        r += dr;
    }
    return 0;
}

int rdf_get_g_values(double* g_vals, int* bins, int N_bins, int N_particles, double density, double r_min, double r_max) {
    double dr = get_dr(r_min, r_max, N_bins);
    double r1, r2;
    int i;
    for (i = 0; i < N_bins; i++) {
        r1 = dr * (i);
        r2 = dr * (i + 1);
        g_vals[i] = ((double) bins[i]) / (2.*PI/3. * (r2*r2*r2 - r1*r1*r1) * density * (double) N_particles);
    }
    return 0;
}

/**********************************************************************************************************************/

int rdf_compute_from_system(MDSystem* sys, int N_bins, double* r_vals, double* g_vals) {
    if (N_bins < 1) {
        return -5;
    }
    int status = 0;
    int* bins = malloc(N_bins * sizeof(int));
    if (!bins) {
        return -10;
    }
    double N_pairs = rdf_count_pairs(sys);
    double* distances = malloc(N_pairs * sizeof(double));
    if (!distances) {
        free(bins);
        return -10;
    }
    // Get all periodic pair distances
    status += rdf_pair_distances(sys, distances, N_pairs);

    // Sort pair distances into bins
    double r_min = 0;
    double r_max = min3(sys->boxlen.x, sys->boxlen.y, sys->boxlen.z)/2;
    status += rdf_sort_distances(distances, N_pairs, bins, N_bins, r_min, r_max);

    // Get g(r) values from bin counts
    double density = calc_system_density(sys);
    status += rdf_get_g_values(g_vals, bins, N_bins, sys->N_particles, density, r_min, r_max);

    // Get r values at bin centers
    status += rdf_get_bin_centers(r_vals, N_bins, r_min, r_max);
    free(bins);
    free(distances);
	return status;
}

int rdf_export_from_system(MDSystem* sys, int N_bins, char* filename) {
    double* r_vals = malloc(N_bins * sizeof(double));
    if (!r_vals) {
        return -10;
    }
    double* g_vals = malloc(N_bins * sizeof(double));
    if (!g_vals) {
        free(r_vals);
        return -10;
    }

    int status = rdf_compute_from_system(sys, N_bins, r_vals, g_vals);
    FILE* f = fopen(filename, "w");
    if (!f) { 
        free(r_vals);
        free(g_vals);
        return -10;
    }
    int i;
    for (i = 0; i < N_bins; i++) {
        fprintf(f, "%f %f\n", r_vals[i], g_vals[i]);
    }

    fclose(f);
    free(r_vals);
    free(g_vals);
    return status;
}