#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "md_all.h"

#include "../examples/custom_potential_ex.h"

int main(int argc, char** argv) {

	custom_potential_ex();

	/*
	// Export pair potential and force for a range of radii
	MDSystem sys;
	sys_init(&sys);
	check_potential(&sys, "potential.log");
	sys_destroy(&sys);
	printf("Potential check: Success.\n\n");

	// Run a Monte Carlo simulation
	mc_simulation("liquid256.txt", "mc_results.log", "mc_snapshot.xyz",
	              5, 6.8, 0.83, 0.1,
				  100, 0, 0);
	printf("MC: Success.\n\n");

	// Run a molecular dynamics simulation
	md_simulation("liquid256.txt", "md_results.log", "md_trajectory.xyz",
	              5, 10, 0.002,
				  6.8, 0.83, 1.77, 0.05,
				  100, 100, 10, 0);
	printf("MD: Success.\n\n");
	*/
	return 0;
}
