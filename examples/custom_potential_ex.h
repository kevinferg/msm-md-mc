
double U_morse(double r, double* params) {
    double D_e = params[0];  // Well depth
    double a = params[1];    // Well width
    double r_e = params[2];  // Equilibrium bond distance
    double part = 1 - exp(a * (r_e - r));
    return D_e * part * part;
}

double F_morse(double r, double* params) {
    double D_e = params[0];  // Well depth
    double a   = params[1];  // Well width
    double r_e = params[2];  // Equilibrium bond distance
    double part = exp(a * (r_e - r));
    return -2 * a * D_e * part * (1 - part);
}

int custom_potential_ex(void) {
    MDSystem sys;                      // Create system
    sys_init(&sys);                    // Initialize system 
                                       // (with default LJ potential)
    
    Potential morse;                   // Create a new potential
    potential_init(&morse,             // Assign the potential:
                   2.5,                // - Cutoff radius: 2.5
                   &U_morse, &F_morse, // - Energy and Force functions
      (double[3])  {1, 4, 1.1});       // - Morse params: D_e, a, r_e

    sys.potential = &morse;            // Apply our new Morse potential
                                       // to each pair in the system

    /* Export [r, U, F] table to potential.log */
	check_potential(&sys, "potential.log");
    return 0;
}