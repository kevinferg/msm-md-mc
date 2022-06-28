
double U_morse(double r, double* params) {
    double De = params[0];   // Well depth
    double a = params[1];    // Well width
    double r_eq = params[2]; // Equilibrium bond distance
    double part = 1 - exp(a * (r_eq - r));
    return De * part * part;
}

double F_morse(double r, double* params) {
    double De = params[0];   // Well depth
    double a = params[1];    // Well width
    double r_eq = params[2]; // Equilibrium bond distance
    double part = exp(a * (r_eq - r));
    return -2 * a * De * part * (1 - part);
}

int custom_potential_ex(void) {
    MDSystem sys;                      // Create system
    sys_init(&sys);                    // Initialize system 
                                       // (with default LJ potential)
    
    Potential morse;                   // Create a new potential
    potential_init(&morse,             // Assign the potential:
                   2.5,                // - Cutoff radius: 2.5
                   &U_morse, &F_morse, // - Energy and Force functions
      (double[3])  {1, 4, 1.1});       // - Morse params: De, a, r_eq

    sys.potential = &morse;            // Apply our new Morse potential
                                       // to each pair in the system

    /* Export [r, U, F] table to potential.log */
	check_potential(&sys, "potential.log");
}