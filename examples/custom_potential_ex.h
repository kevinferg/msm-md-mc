
typedef struct MorseParameters {
   double D_e; // Well depth
   double   a; // Well width (larger a = wider well)
   double r_e; // Equilibrium distance
} MorseParameters;

double U_morse(double r, const void* params) {
	const MorseParameters* P = (const MorseParameters*) params;
    double part = 1 - exp(P->a * (P->r_e - r));
    return P->D_e * part * part;
}
double F_morse(double r, const void* params) {
    const MorseParameters* P = (const MorseParameters*) params;
    double part = exp(P->a * (P->r_e - r));
    return -2 * P->a * P->D_e * part * (1 - part);
}

int custom_potential_ex(void) {
    MDSystem sys;                       // Create system
    sys_init(&sys);                     // Initialize system 
                                        // (with default LJ potential)
    MorseParameters params = {          // Define Morse parameters
        .D_e = 1, .a = 4, .r_e = 1.1
    };
    Potential morse;                    // Create a new potential
    potential_init(&morse,              // Assign the potential:
                    2.5,                // - Cutoff radius: 2.5
                    &U_morse, &F_morse, // - Energy and Force functions
                    &params);           // - Morse params: D_e, a, r_e

    sys.potential = &morse;             // Apply our new Morse potential
                                        // to each pair in the system

    /* Export [r, U, F] table to potential.log */
    check_potential(&sys, "potential.log");
    return 0;
}