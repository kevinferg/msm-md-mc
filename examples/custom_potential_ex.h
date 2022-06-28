
double F_morse(double r, double* params) {
    double De = params[0];   // Well depth
    double a = params[1];    // Well width
    double r_eq = params[2]; // Equilibrium bond distance
    double part = exp(a * (r_eq - r));
    return -2 * a * De * part * (1 - part);
}

double U_morse(double r, double* params) {
    double De = params[0];   // Well depth
    double a = params[1];    // Well width
    double r_eq = params[2]; // Equilibrium bond distance
    double part = 1 - exp(a * (r_eq - r));
    return De * part * part;
}

int custom_potential_ex(void) {
    MDSystem sys;
    sys_init(&sys);
    Potential morse;
    potential_init(&morse, 2.5, 
                   &U_morse, &F_morse,
                   (double[3]) {1, 4, 1.1});
    sys.potential = &morse;
	check_potential(&sys, "potential.log");
}