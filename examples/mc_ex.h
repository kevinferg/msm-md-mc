#define MC_EXAMPLE_STEPS 100

int mc_ex(void) {
    MDSystem sys;                          // Create system
    sys_init(&sys);                        // Initialize system with default LJ potential

    io_load_txt(&sys,"liquid256.txt");     // Load particles locations from file "liquid256.txt"
    sys_set_boxlen(&sys, 6.8);             // Periodic boundary with box side lengths 6.8
    log_init(&sys, 1, "log.txt");          // Output material properties to "log.txt" every 1 step

    sys_run_mc(&sys,                       // Run a Monte Carlo simulation...
            0.831716,                      //   Dimensionless temperature:  kT = 0.831716
            0.1,                           //   Max particle perturbation:  dr = 0.1
            MC_EXAMPLE_STEPS,              //   MC_EXAMPLE_STEPS steps
            0);                            //   0 --> Do not print progress
                                                
    io_export_xyz(&sys, "snapshot.xyz");   // Output final particle locations
    sys_destroy(&sys);                     // Free the system and close log file

    return 0;
}