int md_ex(void) {
    MDSystem sys;                            // Create system
    sys_init(&sys);                          // Initialize system with default Lennard-Jones potential

    io_load_txt(&sys, "liquid256.txt");      // Load particles locations from file "liquid256.txt"
    sys_set_boxlen(&sys, 6.8);               // Periodic boundary with box side lengths 6.8

    anim_init(&sys, 100, "trajectory.xyz");  // Output locations to "trajectory.xyz" every 100 steps
    log_init(&sys, 10, "log.txt");           // Output material properties to "log.txt" every 10 steps

    sys_random_velocities(&sys, 1.773);      // Randomize particle speeds to about 1.773 on average 
                                            // (but keep zero overall momentum)

    sys_set_dt(&sys, 0.002);                 // Set time step to 0.002 time units
    sys_nvt_ensemble(&sys, temp, 0.05);      // Set system ensemble to NVT with tau = 0.05 time units
    sys_run(&sys, 50000, 50);                // Run NVT 50,000 steps; print progress every 50 steps

    sys_zero_trupos(&sys);                   // Reset mean-squared-displacement particle locations
    sys_nve_ensemble(&sys);                  // Set system ensemble to NVE
    sys_run(&sys, 50000, 50);                // Run NVE 50,000 steps; print progress every 50 steps

    sys_destroy(&sys);                       // Free the system and close log/trajectory files

    return 0;
}