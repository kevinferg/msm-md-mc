int rdf_ex(void) {
    MDSystem sys;                          // Create system
    sys_init(&sys);                        // Initialize system with default LJ potential

    io_load_txt(&sys,"liquid256.txt");     // Load particle locations from file "liquid256.txt"
    sys_set_boxlen(&sys, 6.8);             // Periodic boundary with box side lengths 6.8

    rdf_export_from_system(&sys,           // Print g(r) for a system
                           50,             // 50 bins
                           "rdf.txt");     // Print r and g(r) to "rdf.txt"

    sys_destroy(&sys);                     // Free the system and close log file
    return 0;
}