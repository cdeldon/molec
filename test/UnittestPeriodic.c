/*                   _
 *   _ __ ___   ___ | | ___  ___
 *  | '_ ` _ \ / _ \| |/ _ \/ __|
 *  | | | | | | (_) | |  __/ (__
 *  |_| |_| |_|\___/|_|\___|\___| - Molecular Dynamics Framework
 *
 *  Copyright (C) 2016  Carlo Del Don  (deldonc@student.ethz.ch)
 *                      Michel Breyer  (mbreyer@student.ethz.ch)
 *                      Florian Frei   (flofrei@student.ethz.ch)
 *                      Fabian Thuring (thfabian@student.ethz.ch)
 *
 *  This file is distributed under the MIT Open Source License.
 *  See LICENSE.txt for details.
 */

#include "Unittest.h"
#include <molec/Periodic.h>

/**
 * Test periodic boundary conditions of refrence implementation
 */
TEST_CASE(molec_UnittestPeriodicRefrence)
{
    // Initialize simulation
    molec_Simulation_SOA_t* sim = molec_setup_simulation_SOA();

    const Real L = molec_parameter->L;
    const int N = molec_parameter->N;
    
    // Randomly move the atoms
    for(int i = 0; i < N; ++i)
        sim->x[i] += (rand() / (Real) RAND_MAX) * L;
    
    // Apply periodic boundary conditions
    molec_periodic_refrence(sim->x, N);
    
    // Check if all the atoms are back in the bounding box
    for(int i = 0; i < N; ++i)
        CHECK(sim->x[i] <= L && sim->x[i] >= 0.0)

    molec_teardown_simulation_SOA(sim);
}
