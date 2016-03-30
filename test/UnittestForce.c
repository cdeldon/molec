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
#include <molec/Force.h>
#include <math.h>

/**
 * @brief Test Force calculation implementation
 *
 * Tests the resulting force acting on each particle of the system
 * using different force calculation routines
 * The naive N^2 implementation is to be considered the reference
 */
TEST_CASE(molec_UnittestForce)
{
    const int n_atoms_per_dimension = 8;
    const int r_seed = 42;
    Real Epot;

    molec_NAtoms = pow(n_atoms_per_dimension, 3);

    // Initialize simulation seeding the RNG
    srand(r_seed);
    molec_Simulation_SOA_t* sim = molec_setup_simulation_SOA();


    const int N = molec_parameter->N;
    if(N != molec_NAtoms)
        molec_error("Number of atoms allocated by simulation setup differs from the desired one\n");

    // Compute the reference force that acts on the atoms
    molec_force_N2_refrence(sim, &Epot, N);

    // Store the computed forces as reference



    molec_teardown_simulation_SOA(sim);
}


