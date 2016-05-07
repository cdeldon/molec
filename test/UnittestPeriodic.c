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

/* Global variables to keep track of different periodics */
molec_periodic periodics[32];
int num_periodic = 0;

/**
 * Register a new periodic to be tested
 *
 * @param periodic  pointer to an periodic
 */
void molec_register_periodic(molec_periodic periodic);

/**
 * Compare the outpout of the registered perioidc with the output of the
 * reference implementation (molec_periodic_refrence()).
 */
void molec_run_periodic_test();

TEST_CASE(molec_UnittestPeriodic)
{
    molec_register_periodic(&molec_periodic_refrence);
    molec_register_periodic(&molec_periodic_close);
    molec_register_periodic(&molec_periodic_close4);

    molec_run_periodic_test();
}

void molec_register_periodic(molec_periodic periodic)
{
    periodics[num_periodic] = periodic;
    num_periodic++;
}

void molec_run_periodic_test()
{
    // Initialize simulation
    molec_Simulation_SOA_t* sim = molec_setup_simulation_SOA();

    const float L_x = molec_parameter->L_x;

    const int N = molec_parameter->N;

    for (int i = 0; i < num_periodic; ++i)
    {
        // Randomly move the atoms
        for(int i = 0; i < N; ++i)
            sim->x[i] += (rand() / (float) RAND_MAX) * L_x;

        // Apply periodic boundary conditions
        periodics[i](sim->x, N, L_x);

        // Check if all the atoms are back in the bounding box
        for(int i = 0; i < N; ++i)
        {
            CHECK_LE_FLOAT(sim->x[i], L_x);
            CHECK_GE_FLOAT(sim->x[i], 0.0);
        }
    }

    molec_teardown_simulation_SOA(sim);
}
