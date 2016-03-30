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

int molec_NAtoms = 1000;

 Real* molec_init_vector(const int N)
 {
     Real* vec;
     MOLEC_MALLOC(vec, sizeof(Real) * N);
     return vec;
 }

Real* molec_random_vector(const int N)
{
    Real* vec;
    MOLEC_MALLOC(vec, sizeof(Real) * N);
    for (int i=0; i < N; i++) vec[i] = MOLEC_RANDOM;
    return vec;
}

Real* molec_copy_vector(const Real* vec, const int N)
{
    Real* vec_cpy = malloc(sizeof(Real) * N);
    memcpy(vec_cpy, vec, sizeof(Real) * N);
    return vec_cpy;
}

void molec_free_vector(Real* vec)
{
    MOLEC_FREE(vec);
}

molec_Simulation_SOA_t* molec_setup_simulation_SOA()
{
    const int N = molec_NAtoms;

    // Set parameters
    molec_parameter_init(N);

    // Allocate simulation struct and arrays
    molec_Simulation_SOA_t* sim = molec_init_simulation_SOA();

    // Set initial conditions
    molec_set_initial_condition(sim);

    return sim;
}

void molec_teardown_simulation_SOA(molec_Simulation_SOA_t* sim)
{
    molec_free_simulation_SOA(sim);
    MOLEC_FREE(molec_parameter);
}
