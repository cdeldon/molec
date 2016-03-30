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
#include <molec/Integrator.h>

#include <math.h>

/* Global variables to keep track of different integrators */
molec_integrator integrators[32];
int num_integrators = 0;

/**
 * Register a new integrator to be tested
 *
 * @param integrator  pointer to an integrator
 */
void molec_register_integrator(molec_integrator integrator);

/**
 * Verify new integrators by comparing their output to the reference
 * implementation molec_integrator_leapfrog_refrence()
 */
TEST_CASE(molec_UnittestIntegrator)
{
    // register integrators
    molec_register_integrator(&molec_integrator_leapfrog_unroll_2);

    //  initialize simulation parameters
    const int N = 1000;
    molec_parameter_init(N);

    // generate random initial position, velocity and force vectors
    Real* x_init = molec_random_vector(N);
    Real* v_init = molec_random_vector(N);
    Real* f = molec_random_vector(N);

    // integrate using the reference implementation
    Real* x_ref = molec_copy_vector(x_init, N);
    Real* v_ref = molec_copy_vector(v_init, N);
    Real Ekin_ref;

    molec_integrator_leapfrog_refrence(x_ref, v_ref, f, &Ekin_ref, N);

    // integrage and compare with the reference implementation
    for (int i = 0; i < num_integrators; ++i)
    {
        Real* x = molec_copy_vector(x_init, N);
        Real* v = molec_copy_vector(v_init, N);
        Real Ekin;
        integrators[i](x, v, f, &Ekin, N);

        ALLCLOSE_DOUBLE(x, x_ref, N, 1e-08, 1e-05);
        ALLCLOSE_DOUBLE(v, v_ref, N, 1e-08, 1e-05);

        CLOSE_DOUBLE(Ekin, Ekin_ref, 1e-08);

        molec_free_vector(x);
        molec_free_vector(v);
    }

    molec_free_vector(x_init);
    molec_free_vector(v_init);
    molec_free_vector(f);

    molec_free_vector(x_ref);
    molec_free_vector(v_ref);
}

void molec_register_integrator(molec_integrator integrator)
{
    integrators[num_integrators] = integrator;
    num_integrators++;
}
