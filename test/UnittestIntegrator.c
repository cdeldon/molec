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
molec_force_integration integrators[32];
int num_integrators = 0;

/**
 * Register a new integrator to be tested
 *
 * @param integrator  pointer to an integrator
 */
void molec_register_integrator(molec_force_integration integrator);

/**
 * Compare the outpout of the registered integrators with the output of the
 * reference implementation (molec_integrator_leapfrog_refrence()).
 */
void molec_run_integrator_test();

TEST_CASE(molec_UnittestIntegrator)
{
    molec_register_integrator(&molec_integrator_leapfrog_unroll_2);
    molec_register_integrator(&molec_integrator_leapfrog_unroll_4);
    molec_register_integrator(&molec_integrator_leapfrog_unroll_8);

    //  initialize simulation parameters and run test
    molec_NAtoms = 1000;
    molec_parameter_init(molec_NAtoms);
    molec_run_integrator_test();
}

void molec_register_integrator(molec_force_integration integrator)
{
    integrators[num_integrators] = integrator;
    num_integrators++;
}

void molec_run_integrator_test()
{
    const int N = molec_parameter->N;

    // generate random initial position, velocity and force vectors
    float* x_init = molec_random_vector(N);
    float* v_init = molec_random_vector(N);
    float* f = molec_random_vector(N);

    // integrate using the reference implementation
    float* x_ref = molec_copy_vector(x_init, N);
    float* v_ref = molec_copy_vector(v_init, N);
    float Ekin_ref;

    molec_integrator_leapfrog_refrence(x_ref, v_ref, f, &Ekin_ref, N);

    // integrage and compare with the reference implementation
    for (int i = 0; i < num_integrators; ++i)
    {
        float* x = molec_copy_vector(x_init, N);
        float* v = molec_copy_vector(v_init, N);
        float Ekin;
        integrators[i](x, v, f, &Ekin, N);

        ALLCLOSE_FLOAT(x, x_ref, N, MOLEC_ATOL, MOLEC_RTOL);
        ALLCLOSE_FLOAT(v, v_ref, N, MOLEC_ATOL, MOLEC_RTOL);

        CLOSE_FLOAT(Ekin, Ekin_ref, MOLEC_ATOL);

        molec_free_vector(x);
        molec_free_vector(v);
    }

    molec_free_vector(x_init);
    molec_free_vector(v_init);
    molec_free_vector(f);

    molec_free_vector(x_ref);
    molec_free_vector(v_ref);
}
