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

molec_Simulation_SOA_t* molec_setup_simulation_SOA()
{
    const int N = molec_NAtoms;

    // Set parameters
    molec_parameter_init(N);

    // Allocate simulation struct and arrays
    molec_Simulation_SOA_t* sim = malloc(sizeof(molec_Simulation_SOA_t));

    MOLEC_MALLOC(sim->x, sizeof(Real) * N);
    MOLEC_MALLOC(sim->y, sizeof(Real) * N);
    MOLEC_MALLOC(sim->z, sizeof(Real) * N);

    MOLEC_MALLOC(sim->v_x, sizeof(Real) * N);
    MOLEC_MALLOC(sim->v_y, sizeof(Real) * N);
    MOLEC_MALLOC(sim->v_z, sizeof(Real) * N);

    MOLEC_MALLOC(sim->f_x, sizeof(Real) * N);
    MOLEC_MALLOC(sim->f_y, sizeof(Real) * N);
    MOLEC_MALLOC(sim->f_z, sizeof(Real) * N);

    // Set initial conditions
    molec_set_initial_condition(sim);

    return sim;
}

void molec_teardown_simulation_SOA(molec_Simulation_SOA_t* sim)
{
    MOLEC_FREE(sim->x);
    MOLEC_FREE(sim->y);
    MOLEC_FREE(sim->z);

    MOLEC_FREE(sim->v_x);
    MOLEC_FREE(sim->v_y);
    MOLEC_FREE(sim->v_z);

    MOLEC_FREE(sim->f_x);
    MOLEC_FREE(sim->f_y);
    MOLEC_FREE(sim->f_z);

    MOLEC_FREE(sim);
    MOLEC_FREE(molec_parameter);
}

