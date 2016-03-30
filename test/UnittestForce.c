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
#include <molec/Sort.h>
#include <math.h>
#include <string.h>

// arrays storing the reference force computed with the naice N^2 algorithm
static Real* f_x_reference, *f_y_reference, *f_z_reference;


#define MOLEC_MAX_FORCE_ROUTINES 10
static int molec_num_functions = 0;

// arrays of function pointers to force calculation routines
molec_force_calculation force_routines[MOLEC_MAX_FORCE_ROUTINES];
char* force_routines_name[MOLEC_MAX_FORCE_ROUTINES];

void add_function(molec_force_calculation f, char *name)
{
    if(molec_num_functions >= MOLEC_MAX_FORCE_ROUTINES)
    {
        printf("Couldn't register %s, too many functions registered (Max: %d)",
                    name, MOLEC_MAX_FORCE_ROUTINES);
        return;
    }

    force_routines[molec_num_functions] = f;
    force_routines_name[molec_num_functions] = name;

    ++molec_num_functions;
}

void molec_force_test_register_functions()
{
    // Registers slow_filter with the driver
    add_function(&molec_force_N2_refrence, "Naive N^2 implementation");
    add_function(&molec_force_cellList, "Cell list (while loop)");
    add_function(&molec_force_cellList_for, "Cell list (for loop)");
    //add_function(&molec_force_cellList_for_swap, "Cell list with swap");

    // add here functions to be registered
}

void molec_compute_reference_forces(molec_Simulation_SOA_t* sim, const int N)
{
    Real Epot;
    // Compute the reference force that acts on the atoms
    molec_force_cellList_for(sim, &Epot, N);

    if(!f_x_reference)
    {
        // Store the computed forces as reference
        MOLEC_MALLOC(f_x_reference, N * sizeof(Real));
        MOLEC_MALLOC(f_y_reference, N * sizeof(Real));
        MOLEC_MALLOC(f_z_reference, N * sizeof(Real));
    }

    // sort the force vectors, in order to compare with other force routines
    molec_sort_qsort_forces(sim);

    memcpy(f_x_reference, sim->f_x, N * sizeof(Real));
    memcpy(f_y_reference, sim->f_y, N * sizeof(Real));
    memcpy(f_z_reference, sim->f_z, N * sizeof(Real));
}

void molec_check_forces(molec_force_calculation force_routine, molec_Simulation_SOA_t* sim, const int N)
{
    Real Epot;

    // compute forces with routine passed as argument
    force_routine(sim, &Epot, N);

    // sort the molecules according to a common order, so that the forces are comparable
    molec_sort_qsort_forces(sim);

    // check wheter the computed forces are ok
    ALLCLOSE_DOUBLE(sim->f_x, f_x_reference, N, 1e-8, 1e-8)
    ALLCLOSE_DOUBLE(sim->f_y, f_y_reference, N, 1e-8, 1e-8)
    ALLCLOSE_DOUBLE(sim->f_z, f_z_reference, N, 1e-8, 1e-8)
}

/**
 * @brief Test Force calculation implementation
 *
 * Tests the resulting force acting on each particle of the system
 * using different force calculation routines
 * The naive N^2 implementation is to be considered the reference
 */
TEST_CASE(molec_UnittestForce)
{
    const int n_atoms_per_dimension = 15;
    const int r_seed = 42;

    molec_NAtoms = pow(n_atoms_per_dimension, 3);

    // register the functions to be tested
    molec_force_test_register_functions();

    // Initialize simulation seeding the RNG
    srand(r_seed);
    molec_Simulation_SOA_t* sim = molec_setup_simulation_SOA();

    const int N = molec_parameter->N;

    // compute the reference forces
    molec_compute_reference_forces(sim, N);

    for(int r = 0; r < molec_num_functions; ++r)
    {
        // reset the configurations of the atoms
        srand(r_seed);
        molec_set_initial_condition(sim);

        //printf("Running force comparison with %s\n", force_routines_name[r]);

        // compute the forces using the routines specified above
        molec_check_forces(force_routines[r], sim, N);
    }

    molec_teardown_simulation_SOA(sim);
}


