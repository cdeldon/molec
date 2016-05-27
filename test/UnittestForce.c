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

#define TINYTEST_PRINT_ALL
#include "Unittest.h"
#include <math.h>
#include <molec/Force.h>
#include <molec/Sort.h>
#include <string.h>

// arrays storing the reference force computed with the naice N^2 algorithm
static float *f_x_reference, *f_y_reference, *f_z_reference;

#define MOLEC_MAX_FORCE_ROUTINES 10
static int molec_num_functions = 0;

// arrays of function pointers to force calculation routines
molec_force_calculation force_routines[MOLEC_MAX_FORCE_ROUTINES];
char* force_routines_name[MOLEC_MAX_FORCE_ROUTINES];

/**
 * @brief adds a force calculation routine to be checked by the test
 * @param f      function pointer of type @c molec_force_calculation
 * @param name   short description of force routine
 */
void add_function(molec_force_calculation f, char* name)
{
    if(molec_num_functions >= MOLEC_MAX_FORCE_ROUTINES)
    {
        molec_error("Couldn't register %s, too many functions registered (Max: %d)", name,
                    MOLEC_MAX_FORCE_ROUTINES);
    }

    force_routines[molec_num_functions] = f;
    force_routines_name[molec_num_functions] = name;

    ++molec_num_functions;
}

/**
 * @brief called by the test at the begin of the testcase
 *
 * All force calculation routines that have to be checked are to be added inside this function as
 * follows:
 * @code
 * add_function(&<name_of_force_routine>, "this is a short description");
 * @endcode
 */
void molec_force_test_register_functions()
{
    add_function(&molec_force_N2_refrence, "Naive N^2 implementation");
    add_function(&molec_force_cellList_knuth, "Cell list (Knut)");
    add_function(&molec_force_cellList_reference, "Cell list reference");
    add_function(&molec_force_cellList_v1, "Cell list (v1)");
    add_function(&molec_force_cellList_v2, "Cell list (v2)");
    add_function(&molec_force_quadrant, "Quadrant");
    add_function(&molec_force_quadrant_ghost, "Quadrant (ghost)");
    add_function(&molec_force_quadrant_ghost_avx, "Quadrant AVX");
    add_function(&molec_force_quadrant_ghost_fma, "Quadrant FMA");
}

/**
 * @brief computes the reference forces acting on the particles
 *
 * Computes the reference forces acting on the particles using a straight forward N2 implementation
 *
 * @param sim   simulation struct containing particle position
 * @param N     number of particles in the simulation
 */
void molec_compute_reference_forces(molec_Simulation_SOA_t* sim, const int N)
{
    float Epot;
    // Compute the reference force that acts on the atoms
    molec_force_cellList_reference(sim, &Epot, N);

    if(!f_x_reference)
    {
        // Store the computed forces as reference
        MOLEC_MALLOC(f_x_reference, N * sizeof(float));
        MOLEC_MALLOC(f_y_reference, N * sizeof(float));
        MOLEC_MALLOC(f_z_reference, N * sizeof(float));
    }

    // sort the force vectors, in order to compare with other force routines
    molec_sort_qsort_forces(sim);

    memcpy(f_x_reference, sim->f_x, N * sizeof(float));
    memcpy(f_y_reference, sim->f_y, N * sizeof(float));
    memcpy(f_z_reference, sim->f_z, N * sizeof(float));
}

/**
 * @brief Compares the forces computed with the force routine passed as argument with the reference
 * one
 *
 * @param force_routine   function pointer to force calculation routine to be tested
 * @param sim             simulation struct containing particle position and forces
 * @param N               number of particles in the simulation
 * @param description     short description of force calculation routine (used for error messages)
 */
void molec_check_forces(molec_force_calculation force_routine,
                        molec_Simulation_SOA_t* sim,
                        const int N,
                        const char* description)
{
    float Epot;

    // compute forces with routine passed as argument
    force_routine(sim, &Epot, N);

    // sort the molecules according to a common order, so that the forces are comparable
    molec_sort_qsort_forces(sim);

    // check whether the computed forces are ok
    ALLCLOSE_FLOAT_MSG(sim->f_x, f_x_reference, N, MOLEC_ATOL, MOLEC_RTOL, description)
    ALLCLOSE_FLOAT_MSG(sim->f_y, f_y_reference, N, MOLEC_ATOL, MOLEC_RTOL, description)
    ALLCLOSE_FLOAT_MSG(sim->f_z, f_z_reference, N, MOLEC_ATOL, MOLEC_RTOL, description)
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
    const int n_atoms_per_dimension = 25;
    const int r_seed = 40;

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

        // compute the forces using the routines specified above
        molec_check_forces(force_routines[r], sim, N, force_routines_name[r]);
    }

    molec_teardown_simulation_SOA(sim);
}
