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

/**
 * @brief adds a force calculation routine to be checked by the test
 * @param f      function pointer of type @c molec_force_calculation
 * @param name   short description of force routine
 */
void add_function(molec_force_calculation f, char *name)
{
    if(molec_num_functions >= MOLEC_MAX_FORCE_ROUTINES)
    {
        molec_error("Couldn't register %s, too many functions registered (Max: %d)",
                    name, MOLEC_MAX_FORCE_ROUTINES);
    }

    force_routines[molec_num_functions] = f;
    force_routines_name[molec_num_functions] = name;

    ++molec_num_functions;
}

/**
 * @brief called by the test at the begin of the testcase
 *
 * All force calculation routines that have to be checked are
 * to be added inside this function as follows:
 * @code
 * add_function(&<name_of_force_routine>, "this is a short description");
 * @endcode
 */
void molec_force_test_register_functions()
{
    // Registers slow_filter with the driver
    add_function(&molec_force_N2_refrence, "Naive N^2 implementation");
    add_function(&molec_force_cellList, "Cell list (while loop)");
    add_function(&molec_force_cellList_for, "Cell list (for loop)");
    add_function(&molec_force_reference_dp,"Cell List(for loop Flo");
    add_function(&molec_force_cellList_for_swap, "Cell list with swap");
    add_function(&molec_force_cellList_dummy, "dummy cell list implementation");

    // add here functions to be registered
}

/**
 * @brief computes the reference forces acting on the particles
 *
 * Computes the reference forces acting on the particles using a
 * straight forward N2 implementation
 *
 * @param sim   simulation struct containing particle position
 * @param N     number of particles in the simulation
 */
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

/**
 * @brief Compares the forces computed with the force routine passed as argument with the reference one
 *
 * @param force_routine   function pointer to force calculation routine to be tested
 * @param sim             simulation struct containing particle position and forces
 * @param N               number of particles in the simulation
 * @param description     short description of force calculation routine (used for error messages)
 */
void molec_check_forces(molec_force_calculation force_routine, molec_Simulation_SOA_t* sim,
                        const int N, const char* description)
{
    Real Epot;

    // compute forces with routine passed as argument
    force_routine(sim, &Epot, N);

    // sort the molecules according to a common order, so that the forces are comparable
    molec_sort_qsort_forces(sim);

    // check wheter the computed forces are ok
    ALLCLOSE_DOUBLE_MSG(sim->f_x, f_x_reference, N, 1e-08, 1e-05, description)
    ALLCLOSE_DOUBLE_MSG(sim->f_y, f_y_reference, N, 1e-08, 1e-05, description)
    ALLCLOSE_DOUBLE_MSG(sim->f_z, f_z_reference, N, 1e-08, 1e-05, description)
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
    const int n_atoms_per_dimension = 18;
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


