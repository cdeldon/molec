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

#ifndef MOLEC_UNITTEST_H
#define MOLEC_UNITTEST_H

#include <tinytest/tinytest.h>
#include <molec/InitialCondition.h>
#include <molec/Parameter.h>
#include <molec/Simulation.h>

#include <string.h>
#include <stdlib.h>

#define MOLEC_ATOL 1e-04f
#define MOLEC_RTOL 1e-02f

/**
 * Number of atoms used during unittesting (defined in Unittest.c)
 * and derider particle density
 */
extern int molec_NAtoms;
extern float molec_Rho;

/* Macro for generating random numbers */
#define MOLEC_RANDOM (((float)rand()) / RAND_MAX)

/**
 * Allocate memory for an array of size N
 */
 float* molec_init_vector(const int N);

/**
 * Generate a random array of size N
 */
float* molec_random_vector(const int N);

/**
 *  Return a deep copy of the given array
 */
float* molec_copy_vector(const float* vec, const int N);

/**
 * Free array
 */
void molec_free_vector(float*);

/**
 * @brief Setup and initialize the simulation struct using molec_NAtoms atoms
 *
 * This routine will allocate all memory needed during the simulation and set appropriate initial
 * conditions. This routine is meant to be called in every unit test (similar to SetUp and TearDown
 * of Google's Testing Framework):
 *
 * @code{.c}
 *  TEST_CASE(myTest)
 *  {
 *      // Setup simulation
 *      molec_Simulation_SOA_t* sim = molec_setup_simulation_SOA();
 *
 *      // Test ...
 *
 *      // Tear-down simulation
 *      molec_teardown_simulation_SOA(sim);
 *  }
 * @endcode
 *
 * @param return   Simulation struct holding the position, velocity and froce arrays
 */
molec_Simulation_SOA_t* molec_setup_simulation_SOA();

/**
 * Free all memory allocated during molec_setup_simulation()
 */
void molec_teardown_simulation_SOA(molec_Simulation_SOA_t* sim);

#endif
