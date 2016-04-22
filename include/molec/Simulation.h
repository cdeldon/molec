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

#ifndef MOLEC_SIMULATION_H
#define MOLEC_SIMULATION_H

#include <molec/Common.h>

/**
 * @brief Struct of arrays of the position, velocity and force arrays
 *
 * The size of the arrays can be queried using molec_parameter->N.
 */
typedef struct molec_Simulation_SOA
{
    /** Position */
    float* x;
    float* y;
    float* z;

    /** Velocity */
    float* v_x;
    float* v_y;
    float* v_z;

    /** Force */
    float* f_x;
    float* f_y;
    float* f_z;

} molec_Simulation_SOA_t;



/**
 * Return a pointer to an allocated simulation SOA
 */
molec_Simulation_SOA_t* molec_init_simulation_SOA();

/**
 * Free memory of a simulation SOA
 */
void molec_free_simulation_SOA(molec_Simulation_SOA_t* simulation);

/**
 * Print the current state of the simulation
 */
void molec_print_simulation_SOA(const molec_Simulation_SOA_t* simulation);

/**
 * Run the refrence version of the MD-Simulation using N^2 force computation
 *
 */
void molec_run_simulation(void (*molec_compute_force)( molec_Simulation_SOA_t*, float*, int),
                          void (*molec_force_integration)(float*, float*, const float*, float*, const int),
                          void (*molec_periodic)(float*, const int));

#endif
