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

#ifndef MOLEC_FORCE_H
#define MOLEC_FORCE_H

#include <molec/Common.h>
#include <molec/Simulation.h>

/** Define function pointer for force calculation routine */
typedef void (*molec_force_calculation)(molec_Simulation_SOA_t*, float*, const int);

/**
 * @brief Calculate short-range interaction force using N^2 approach
 *
 * This function computes the force between all the particiles using a straightforward N^2 
 * approach. In addition the potential energy is computed.
 * @see http://polymer.bu.edu/Wasser/robert/work/node8.html for derivation of the forces.
 *
 * @param sim   Simulation holding the position, velocity and force arrays
 * @param Epot  Real scalar to store the potential energy
 * @param N     Size of arrays
 */
void molec_force_N2_refrence(molec_Simulation_SOA_t* sim, float* Epot, const int N);

/**
 * @brief Calculate short-range interaction force using cell-list approach [1]
 *
 * This function computes the force between all the particiles exploting the short
 * range interaction form of the Lennard-Jones potential using a cell list
 * data structure which allows neighbourhood queries in constant time.
 *
 * @see [1] http://cacs.usc.edu/education/cs596/01-1LinkedListCell.pdf
 *
 * @param sim   Simulation holding the position, velocity and force arrays
 * @param Epot  Real scalar to store the potential energy
 * @param N     Size of arrays
 */
void molec_force_cellList_knuth(molec_Simulation_SOA_t* sim, float* Epot, const int N);

/**
 * @brief Calculate short-range interaction using table of particle indices per cell
 *
 * Cell list implementation using double pointers instead of a one dimensional array
 *
 * @param sim   Simulation holding the position, velocity and force arrays
 * @param Epot  Real scalar to store the potential energy
 * @param N     Size of arrays
 */
void molec_force_cellList_reference(molec_Simulation_SOA_t *sim, float *Epot, const int N);
void molec_force_cellList_v1(molec_Simulation_SOA_t *sim, float *Epot, const int N);
void molec_force_cellList_v2(molec_Simulation_SOA_t *sim, float *Epot, const int N);
void molec_force_quadrant(molec_Simulation_SOA_t *sim, float *Epot, const int N);
void molec_force_quadrant_avx(molec_Simulation_SOA_t* sim, float* Epot, const int N);

#if 0
/**
 * @brief Calculation short-range interaction with double index arrays
 *
 * @param sim   Simulation holding the position, velocity and force arrays
 * @param Epot  Real scalar to store the potential energy
 * @param N     Size of arrays
 */
void molec_force_cellList_double_pointer_v2(molec_Simulation_SOA_t* sim, float* Epot, const int N);

/**
 * @brief Calculate short-range interaction force using cell-list approach [1]
 *
 * the internal while loop is replaced by a for loop
 *
 * @see [1] http://cacs.usc.edu/education/cs596/01-1LinkedListCell.pdf
 *
 * @param sim   Simulation holding the position, velocity and force arrays
 * @param Epot  Real scalar to store the potential energy
 * @param N     Size of arrays
 */
void molec_force_cellList_for(molec_Simulation_SOA_t* sim, float* Epot, const int N);

/**
 * @brief Calculate short-range interaction force using cell-list approach [1]
 *
 * From the implementation point of view, this function is different from
 * @c molec_force_cellList_for as the position and velocity arrays have a local swapping copy
 * that is used in the sorting of the particles
 *
 * @see [1] http://cacs.usc.edu/education/cs596/01-1LinkedListCell.pdf
 *
 * @param sim   Simulation holding the position, velocity and force arrays
 * @param Epot  Real scalar to store the potential energy
 * @param N     Size of arrays
 */
void molec_force_cellList_for_swap(molec_Simulation_SOA_t* sim, float* Epot, const int N);
#endif

#endif

