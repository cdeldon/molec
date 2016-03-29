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


// Define function pointer for force calculation routine
typedef void (*molec_force_calculation)(molec_Simulation_SOA_t*, Real*, const int);

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
void molec_force_N2_refrence(molec_Simulation_SOA_t* sim, Real* Epot, const int N);

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
void molec_force_cellList(molec_Simulation_SOA_t* sim, Real* Epot, const int N);

/**
 * @brief Calculate short-range interaction force using cell-list approach [1]
 *
 * This function computes the force between all the particiles exploting the short
 * range interaction form of the Lennard-Jones potential using a cell list
 * data structure which allows neighbourhood queries in constant time.
 *
 * From the implementation point of view, this function is more prone to
 * be adapted for the Gonnet algorithm (@c molec_force_gonnet) as the
 * cell list data structure is transformed in local vectors at each cell iteration
 *
 * @see [1] http://cacs.usc.edu/education/cs596/01-1LinkedListCell.pdf
 *
 * @param sim   Simulation holding the position, velocity and force arrays
 * @param Epot  Real scalar to store the potential energy
 * @param N     Size of arrays
 */
void molec_force_cellList_for(molec_Simulation_SOA_t* sim, Real* Epot, const int N);

/**
 * @brief Calculate short-range interaction force using Gonnet
 *
 * @param sim   Simulation holding the position, velocity and force arrays
 * @param Epot  Real scalar to store the potential energy
 * @param N     Size of arrays
 */
void molec_force_gonnet(molec_Simulation_SOA_t* sim, Real* Epot, const int N);

#endif

