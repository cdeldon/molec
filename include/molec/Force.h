/*                   _
 *   _ __ ___   ___ | | ___  ___
 *  | '_ ` _ \ / _ \| |/ _ \/ __|
 *  | | | | | | (_) | |  __/ (__
 *  |_| |_| |_|\___/|_|\___|\___| - Molecular Dynamics Framework
 *
 *  Copyright (C) 2016  Carlo Del Don  (deldonc@student.ethz.ch)
 *                      Michel Breyer  (mbreyer@student.ethz.ch)
 *                      Florian Frei   (frofrei@student.ethz.ch)
 *                      Fabian Thuring (thfabian@student.ethz.ch)
 *
 *  This file is distributed under the MIT Open Source License.
 *  See LICENSE.txt for details.
 */

#ifndef MOLEC_FORCE_H
#define MOLEC_FORCE_H

#include <molec/Common.h>
#include <molec/Simulation.h>

/**
 * @brief Calculate short-range interaction force using the Lenard-Jones potential
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
 * @brief Calculate short-range interaction force using the Lenard-Jones potential
 *
 * This function computes the force between all the particiles exploting the short
 * range interaction form of the Lennard-Jones potential, introducing a cell list
 * data structure which allows to query neighbourhood queries in constant time.
 * @see http://cacs.usc.edu/education/cs596/01-1LinkedListCell.pdf for an algorithmic overview
 *
 * @param sim   Simulation holding the position, velocity and force arrays
 * @param Epot  Real scalar to store the potential energy
 * @param N     Size of arrays
 */
void molec_force_cellList(molec_Simulation_SOA_t* sim, Real* Epot, const int N);

#endif

