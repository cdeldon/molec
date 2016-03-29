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

#ifndef MOLEC_SORT_H
#define MOLEC_SORT_H

#include <molec/Common.h>
#include <molec/Simulation.h>

/**
  * @brief Structure that helps in defining a custom sorting algorithm
  *
  * @param key      contains the value of the particle to be compared with the other particles to 
  *                 be sorted
  * @param value    integer index used to sort other arrays accordingly
  */
typedef struct molec_Sort_Pair {
    Real key;
    int value;
} molec_Sort_Pair_t;


/**
  * @brief Compares two @c molec_Sort_Pair_t
  *
  * The return value is determined by:
  * -1	pair1 goes before pair2
  *  1	pair1 goes after pair2
  *
  */
int molec_compare(const void* pair1, const void* pair2);

/**
 * @brief Routine that sorts the particles in an increasing order according to the x component
 *
 * Routine that sorts the particles in an increasing order according to the x component
 * using the standard built in 'qsort' function *
 */
void molec_sort_qsort(molec_Simulation_SOA_t* sim);

#endif 
