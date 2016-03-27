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
 * @brief Routine that sorts the particles in an increasing order according to the x component
 *
 * Routine that sorts the particles in an increasing order according to the x component
 * using the standard built in 'qsort' function *
 */
void molec_sort_qsort(molec_Simulation_SOA_t* sim);

#endif // MOLEC_SORT_H
