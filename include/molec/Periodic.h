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

#ifndef MOLEC_PERIODIC_H
#define MOLEC_PERIODIC_H

#include <molec/Common.h>

/** Define function pointer for periodic routine */
typedef void (*molec_periodic)(float*, const int, const float);

/**
 * @brief Apply periodic boundary conditions to the given position array
 *
 * After the operation every position in @c x will be between [0, L_i] where L_i is the extend of the
 * simulation along the i-th dimension (x,y,z) as given in molec_Parameter.
 *
 * @param x     Array of lenght N representing the position
 * @param N     Size of arrays
 * @param L     Size of bounding box in the corresponding dimension
 */
void molec_periodic_refrence(float* x, const int N, const float L);
void molec_periodic_close4(float* x, const int N, const float L);

#endif
