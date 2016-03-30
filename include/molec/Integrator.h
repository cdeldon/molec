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

#ifndef MOLEC_INTEGRATOR_H
#define MOLEC_INTEGRATOR_H

#include <molec/Common.h>

/**
 * @brief Leapfrog integration scheme (refrence version)
 *
 * Integrate the position and velocity array @c x and @c v using a leapfrog integration scheme.
 * In addition the kinect energy @c Ekin is going to be calculated.
 * @see https://en.wikipedia.org/wiki/Leapfrog_integration
 *
 * @param x     Array of lenght N representing the position (stored at t_i)
 * @param v     Array of lenght N representing the velocity (stored at t_{i+1/2})
 * @param f     Array of lenght N representing the force (stored at t_i)
 * @param Ekin  Real scalar to store the kinetic energy
 * @param N     Size of arrays
 */
void molec_integrator_leapfrog_refrence(Real* x, Real* v, const Real* f, Real* Ekin, const int N);
void molec_integrator_leapfrog_unroll_2(Real* x, Real* v, const Real* f, Real* Ekin, const int N);

#endif

