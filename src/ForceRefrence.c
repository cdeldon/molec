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

#include <molec/Force.h>
#include <molec/Parameter.h>
#include <molec/Inline.h>


void molec_force_N2_refrence(molec_Simulation_SOA_t* sim, Real* Epot, const int N)
{
    const Real L = molec_parameter->L;
    // Local aliases
    const Real* x = sim->x;
    const Real* y = sim->y;
    const Real* z = sim->z;
    Real* f_x = sim->f_x;
    Real* f_y = sim->f_y;
    Real* f_z = sim->f_z;

    Real Epot_ = 0;

    // Reset forces
    for(int i = 0; i < N; ++i)
        f_x[i] = f_y[i] = f_z[i] = 0.0;

    for(int i = 0; i < N; ++i)
    {

        const Real xi = x[i];
        const Real yi = y[i];
        const Real zi = z[i];

        Real f_xi = f_x[i];
        Real f_yi = f_y[i];
        Real f_zi = f_z[i];

        for(int j = i + 1; j < N; ++j)
        {
        const Real xij = dist(xi, x[j], L);
        const Real yij = dist(yi, y[j], L);
        const Real zij = dist(zi, z[j], L);

        Real fr=update(xij,yij,zij,&f_xi,&f_yi,&f_zi,&Epot_);

        f_x[j] -= fr * xij;
        f_y[j] -= fr * yij;
        f_z[j] -= fr * zij;

        }

        f_x[i] = f_xi;
        f_y[i] = f_yi;
        f_z[i] = f_zi;
    }

    *Epot = Epot_;
}
