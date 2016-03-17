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

#include <molec/Force.h>
#include <molec/Parameter.h>

/**
 * Calculate distance between x and y taking periodic boundaries into account
 */
MOLEC_INLINE Real dist(Real x, Real y, Real L)
{
    Real r = x - y;
    if(r < -L / 2)
        r += L;
    else if(r > L / 2)
        r -= L;
    return r;
}

void molec_force_N2_refrence(molec_Simulation_SOA_t* sim, Real* Epot, const int N)
{
    assert(molec_parameter);
    const Real sigLJ = molec_parameter->sigLJ;
    const Real epsLJ = molec_parameter->epsLJ;
    const Real L = molec_parameter->L;
    const Real Rcut = molec_parameter->Rcut;

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

        for(int j = i + 1; j < N; ++j)
        {
            const Real xij = dist(xi, x[j], L);
            const Real yij = dist(yi, y[j], L);
            const Real zij = dist(zi, z[j], L);

            const Real r2 = xij * xij + yij * yij + zij * zij;

            if(r2 < Rcut * Rcut)
            {
                // V(s) = 4 * eps * (s^12 - s^6) with  s = sig/r
                const Real s2 = (sigLJ * sigLJ) / r2;
                const Real s6 = s2 * s2 * s2;

                Epot_ += 4 * epsLJ * (s6 * s6 - s6);

                const Real fr = 24 * epsLJ / r2 * (2 * s6 * s6 - s6);

                f_x[i] += fr * xij;
                f_y[i] += fr * yij;
                f_z[i] += fr * zij;

                f_x[j] -= fr * xij;
                f_y[j] -= fr * yij;
                f_z[j] -= fr * zij;
            }
        }
    }

    *Epot = Epot_;
}

