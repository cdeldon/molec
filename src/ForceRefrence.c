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

/**
 * Calculate distance between x and y taking periodic boundaries into account
 */
MOLEC_INLINE float dist(float x, float y, float L)
{
    float r = x - y;
    if(r < -L / 2)
        r += L;
    else if(r > L / 2)
        r -= L;
    return r;
}

void molec_force_N2_refrence(molec_Simulation_SOA_t* sim, float* Epot, const int N)
{
    assert(molec_parameter);
    const float sigLJ = molec_parameter->sigLJ;
    const float epsLJ = molec_parameter->epsLJ;
    const float L = molec_parameter->L;
    const float Rcut2 = molec_parameter->Rcut2;

    // Local aliases
    const float* x = sim->x;
    const float* y = sim->y;
    const float* z = sim->z;
    float* f_x = sim->f_x;
    float* f_y = sim->f_y;
    float* f_z = sim->f_z;

    float Epot_ = 0;

    // Reset forces
    for(int i = 0; i < N; ++i)
        f_x[i] = f_y[i] = f_z[i] = 0.0;

    for(int i = 0; i < N; ++i)
    {
        const float xi = x[i];
        const float yi = y[i];
        const float zi = z[i];

        float f_xi = f_x[i];
        float f_yi = f_y[i];
        float f_zi = f_z[i];

        for(int j = i + 1; j < N; ++j)
        {
            const float xij = dist(xi, x[j], L);
            const float yij = dist(yi, y[j], L);
            const float zij = dist(zi, z[j], L);

            const float r2 = xij * xij + yij * yij + zij * zij;

            if(r2 < Rcut2)
            {
                // V(s) = 4 * eps * (s^12 - s^6) with  s = sig/r
                const float s2 = (sigLJ * sigLJ) / r2;
                const float s6 = s2 * s2 * s2;

                Epot_ += 4 * epsLJ * (s6 * s6 - s6);

                const float fr = 24 * epsLJ / r2 * (2 * s6 * s6 - s6);

                f_xi += fr * xij;
                f_yi += fr * yij;
                f_zi += fr * zij;

                f_x[j] -= fr * xij;
                f_y[j] -= fr * yij;
                f_z[j] -= fr * zij;
            }
        }

        f_x[i] = f_xi;
        f_y[i] = f_yi;
        f_z[i] = f_zi;
    }

    *Epot = Epot_;
}
