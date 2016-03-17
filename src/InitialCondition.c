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

#include <molec/InitialCondition.h>
#include <molec/Parameter.h>
#include <stdlib.h>
#include <math.h>

void molec_set_initial_condition(molec_Simulation_SOA_t* sim)
{
    const int N = molec_parameter->N;
    const Real L = molec_parameter->L;
    const Real scaling = molec_parameter->scaling;

    // Set initial positions (regular grid + random offset)
    const int num_atom_along_axis = (int) ceil(pow(N, 1.0 / 3.0));
    const Real spread = L / num_atom_along_axis;
    const Real dx = scaling * spread / 2.0;

    int atom_idx = 0;
    for(int iz = 0; iz < num_atom_along_axis; ++iz)
        for(int iy = 0; iy < num_atom_along_axis; ++iy)
            for(int ix = 0; ix < num_atom_along_axis && atom_idx < N; ++ix, ++atom_idx)
            {
                sim->x[atom_idx] = (0.5 + ix) * spread
                                          + dx * (2 * (Real) rand() / RAND_MAX - 1);
                sim->y[atom_idx] = (0.5 + iy) * spread
                                          + dx * (2 * (Real) rand() / RAND_MAX - 1);
                sim->z[atom_idx] = (0.5 + iz) * spread
                                          + dx * (2 * (Real) rand() / RAND_MAX - 1);
            }

    // Set initial velocities
    for(int i = 0; i < N; ++i)
        sim->v_x[i] = sim->v_y[i] = sim->v_z[i] = 0.0;

    // Set initial forces
    for(int i = 0; i < N; ++i)
        sim->f_x[i] = sim->f_y[i] = sim->f_z[i] = 0.0;
}

