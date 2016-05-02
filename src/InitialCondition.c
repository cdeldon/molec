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

#include <math.h>
#include <molec/InitialCondition.h>
#include <molec/Parameter.h>
#include <stdlib.h>

void molec_set_initial_condition(molec_Simulation_SOA_t* sim)
{
    const int N = molec_parameter->N;
    const float L_x = molec_parameter->L_x;
    const float L_y = molec_parameter->L_y;
    const float L_z = molec_parameter->L_z;
    const float scaling = molec_parameter->scaling;

    // Set initial positions (regular grid + random offset)
    float fraction = (L_x * L_x) / (L_y * L_z);

    const int num_atom_along_axis_x = (int) ceil(pow(fraction * N, 1./3));
    const int num_atom_along_axis_y = (int) ceil(num_atom_along_axis_x * L_y / L_x);
    const int num_atom_along_axis_z = (int) ceil(num_atom_along_axis_x * L_z / L_x);

    const float spread_x = L_x / num_atom_along_axis_x;
    const float spread_y = L_y / num_atom_along_axis_y;
    const float spread_z = L_z / num_atom_along_axis_z;
    const float dx = scaling * spread_x / 2.0;
    const float dy = scaling * spread_y / 2.0;
    const float dz = scaling * spread_z / 2.0;

    int atom_idx = 0;
    for(int iz = 0; iz < num_atom_along_axis_z; ++iz)
        for(int iy = 0; iy < num_atom_along_axis_y; ++iy)
            for(int ix = 0; ix < num_atom_along_axis_x && atom_idx < N; ++ix, ++atom_idx)
            {
                sim->x[atom_idx] = (0.5 + ix) * spread_x + dx * (2 * (float) rand() / RAND_MAX - 1);
                sim->y[atom_idx] = (0.5 + iy) * spread_y + dy * (2 * (float) rand() / RAND_MAX - 1);
                sim->z[atom_idx] = (0.5 + iz) * spread_z + dz * (2 * (float) rand() / RAND_MAX - 1);
            }

    // Set initial velocities
    for(int i = 0; i < N; ++i)
        sim->v_x[i] = sim->v_y[i] = sim->v_z[i] = 0.0;

    // Set initial forces
    for(int i = 0; i < N; ++i)
        sim->f_x[i] = sim->f_y[i] = sim->f_z[i] = 0.0;
}
