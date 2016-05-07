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
#include <molec/CellVector.h>
#include <string.h>
#include <math.h>

/**
 * Calculate distance between x and y taking periodic boundaries into account
 */
MOLEC_INLINE float distL2(float x, float y, float L, float L2)
{
    float r = x - y;
    if(r < -L2)
        r += L;
    else if(r > L2)
        r -= L;
    return r;
}

/**
 * Calculate the positive modulo between two integers, used for periodic BC
 */
MOLEC_INLINE int mod(int b, int m)
{
    return (b + m) % m;
}


/**
 * @brief neighbor_cells[i] constains the 27 indices of its neighbors
 */
static int ** neighbor_cells;

void molec_force_quadrant(molec_Simulation_SOA_t* sim, float* Epot, const int N)
{
    assert(molec_parameter);
    const float sigLJ = molec_parameter->sigLJ;
    const float epsLJ = molec_parameter->epsLJ;
    const float L_x = molec_parameter->L_x;
    const float L_y = molec_parameter->L_y;
    const float L_z = molec_parameter->L_z;
    const float L_x2 = L_x * 0.5f;
    const float L_y2 = L_y * 0.5f;
    const float L_z2 = L_z * 0.5f;
    const float Rcut2 = molec_parameter->Rcut2;

    molec_CellList_Parameter_t cellList_parameter = molec_parameter->cellList;

    // Local aliases
    const float* x = sim->x;
    const float* y = sim->y;
    const float* z = sim->z;
    float* f_x = sim->f_x;
    float* f_y = sim->f_y;
    float* f_z = sim->f_z;

    float Epot_ = 0;

    // Build the neighbor_cell array only if not initialized before
    if(neighbor_cells == NULL)
    {
        MOLEC_MALLOC(neighbor_cells, cellList_parameter.N * sizeof(int *));
        for(int i = 0; i < cellList_parameter.N; ++i)
            MOLEC_MALLOC(neighbor_cells[i], 27 * sizeof(int));

        // build the cell-neighborhood
        for(int idx_z = 0; idx_z < cellList_parameter.N_z; ++idx_z)
            for(int idx_y = 0; idx_y < cellList_parameter.N_y; ++idx_y)
                for(int idx_x = 0; idx_x < cellList_parameter.N_x; ++idx_x)
                {
                    // compute scalar cell index
                    const int idx = idx_x
                                    + cellList_parameter.N_x * (idx_y + cellList_parameter.N_y * idx_z);

                    int neighbor_number = 0;
                    // loop over neighbor cells
                    for(int d_z = -1; d_z <= 1; ++d_z)
                        for(int d_y = -1; d_y <= 1; ++d_y)
                            for(int d_x = -1; d_x <= 1; ++d_x)
                            {
                                // compute cell index considering periodic BC
                                int n_idx_z = mod(idx_z + d_z, cellList_parameter.N_z);
                                int n_idx_y = mod(idx_y + d_y, cellList_parameter.N_y);
                                int n_idx_x = mod(idx_x + d_x, cellList_parameter.N_x);

                                // linear index
                                int n_idx = n_idx_x
                                            + cellList_parameter.N_x
                                                  * (n_idx_y + cellList_parameter.N_y * n_idx_z);

                                // store the neighbor index
                                neighbor_cells[idx][neighbor_number++] = n_idx;
                            }
                }
    }

    // for each particle compute cell index
    int* c_idx;
    MOLEC_MALLOC(c_idx, sizeof(int) * N);

    /* Note that the cell list will be traversed in the following form:
     *    for(z=0;z<cellList_parameter.N_z;++z)
     *      for(y=0;y<cellList_parameter.N_y;++y)
     *          for(x=0;x<cellList_parameter.N_x;++x)
     *
     * the fastest running index is x, while the slowest is z
     */
    for(int i = 0; i < N; ++i)
    {
        // linear one dimensional index of cell associated to i-th particle
        int idx_x = x[i] / cellList_parameter.c_x;
        int idx_y = y[i] / cellList_parameter.c_y;
        int idx_z = z[i] / cellList_parameter.c_z;

        // linear index of cell
        int idx = idx_x + cellList_parameter.N_x * (idx_y + cellList_parameter.N_y * idx_z);
        c_idx[i] = idx;
    }

    // count number of particles in each cell
    int* particles_in_cell_idx;
    MOLEC_MALLOC(particles_in_cell_idx, sizeof(int) * cellList_parameter.N);
    for(int idx = 0; idx < cellList_parameter.N; ++idx)
        particles_in_cell_idx[idx] = 0;

    for(int i = 0; i < N; ++i)
        particles_in_cell_idx[c_idx[i]] += 1;

    // count max number of particles per cell
    int max_particles_per_cell = 0;
    for(int idx = 0; idx < cellList_parameter.N; ++idx)
        max_particles_per_cell = fmax(max_particles_per_cell, particles_in_cell_idx[idx]);

    // generate cell list, cellList[idx][k] is the k-th particle of cell idx
    int** cellList;
    MOLEC_MALLOC(cellList, sizeof(int*) * cellList_parameter.N);
    for(int idx = 0; idx < cellList_parameter.N; ++idx)
        MOLEC_MALLOC(cellList[idx], sizeof(int) * max_particles_per_cell);

    // set of counters for next for loop
    int* cell_counter;
    MOLEC_MALLOC(cell_counter, sizeof(int) * cellList_parameter.N);
    for(int idx = 0; idx < cellList_parameter.N; ++idx)
        cell_counter[idx] = 0;

    for(int i = 0; i < N; ++i)
    {
        int idx = c_idx[i];
        cellList[idx][cell_counter[idx]++] = i;
    }

    // Reset forces
    for(int i = 0; i < N; ++i)
        f_x[i] = f_y[i] = f_z[i] = 0.0;

    // Loop over the cells
    for(int idx_z = 0; idx_z < cellList_parameter.N_z; ++idx_z)
        for(int idx_y = 0; idx_y < cellList_parameter.N_y; ++idx_y)
            for(int idx_x = 0; idx_x < cellList_parameter.N_x; ++idx_x)
            {
                // compute scalar cell index
                const int idx = idx_x
                                + cellList_parameter.N_x * (idx_y + cellList_parameter.N_y * idx_z);

                int* particles_idx = cellList[idx];

                // loop over neighbor cells
                for(int neighbor_number = 0; neighbor_number < 27; ++neighbor_number)
                {
                    int n_idx = neighbor_cells[idx][neighbor_number];

                    // only compute interaction if neighbor cell n_idx is smaller than
                    // the current cell idx
                    if(idx > n_idx)
                    {
                        int* particles_n_idx = cellList[n_idx];

                        // all particles need to interact
                        int num_particles_in_cell_idx = particles_in_cell_idx[idx];
                        int num_particles_in_cell_n_idx = particles_in_cell_idx[n_idx];
                        for(int k_idx = 0; k_idx < num_particles_in_cell_idx; ++k_idx)
                        {
                            int i = particles_idx[k_idx];

                            // local aliases for particle i in cell idx
                            const float xi = x[i];
                            const float yi = y[i];
                            const float zi = z[i];

                            float f_xi = f_x[i];
                            float f_yi = f_y[i];
                            float f_zi = f_z[i];

                            for(int k_n_idx = 0; k_n_idx < num_particles_in_cell_n_idx;
                                ++k_n_idx)
                            {
                                int j = particles_n_idx[k_n_idx];

                                const float xij = distL2(xi, x[j], L_x, L_x2);
                                const float yij = distL2(yi, y[j], L_y, L_y2);
                                const float zij = distL2(zi, z[j], L_z, L_z2);

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
                    }
                    // handle case of same cell
                    else if(idx == n_idx)
                    {
                        int* particles_idx = cellList[idx];

                        // all particles need to interact
                        int num_particles_in_cell_idx = particles_in_cell_idx[idx];

                        for(int k_idx1 = 0; k_idx1 < num_particles_in_cell_idx; ++k_idx1)
                        {
                            int i = particles_idx[k_idx1];

                            // local aliases for particle i in cell idx
                            const float xi = x[i];
                            const float yi = y[i];
                            const float zi = z[i];

                            float f_xi = f_x[i];
                            float f_yi = f_y[i];
                            float f_zi = f_z[i];

                            for(int k_idx2 = k_idx1 + 1; k_idx2 < num_particles_in_cell_idx;
                                ++k_idx2)
                            {
                                int j = particles_idx[k_idx2];

                                const float xij = distL2(xi, x[j], L_x, L_x2);
                                const float yij = distL2(yi, y[j], L_y, L_y2);
                                const float zij = distL2(zi, z[j], L_z, L_z2);

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
                    } // end idx == n_idx
                } // end loop over neighbor cells

            } // end loop over cells idx

    // free memory
    MOLEC_FREE(c_idx);
    MOLEC_FREE(particles_in_cell_idx);
    for(int idx = 0; idx < cellList_parameter.N; ++idx)
    {
        int* ptr = cellList[idx];
        MOLEC_FREE(ptr);
    }

    MOLEC_FREE(cellList);
    MOLEC_FREE(cell_counter);

    *Epot = Epot_;
}
