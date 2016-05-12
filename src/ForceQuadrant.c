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
#include <molec/Quadrant.h>
#include <molec/Parameter.h>

#include <string.h>

/**
 * Calculate distance between x and y taking periodic boundaries into account
 *
 * @param      x  coordinate of first particle
 * @param      y  coordinate of second particle
 * @param      L  bounding box size
 * @param      L  bounding box size * 0.5
 *
 * @return     distance between two particles taking into
 *             account periodicity of bounding box
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
 * @brief neighbor_cells[i] constains the 27 indices of its neighbor cells
 */
static int** neighbor_cells;

void molec_force_quadrant(molec_Simulation_SOA_t* sim, float* Epot, const int N)
{
    assert(molec_parameter);

    // Local aliases

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

    float* x = sim->x;
    float* y = sim->y;
    float* z = sim->z;
    float* v_x = sim->v_x;
    float* v_y = sim->v_y;
    float* v_z = sim->v_z;
    float* f_x = sim->f_x;
    float* f_y = sim->f_y;
    float* f_z = sim->f_z;

    float Epot_ = 0;

    // Build the neighbor_cell array only if not initialized before
    if(neighbor_cells == NULL)
    {
        MOLEC_MALLOC(neighbor_cells, cellList_parameter.N * sizeof(int*));
        for(int i = 0; i < cellList_parameter.N; ++i)
            MOLEC_MALLOC(neighbor_cells[i], 27 * sizeof(int));

        // build the cell-neighborhood
        for(int idx_z = 0; idx_z < cellList_parameter.N_z; ++idx_z)
            for(int idx_y = 0; idx_y < cellList_parameter.N_y; ++idx_y)
                for(int idx_x = 0; idx_x < cellList_parameter.N_x; ++idx_x)
                {
                    // compute scalar cell index
                    const int idx
                        = idx_x + cellList_parameter.N_x * (idx_y + cellList_parameter.N_y * idx_z);

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

    molec_Quadrant_t* quadrants = malloc(cellList_parameter.N * sizeof(molec_Quadrant_t));

    for (int i = 0; i < cellList_parameter.N; ++i)
        quadrants[i].N = 0;


    // for each particle compute cell index and count number of particles inside each cell
    int* c_idx = malloc(sizeof(int) * N);

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
        quadrants[idx].N += 1;
    }

    for (int i = 0; i < cellList_parameter.N; ++i)
    {
        MOLEC_MALLOC(quadrants[i].x, quadrants[i].N * sizeof(float));
        MOLEC_MALLOC(quadrants[i].y, quadrants[i].N * sizeof(float));
        MOLEC_MALLOC(quadrants[i].z, quadrants[i].N * sizeof(float));

        MOLEC_MALLOC(quadrants[i].v_x, quadrants[i].N * sizeof(float));
        MOLEC_MALLOC(quadrants[i].v_y, quadrants[i].N * sizeof(float));
        MOLEC_MALLOC(quadrants[i].v_z, quadrants[i].N * sizeof(float));

        MOLEC_MALLOC(quadrants[i].f_x, quadrants[i].N * sizeof(float));
        MOLEC_MALLOC(quadrants[i].f_y, quadrants[i].N * sizeof(float));
        MOLEC_MALLOC(quadrants[i].f_z, quadrants[i].N * sizeof(float));
    }

    // array that counts in which position store the new particles (initialized to 0)
    // and copy the particles into the quadrants
    int *current_particle_number_of_quadrant = calloc(cellList_parameter.N, sizeof(int));
    for (int i = 0; i < N; ++i)
    {
        int idx = c_idx[i];
        int pos = current_particle_number_of_quadrant[idx];

        quadrants[idx].x[pos] = x[i];
        quadrants[idx].y[pos] = y[i];
        quadrants[idx].z[pos] = z[i];

        quadrants[idx].v_x[pos] = v_x[i];
        quadrants[idx].v_y[pos] = v_y[i];
        quadrants[idx].v_z[pos] = v_z[i];

        quadrants[idx].f_x[pos] = 0.f;
        quadrants[idx].f_y[pos] = 0.f;
        quadrants[idx].f_z[pos] = 0.f;

        current_particle_number_of_quadrant[idx] += 1;
    }
    free(current_particle_number_of_quadrant);


    // loop over the quadrants
    for (int idx = 0; idx < cellList_parameter.N; ++idx)
    {
        // load the quadrant associated to cell idx
        molec_Quadrant_t q_idx = quadrants[idx];
        int N_idx = q_idx.N;

        for (int neighbor_cell = 0; neighbor_cell < 27; ++neighbor_cell)
        {
            int n_idx = neighbor_cells[idx][neighbor_cell];

            if (idx > n_idx)
            {
                // load the quadrant associated to cell n_idx
                molec_Quadrant_t q_n_idx = quadrants[n_idx];
                int N_n_idx = q_n_idx.N;

                for(int i = 0; i < N_idx; ++i)
                {
                    float xi = q_idx.x[i];
                    float yi = q_idx.y[i];
                    float zi = q_idx.z[i];

                    float f_xi = q_idx.f_x[i];
                    float f_yi = q_idx.f_y[i];
                    float f_zi = q_idx.f_z[i];

                    for(int j = 0; j < N_n_idx; ++j)
                    {
                        const float xij = distL2(xi, q_n_idx.x[j], L_x, L_x2);
                        const float yij = distL2(yi, q_n_idx.y[j], L_y, L_y2);
                        const float zij = distL2(zi, q_n_idx.z[j], L_z, L_z2);

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

                            q_n_idx.f_x[j] -= fr * xij;
                            q_n_idx.f_y[j] -= fr * yij;
                            q_n_idx.f_z[j] -= fr * zij;
                        }

                    }
                    q_idx.f_x[i] = f_xi;
                    q_idx.f_y[i] = f_yi;
                    q_idx.f_z[i] = f_zi;
                }
                quadrants[n_idx] = q_n_idx;
            }
            else if (idx == n_idx)
            {
                for(int i = 0; i < N_idx; ++i)
                {
                    float xi = q_idx.x[i];
                    float yi = q_idx.y[i];
                    float zi = q_idx.z[i];

                    float f_xi = q_idx.f_x[i];
                    float f_yi = q_idx.f_y[i];
                    float f_zi = q_idx.f_z[i];

                    for(int j = i+1; j < N_idx; ++j)
                    {

                        const float xij = distL2(xi, q_idx.x[j], L_x, L_x2);
                        const float yij = distL2(yi, q_idx.y[j], L_y, L_y2);
                        const float zij = distL2(zi, q_idx.z[j], L_z, L_z2);

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

                            q_idx.f_x[j] -= fr * xij;
                            q_idx.f_y[j] -= fr * yij;
                            q_idx.f_z[j] -= fr * zij;
                        }

                    }
                    q_idx.f_x[i] = f_xi;
                    q_idx.f_y[i] = f_yi;
                    q_idx.f_z[i] = f_zi;
                }
            } // end idx == n_idx
        } // end loop over neighbor cells
        quadrants[idx] = q_idx;
    } // end loop over quadrants


    // copy back data from quadrants to 1D arrays
    int n_1D = 0;
    for (int idx = 0; idx < cellList_parameter.N; ++idx)
    {
        molec_Quadrant_t q_idx = quadrants[idx];
        int N_idx = q_idx.N;

        memcpy(x + n_1D, q_idx.x, N_idx * sizeof(float));
        memcpy(y + n_1D, q_idx.y, N_idx * sizeof(float));
        memcpy(z + n_1D, q_idx.z, N_idx * sizeof(float));

        memcpy(v_x + n_1D, q_idx.v_x, N_idx * sizeof(float));
        memcpy(v_y + n_1D, q_idx.v_y, N_idx * sizeof(float));
        memcpy(v_z + n_1D, q_idx.v_z, N_idx * sizeof(float));

        memcpy(f_x + n_1D, q_idx.f_x, N_idx * sizeof(float));
        memcpy(f_y + n_1D, q_idx.f_y, N_idx * sizeof(float));
        memcpy(f_z + n_1D, q_idx.f_z, N_idx * sizeof(float));

        n_1D += N_idx;
    }


    // free memory
    free(c_idx);

    for (int i = 0; i < cellList_parameter.N; ++i)
    {
        MOLEC_FREE(quadrants[i].x);
        MOLEC_FREE(quadrants[i].y);
        MOLEC_FREE(quadrants[i].z);
        MOLEC_FREE(quadrants[i].v_x);
        MOLEC_FREE(quadrants[i].v_y);
        MOLEC_FREE(quadrants[i].v_z);
        MOLEC_FREE(quadrants[i].f_x);
        MOLEC_FREE(quadrants[i].f_y);
        MOLEC_FREE(quadrants[i].f_z);
    }

    free(quadrants);

    *Epot = Epot_;
}
