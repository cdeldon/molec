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

#include <molec/Quadrant.h>

void molec_build_cell_neighbors(int** neighbor_cells,
                                molec_CellList_Parameter_t cellList_parameter)
{
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


molec_Quadrant_t* molec_quadrant_init(const int N,
                                      molec_CellList_Parameter_t cellList_parameter,
                                      molec_Simulation_SOA_t* sim)
{
    const float* x = sim->x;
    const float* y = sim->y;
    const float* z = sim->z;

    const float* v_x = sim->v_x;
    const float* v_y = sim->v_y;
    const float* v_z = sim->v_z;


    molec_Quadrant_t* quadrants = malloc(cellList_parameter.N * sizeof(molec_Quadrant_t));

    for(int i = 0; i < cellList_parameter.N; ++i)
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

    for(int i = 0; i < cellList_parameter.N; ++i)
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
    int* current_particle_number_of_quadrant = calloc(cellList_parameter.N, sizeof(int));

    for(int i = 0; i < N; ++i)
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

    free(c_idx);
    free(current_particle_number_of_quadrant);

    return quadrants;
}


void molec_quadrants_finalize(molec_Quadrant_t* quadrants,
                         molec_CellList_Parameter_t cellList_parameter,
                         molec_Simulation_SOA_t* sim)
{
    // copy back data from quadrants to 1D arrays
    int n_1D = 0;
    for(int idx = 0; idx < cellList_parameter.N; ++idx)
    {
        molec_Quadrant_t q_idx = quadrants[idx];
        int N_idx = q_idx.N;

        memcpy(sim->x + n_1D, q_idx.x, N_idx * sizeof(float));
        memcpy(sim->y + n_1D, q_idx.y, N_idx * sizeof(float));
        memcpy(sim->z + n_1D, q_idx.z, N_idx * sizeof(float));

        memcpy(sim->v_x + n_1D, q_idx.v_x, N_idx * sizeof(float));
        memcpy(sim->v_y + n_1D, q_idx.v_y, N_idx * sizeof(float));
        memcpy(sim->v_z + n_1D, q_idx.v_z, N_idx * sizeof(float));

        memcpy(sim->f_x + n_1D, q_idx.f_x, N_idx * sizeof(float));
        memcpy(sim->f_y + n_1D, q_idx.f_y, N_idx * sizeof(float));
        memcpy(sim->f_z + n_1D, q_idx.f_z, N_idx * sizeof(float));

        n_1D += N_idx;
    }


    // free memory

    for(int i = 0; i < cellList_parameter.N; ++i)
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
}
