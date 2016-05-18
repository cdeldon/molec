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

void molec_build_cell_neighbors(int** neighbor_cells, molec_CellList_Parameter_t cellList_parameter)
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

    // for each particle compute cell index and the size of each quadrant

    int* cell_idx = malloc(sizeof(int) * N);

    for(int i = 0; i < N; ++i)
    {
        // linear one dimensional index of cell associated to i-th particle
        int idx_x = x[i] / cellList_parameter.c_x;
        int idx_y = y[i] / cellList_parameter.c_y;
        int idx_z = z[i] / cellList_parameter.c_z;

        // linear index of cell
        int idx = idx_x + cellList_parameter.N_x * (idx_y + cellList_parameter.N_y * idx_z);

        cell_idx[i] = idx;
        quadrants[idx].N += 1;
    }

    for(int i = 0; i < cellList_parameter.N; ++i)
    {
        int pad = quadrants[i].N % 8;
        quadrants[i].N_pad = quadrants[i].N + pad;
    }

    // allocate memory knowing the size of each quadrant
    for(int i = 0; i < cellList_parameter.N; ++i)
    {
        MOLEC_MALLOC(quadrants[i].x, quadrants[i].N_pad * sizeof(float));
        MOLEC_MALLOC(quadrants[i].y, quadrants[i].N_pad * sizeof(float));
        MOLEC_MALLOC(quadrants[i].z, quadrants[i].N_pad * sizeof(float));

        MOLEC_MALLOC(quadrants[i].v_x, quadrants[i].N_pad * sizeof(float));
        MOLEC_MALLOC(quadrants[i].v_y, quadrants[i].N_pad * sizeof(float));
        MOLEC_MALLOC(quadrants[i].v_z, quadrants[i].N_pad * sizeof(float));

        MOLEC_MALLOC(quadrants[i].f_x, quadrants[i].N_pad * sizeof(float));
        MOLEC_MALLOC(quadrants[i].f_y, quadrants[i].N_pad * sizeof(float));
        MOLEC_MALLOC(quadrants[i].f_z, quadrants[i].N_pad * sizeof(float));
    }

    // for each particle copy position, velocity and force inside their corresponding quadrant

    // array to keep track of the index inside each quadrant where to put particles in
    int* current_particle_number_of_quadrant = calloc(cellList_parameter.N, sizeof(int));

    for(int i = 0; i < N; ++i)
    {
        int idx = cell_idx[i];
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

    float pad_value = -10 * molec_parameter->Rcut;
    for(int i = 0; i < cellList_parameter.N; ++i)
    {
        for(int j = quadrants[i].N; j < quadrants[i].N_pad; ++j)
        {
            quadrants[i].x[j] = pad_value;
            quadrants[i].y[j] = pad_value;
            quadrants[i].z[j] = pad_value;
        }
    }

    free(cell_idx);
    free(current_particle_number_of_quadrant);

    return quadrants;
}


void molec_quadrants_finalize(molec_Quadrant_t* quadrants,
                              molec_CellList_Parameter_t cellList_parameter,
                              molec_Simulation_SOA_t* sim)
{
    // copy data from quadrants to 1D arrays
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

/**************************************************************************************************/

void molec_build_cell_neighbors_ghost(int** neighbor_cells,
                                      molec_CellList_Parameter_t cellList_parameter)
{
    // Number of cells in each direction considering also ghost cells
    const int N_x_ghost = cellList_parameter.N_x + 2;
    const int N_y_ghost = cellList_parameter.N_y + 2;

    // Traverse internal cellList considering also ghost cells
    for(int idx_z = 1; idx_z <= cellList_parameter.N_z; ++idx_z)
        for(int idx_y = 1; idx_y <= cellList_parameter.N_y; ++idx_y)
            for(int idx_x = 1; idx_x <= cellList_parameter.N_x; ++idx_x)
            {
                // compute scalar cell index
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);

                int neighbor_number = 0;

                // loop over neighbor cells
                for(int d_z = -1; d_z <= 1; ++d_z)
                    for(int d_y = -1; d_y <= 1; ++d_y)
                        for(int d_x = -1; d_x <= 1; ++d_x)
                        {
                            // compute cell index (no need to consider BC as we looped over internal
                            // cells only)
                            int n_idx_z = idx_z + d_z;
                            int n_idx_y = idx_y + d_y;
                            int n_idx_x = idx_x + d_x;

                            // linear index
                            int n_idx = n_idx_x + N_x_ghost * (n_idx_y + N_y_ghost * n_idx_z);

                            // store the neighbor index
                            neighbor_cells[idx][neighbor_number++] = n_idx;
                        }
            }
}


molec_Quadrant_t* molec_quadrant_init_ghost(const int N,
                                            molec_CellList_Parameter_t cellList_parameter,
                                            molec_Simulation_SOA_t* sim)
{
    const float* x = sim->x;
    const float* y = sim->y;
    const float* z = sim->z;

    const float* v_x = sim->v_x;
    const float* v_y = sim->v_y;
    const float* v_z = sim->v_z;

    // Number of cells in each direction considering also ghost cells
    const int N_x_ghost = cellList_parameter.N_x + 2;
    const int N_y_ghost = cellList_parameter.N_y + 2;
    const int N_z_ghost = cellList_parameter.N_z + 2;

    const int N_ghost = N_x_ghost * N_y_ghost * N_z_ghost;

    molec_Quadrant_t* quadrants = malloc(N_ghost * sizeof(molec_Quadrant_t));

    for(int i = 0; i < N_ghost; ++i)
        quadrants[i].N = 0;


    // for each particle compute cell index and the size of each quadrant considering also ghost
    // cells in the indices
    int* cell_idx = malloc(sizeof(int) * N);

    for(int i = 0; i < N; ++i)
    {
        // linear one dimensional index of cell associated to i-th particle considering boundary
        // ghost quadrants
        int idx_x = 1 + x[i] / cellList_parameter.c_x;
        int idx_y = 1 + y[i] / cellList_parameter.c_y;
        int idx_z = 1 + z[i] / cellList_parameter.c_z;

        // linear index of cell
        int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);

        cell_idx[i] = idx;
        quadrants[idx].N += 1;
    }

    for(int i = 0; i < N_ghost; ++i)
    {
        int pad = quadrants[i].N % 8;
        quadrants[i].N_pad = quadrants[i].N + pad;
    }

    // allocate memory knowing the size of each internal quadrant
    for(int idx_z = 1; idx_z <= cellList_parameter.N_z; ++idx_z)
        for(int idx_y = 1; idx_y <= cellList_parameter.N_y; ++idx_y)
            for(int idx_x = 1; idx_x <= cellList_parameter.N_x; ++idx_x)
            {
                // linear index of cell
                int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
                quadrants[idx].idx = idx;

                MOLEC_MALLOC(quadrants[idx].x, quadrants[idx].N_pad * sizeof(float));
                MOLEC_MALLOC(quadrants[idx].y, quadrants[idx].N_pad * sizeof(float));
                MOLEC_MALLOC(quadrants[idx].z, quadrants[idx].N_pad * sizeof(float));

                MOLEC_MALLOC(quadrants[idx].v_x, quadrants[idx].N_pad * sizeof(float));
                MOLEC_MALLOC(quadrants[idx].v_y, quadrants[idx].N_pad * sizeof(float));
                MOLEC_MALLOC(quadrants[idx].v_z, quadrants[idx].N_pad * sizeof(float));

                MOLEC_MALLOC(quadrants[idx].f_x, quadrants[idx].N_pad * sizeof(float));
                MOLEC_MALLOC(quadrants[idx].f_y, quadrants[idx].N_pad * sizeof(float));
                MOLEC_MALLOC(quadrants[idx].f_z, quadrants[idx].N_pad * sizeof(float));
            }

    // for each particle copy position, velocity and force inside their corresponding quadrant

    // array to keep track of the index inside each quadrant where to put particles in
    int* current_particle_number_of_quadrant = calloc(N_ghost, sizeof(int));

    for(int i = 0; i < N; ++i)
    {
        int idx = cell_idx[i];
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

    // perform linking between ghost cells and original inner cells by linking pointers between them
    int idx_x, idx_y, idx_z;       // indices of ghost cell beeing updated
    int idx_x_m, idx_y_m, idx_z_m; // indices of mirror cell


    // internal 'flat' ghost quadrants shrinked-by-1 faces
    {
        // idx_z = 0 --> mirror original quadrant: idx_z_m = cellList_parameter.N_z
        idx_z = 0;
        idx_z_m = cellList_parameter.N_z;
        for(idx_y = 1; idx_y <= cellList_parameter.N_y; ++idx_y)
            for(idx_x = 1; idx_x <= cellList_parameter.N_x; ++idx_x)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
                // linear index of mirror cell
                idx_x_m = idx_x;
                idx_y_m = idx_y;
                const int idx_m = idx_x_m + N_x_ghost * (idx_y_m + N_y_ghost * idx_z_m);

                quadrants[idx].N = quadrants[idx_m].N;
                quadrants[idx].N_pad = quadrants[idx_m].N_pad;

                // associate the index of the ghost quadrant to the real index of the mirror
                quadrants[idx].idx = quadrants[idx_m].idx;

                quadrants[idx].f_x = quadrants[idx_m].f_x;
                quadrants[idx].f_y = quadrants[idx_m].f_y;
                quadrants[idx].f_z = quadrants[idx_m].f_z;

                quadrants[idx].x = quadrants[idx_m].x;
                quadrants[idx].y = quadrants[idx_m].y;

                // copy-shift the z coordinats of the mirror cell by -L_z
                MOLEC_MALLOC(quadrants[idx].z, quadrants[idx_m].N_pad * sizeof(float));
                for(int k = 0; k < quadrants[idx_m].N_pad; ++k)
                    quadrants[idx].z[k] = quadrants[idx_m].z[k] - molec_parameter->L_z;
            }

        // idx_z = cellList_parameter.N_z + 1 --> mirror original quadrant: idx_z_m = 1
        idx_z = cellList_parameter.N_z + 1;
        idx_z_m = 1;
        for(idx_y = 1; idx_y <= cellList_parameter.N_y; ++idx_y)
            for(idx_x = 1; idx_x <= cellList_parameter.N_x; ++idx_x)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
                // linear index of mirror cell
                idx_x_m = idx_x;
                idx_y_m = idx_y;
                const int idx_m = idx_x_m + N_x_ghost * (idx_y_m + N_y_ghost * idx_z_m);

                quadrants[idx].N = quadrants[idx_m].N;
                                quadrants[idx].N_pad = quadrants[idx_m].N_pad;

                // associate the index of the ghost quadrant to the real index of the mirror
                quadrants[idx].idx = quadrants[idx_m].idx;

                quadrants[idx].f_x = quadrants[idx_m].f_x;
                quadrants[idx].f_y = quadrants[idx_m].f_y;
                quadrants[idx].f_z = quadrants[idx_m].f_z;

                quadrants[idx].x = quadrants[idx_m].x;
                quadrants[idx].y = quadrants[idx_m].y;

                // copy-shift the z coordinats of the mirror cell by +L_z
                MOLEC_MALLOC(quadrants[idx].z, quadrants[idx_m].N_pad * sizeof(float));
                for(int k = 0; k < quadrants[idx_m].N_pad; ++k)
                    quadrants[idx].z[k] = quadrants[idx_m].z[k] + molec_parameter->L_z;
            }

        // idx_y = 0 --> mirror original quadrant: idx_y_m = cellList_parameter.N_y
        idx_y = 0;
        idx_y_m = cellList_parameter.N_y;
        for(idx_z = 1; idx_z <= cellList_parameter.N_z; ++idx_z)
            for(idx_x = 1; idx_x <= cellList_parameter.N_x; ++idx_x)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
                // linear index of mirror cell
                idx_x_m = idx_x;
                idx_z_m = idx_z;
                const int idx_m = idx_x_m + N_x_ghost * (idx_y_m + N_y_ghost * idx_z_m);

                quadrants[idx].N = quadrants[idx_m].N;
                                quadrants[idx].N_pad = quadrants[idx_m].N_pad;

                                // associate the index of the ghost quadrant to the real index of the mirror
                quadrants[idx].idx = quadrants[idx_m].idx;


                quadrants[idx].f_x = quadrants[idx_m].f_x;
                quadrants[idx].f_y = quadrants[idx_m].f_y;
                quadrants[idx].f_z = quadrants[idx_m].f_z;

                quadrants[idx].x = quadrants[idx_m].x;
                quadrants[idx].z = quadrants[idx_m].z;

                // copy-shift the y coordinats of the mirror cell by -L_y
                MOLEC_MALLOC(quadrants[idx].y, quadrants[idx_m].N_pad * sizeof(float));
                for(int k = 0; k < quadrants[idx_m].N_pad; ++k)
                    quadrants[idx].y[k] = quadrants[idx_m].y[k] - molec_parameter->L_y;
            }

        // idx_y = cellList_parameter.N_y + 1 --> mirror original quadrant: idx_y_m = 1
        idx_y = cellList_parameter.N_y + 1;
        idx_y_m = 1;
        for(idx_z = 1; idx_z <= cellList_parameter.N_z; ++idx_z)
            for(idx_x = 1; idx_x <= cellList_parameter.N_x; ++idx_x)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
                // linear index of mirror cell
                idx_x_m = idx_x;
                idx_z_m = idx_z;
                const int idx_m = idx_x_m + N_x_ghost * (idx_y_m + N_y_ghost * idx_z_m);

                quadrants[idx].N = quadrants[idx_m].N;
                                quadrants[idx].N_pad = quadrants[idx_m].N_pad;

                // associate the index of the ghost quadrant to the real index of the mirror
                quadrants[idx].idx = quadrants[idx_m].idx;

                quadrants[idx].f_x = quadrants[idx_m].f_x;
                quadrants[idx].f_y = quadrants[idx_m].f_y;
                quadrants[idx].f_z = quadrants[idx_m].f_z;

                quadrants[idx].x = quadrants[idx_m].x;
                quadrants[idx].z = quadrants[idx_m].z;

                // copy-shift the y coordinats of the mirror cell by +L_y
                MOLEC_MALLOC(quadrants[idx].y, quadrants[idx_m].N_pad * sizeof(float));
                for(int k = 0; k < quadrants[idx_m].N_pad; ++k)
                    quadrants[idx].y[k] = quadrants[idx_m].y[k] + molec_parameter->L_y;
            }

        // idx_x = 0 --> mirror original quadrant: idx_x_m = cellList_parameter.N_x
        idx_x = 0;
        idx_x_m = cellList_parameter.N_x;
        for(idx_z = 1; idx_z <= cellList_parameter.N_z; ++idx_z)
            for(idx_y = 1; idx_y <= cellList_parameter.N_y; ++idx_y)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
                // linear index of mirror cell
                idx_y_m = idx_y;
                idx_z_m = idx_z;
                const int idx_m = idx_x_m + N_x_ghost * (idx_y_m + N_y_ghost * idx_z_m);

                quadrants[idx].N = quadrants[idx_m].N;
                                quadrants[idx].N_pad = quadrants[idx_m].N_pad;

                // associate the index of the ghost quadrant to the real index of the mirror
                quadrants[idx].idx = quadrants[idx_m].idx;

                quadrants[idx].f_x = quadrants[idx_m].f_x;
                quadrants[idx].f_y = quadrants[idx_m].f_y;
                quadrants[idx].f_z = quadrants[idx_m].f_z;

                quadrants[idx].y = quadrants[idx_m].y;
                quadrants[idx].z = quadrants[idx_m].z;

                // copy-shift the y coordinats of the mirror cell by -L_x
                MOLEC_MALLOC(quadrants[idx].x, quadrants[idx_m].N_pad * sizeof(float));
                for(int k = 0; k < quadrants[idx_m].N_pad; ++k)
                    quadrants[idx].x[k] = quadrants[idx_m].x[k] - molec_parameter->L_x;
            }

        // idx_x = cellList_parameter.N_x + 1 --> mirror original quadrant: idx_x_m = 1
        idx_x = cellList_parameter.N_x + 1;
        idx_x_m = 1;
        for(idx_z = 1; idx_z <= cellList_parameter.N_z; ++idx_z)
            for(idx_y = 1; idx_y <= cellList_parameter.N_y; ++idx_y)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
                // linear index of mirror cell
                idx_y_m = idx_y;
                idx_z_m = idx_z;
                const int idx_m = idx_x_m + N_x_ghost * (idx_y_m + N_y_ghost * idx_z_m);

                quadrants[idx].N = quadrants[idx_m].N;
                                quadrants[idx].N_pad = quadrants[idx_m].N_pad;

                // associate the index of the ghost quadrant to the real index of the mirror
                quadrants[idx].idx = quadrants[idx_m].idx;

                quadrants[idx].f_x = quadrants[idx_m].f_x;
                quadrants[idx].f_y = quadrants[idx_m].f_y;
                quadrants[idx].f_z = quadrants[idx_m].f_z;

                quadrants[idx].y = quadrants[idx_m].y;
                quadrants[idx].z = quadrants[idx_m].z;

                // copy-shift the y coordinats of the mirror cell by +L_x
                MOLEC_MALLOC(quadrants[idx].x, quadrants[idx_m].N_pad * sizeof(float));
                for(int k = 0; k < quadrants[idx_m].N_pad; ++k)
                    quadrants[idx].x[k] = quadrants[idx_m].x[k] + molec_parameter->L_x;
            }
    }

    // edge ghost quadrants
    {
        /*             _________________________
                      / _____________________  /|
                     / / ________4__________/ / |
                    / / /| |               / /  |
                   / / / | |              / / . |
                  /7/ /| | |             /8/ /| |
                 / / / | | |            / / / | |      ^
                / / /  |11 |           / / /|12 |      | Z direction
               / /_/___|_|_|__________/ / / | | |      |
              /_______  _______________/ /  | | |
              | __________3___________ | |  | | |
              | | |    | | |_________| | |__| | |
              | | |    | |___________| | |____| |
              | | |   / / _______2___| | |_  / /
              | | |  / / /           | | |/ / /
              |9| | /5/ /            |10 | /6/
              | | |/ / /             | | |/ /
              | | | / /              | | ' /       ^
              | | |/_/_______________| |  /       / Y dirrection
              | |____________________| | /       /
              |__________1_____________|/


                  --> X direction

        */
        // edges parallel to x direction
        {
            // 1
            // idx_z = 0 && idx_y = 0
            idx_z = 0;
            idx_y = 0;
            idx_z_m = cellList_parameter.N_z;
            idx_y_m = cellList_parameter.N_y;
            for(idx_x = 1; idx_x <= cellList_parameter.N_x; ++idx_x)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
                // linear index of mirror cell
                idx_x_m = idx_x;
                const int idx_m = idx_x_m + N_x_ghost * (idx_y_m + N_y_ghost * idx_z_m);

                quadrants[idx].N = quadrants[idx_m].N;
                                quadrants[idx].N_pad = quadrants[idx_m].N_pad;

                // associate the index of the ghost quadrant to the real index of the mirror
                quadrants[idx].idx = quadrants[idx_m].idx;

                quadrants[idx].f_x = quadrants[idx_m].f_x;
                quadrants[idx].f_y = quadrants[idx_m].f_y;
                quadrants[idx].f_z = quadrants[idx_m].f_z;

                quadrants[idx].x = quadrants[idx_m].x;

                // copy-shift the y and z coordinats of the mirror cell by -L_y and -L_z
                MOLEC_MALLOC(quadrants[idx].y, quadrants[idx_m].N_pad * sizeof(float));
                MOLEC_MALLOC(quadrants[idx].z, quadrants[idx_m].N_pad * sizeof(float));
                for(int k = 0; k < quadrants[idx_m].N_pad; ++k)
                {
                    quadrants[idx].y[k] = quadrants[idx_m].y[k] - molec_parameter->L_y;
                    quadrants[idx].z[k] = quadrants[idx_m].z[k] - molec_parameter->L_z;
                }
            }

            // 2
            // idx_z = 0 && idx_y = cellList_parameter.N_y + 1
            idx_z = 0;
            idx_y = cellList_parameter.N_y + 1;
            idx_z_m = cellList_parameter.N_z;
            idx_y_m = 1;
            for(idx_x = 1; idx_x <= cellList_parameter.N_x; ++idx_x)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
                // linear index of mirror cell
                idx_x_m = idx_x;
                const int idx_m = idx_x_m + N_x_ghost * (idx_y_m + N_y_ghost * idx_z_m);

                quadrants[idx].N = quadrants[idx_m].N;
                                quadrants[idx].N_pad = quadrants[idx_m].N_pad;

                // associate the index of the ghost quadrant to the real index of the mirror
                quadrants[idx].idx = quadrants[idx_m].idx;

                quadrants[idx].f_x = quadrants[idx_m].f_x;
                quadrants[idx].f_y = quadrants[idx_m].f_y;
                quadrants[idx].f_z = quadrants[idx_m].f_z;

                quadrants[idx].x = quadrants[idx_m].x;

                // copy-shift the y and z coordinats of the mirror cell by L_y and -L_z
                MOLEC_MALLOC(quadrants[idx].y, quadrants[idx_m].N_pad * sizeof(float));
                MOLEC_MALLOC(quadrants[idx].z, quadrants[idx_m].N_pad * sizeof(float));
                for(int k = 0; k < quadrants[idx_m].N_pad; ++k)
                {
                    quadrants[idx].y[k] = quadrants[idx_m].y[k] + molec_parameter->L_y;
                    quadrants[idx].z[k] = quadrants[idx_m].z[k] - molec_parameter->L_z;
                }
            }

            // 3
            // idx_z = cellList_parameter.N_z + 1 && idy_y = 0
            idx_z = cellList_parameter.N_z + 1;
            idx_y = 0;
            idx_z_m = 1;
            idx_y_m = cellList_parameter.N_y;
            for(idx_x = 1; idx_x <= cellList_parameter.N_x; ++idx_x)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
                // linear index of mirror cell
                idx_x_m = idx_x;
                const int idx_m = idx_x_m + N_x_ghost * (idx_y_m + N_y_ghost * idx_z_m);

                quadrants[idx].N = quadrants[idx_m].N;
                                quadrants[idx].N_pad = quadrants[idx_m].N_pad;

                // associate the index of the ghost quadrant to the real index of the mirror
                quadrants[idx].idx = quadrants[idx_m].idx;

                quadrants[idx].f_x = quadrants[idx_m].f_x;
                quadrants[idx].f_y = quadrants[idx_m].f_y;
                quadrants[idx].f_z = quadrants[idx_m].f_z;

                quadrants[idx].x = quadrants[idx_m].x;

                // copy-shift the y and z coordinats of the mirror cell by -L_y and +L_z
                MOLEC_MALLOC(quadrants[idx].y, quadrants[idx_m].N_pad * sizeof(float));
                MOLEC_MALLOC(quadrants[idx].z, quadrants[idx_m].N_pad * sizeof(float));
                for(int k = 0; k < quadrants[idx_m].N_pad; ++k)
                {
                    quadrants[idx].y[k] = quadrants[idx_m].y[k] - molec_parameter->L_y;
                    quadrants[idx].z[k] = quadrants[idx_m].z[k] + molec_parameter->L_z;
                }
            }

            // 4
            // idx_z = cellList_parameter.N_z + 1 && idy_y = cellList_parameter.N_y  + 1
            idx_z = cellList_parameter.N_z + 1;
            idx_y = cellList_parameter.N_y + 1;
            idx_z_m = 1;
            idx_y_m = 1;
            for(idx_x = 1; idx_x <= cellList_parameter.N_x; ++idx_x)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
                // linear index of mirror cell
                idx_x_m = idx_x;
                const int idx_m = idx_x_m + N_x_ghost * (idx_y_m + N_y_ghost * idx_z_m);

                quadrants[idx].N = quadrants[idx_m].N;
                                quadrants[idx].N_pad = quadrants[idx_m].N_pad;

                // associate the index of the ghost quadrant to the real index of the mirror
                quadrants[idx].idx = quadrants[idx_m].idx;

                quadrants[idx].f_x = quadrants[idx_m].f_x;
                quadrants[idx].f_y = quadrants[idx_m].f_y;
                quadrants[idx].f_z = quadrants[idx_m].f_z;

                quadrants[idx].x = quadrants[idx_m].x;

                // copy-shift the y and z coordinats of the mirror cell by +L_y and +L_z
                MOLEC_MALLOC(quadrants[idx].y, quadrants[idx_m].N_pad * sizeof(float));
                MOLEC_MALLOC(quadrants[idx].z, quadrants[idx_m].N_pad * sizeof(float));
                for(int k = 0; k < quadrants[idx_m].N_pad; ++k)
                {
                    quadrants[idx].y[k] = quadrants[idx_m].y[k] + molec_parameter->L_y;
                    quadrants[idx].z[k] = quadrants[idx_m].z[k] + molec_parameter->L_z;
                }
            }
        }

        // edges parallel to y direction
        {
            // 5
            // idx_z = 0 && idy_x = 0
            idx_z = 0;
            idx_x = 0;
            idx_z_m = cellList_parameter.N_z;
            idx_x_m = cellList_parameter.N_x;
            for(idx_y = 1; idx_y <= cellList_parameter.N_y; ++idx_y)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
                // linear index of mirror cell
                idx_y_m = idx_y;
                const int idx_m = idx_x_m + N_x_ghost * (idx_y_m + N_y_ghost * idx_z_m);

                quadrants[idx].N = quadrants[idx_m].N;
                                quadrants[idx].N_pad = quadrants[idx_m].N_pad;

                // associate the index of the ghost quadrant to the real index of the mirror
                quadrants[idx].idx = quadrants[idx_m].idx;

                quadrants[idx].f_x = quadrants[idx_m].f_x;
                quadrants[idx].f_y = quadrants[idx_m].f_y;
                quadrants[idx].f_z = quadrants[idx_m].f_z;

                quadrants[idx].y = quadrants[idx_m].y;

                // copy-shift the y and z coordinats of the mirror cell by -L_x and -L_z
                MOLEC_MALLOC(quadrants[idx].x, quadrants[idx_m].N_pad * sizeof(float));
                MOLEC_MALLOC(quadrants[idx].z, quadrants[idx_m].N_pad * sizeof(float));
                for(int k = 0; k < quadrants[idx_m].N_pad; ++k)
                {
                    quadrants[idx].x[k] = quadrants[idx_m].x[k] - molec_parameter->L_x;
                    quadrants[idx].z[k] = quadrants[idx_m].z[k] - molec_parameter->L_z;
                }
            }

            // 6
            // idx_z = 0 && idy_x = cellList_parameter.N_x + 1
            idx_z = 0;
            idx_x = cellList_parameter.N_x + 1;
            idx_z_m = cellList_parameter.N_z;
            idx_x_m = 1;
            for(idx_y = 1; idx_y <= cellList_parameter.N_y; ++idx_y)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
                // linear index of mirror cell
                idx_y_m = idx_y;
                const int idx_m = idx_x_m + N_x_ghost * (idx_y_m + N_y_ghost * idx_z_m);

                quadrants[idx].N = quadrants[idx_m].N;
                                quadrants[idx].N_pad = quadrants[idx_m].N_pad;

                // associate the index of the ghost quadrant to the real index of the mirror
                quadrants[idx].idx = quadrants[idx_m].idx;

                quadrants[idx].f_x = quadrants[idx_m].f_x;
                quadrants[idx].f_y = quadrants[idx_m].f_y;
                quadrants[idx].f_z = quadrants[idx_m].f_z;

                quadrants[idx].y = quadrants[idx_m].y;

                // copy-shift the y and z coordinats of the mirror cell by +L_x and -L_z
                MOLEC_MALLOC(quadrants[idx].x, quadrants[idx_m].N_pad * sizeof(float));
                MOLEC_MALLOC(quadrants[idx].z, quadrants[idx_m].N_pad * sizeof(float));
                for(int k = 0; k < quadrants[idx_m].N_pad; ++k)
                {
                    quadrants[idx].x[k] = quadrants[idx_m].x[k] + molec_parameter->L_x;
                    quadrants[idx].z[k] = quadrants[idx_m].z[k] - molec_parameter->L_z;
                }
            }

            // 7
            // idx_z = cellList_parameter.N_z + 1 && idy_x = 0
            idx_z = cellList_parameter.N_z + 1;
            idx_x = 0;
            idx_z_m = 1;
            idx_x_m = cellList_parameter.N_x;
            for(idx_y = 1; idx_y <= cellList_parameter.N_y; ++idx_y)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
                // linear index of mirror cell
                idx_y_m = idx_y;
                const int idx_m = idx_x_m + N_x_ghost * (idx_y_m + N_y_ghost * idx_z_m);

                quadrants[idx].N = quadrants[idx_m].N;
                                quadrants[idx].N_pad = quadrants[idx_m].N_pad;

                // associate the index of the ghost quadrant to the real index of the mirror
                quadrants[idx].idx = quadrants[idx_m].idx;

                quadrants[idx].f_x = quadrants[idx_m].f_x;
                quadrants[idx].f_y = quadrants[idx_m].f_y;
                quadrants[idx].f_z = quadrants[idx_m].f_z;

                quadrants[idx].y = quadrants[idx_m].y;

                // copy-shift the y and z coordinats of the mirror cell by -L_x and +L_z
                MOLEC_MALLOC(quadrants[idx].x, quadrants[idx_m].N_pad * sizeof(float));
                MOLEC_MALLOC(quadrants[idx].z, quadrants[idx_m].N_pad * sizeof(float));
                for(int k = 0; k < quadrants[idx_m].N_pad; ++k)
                {
                    quadrants[idx].x[k] = quadrants[idx_m].x[k] - molec_parameter->L_x;
                    quadrants[idx].z[k] = quadrants[idx_m].z[k] + molec_parameter->L_z;
                }
            }

            // 8
            // idx_z = cellList_parameter.N_z + 1 && idy_x = cellList_parameter.N_x + 1
            idx_z = cellList_parameter.N_z + 1;
            idx_x = cellList_parameter.N_x + 1;
            idx_z_m = 1;
            idx_x_m = 1;
            for(idx_y = 1; idx_y <= cellList_parameter.N_y; ++idx_y)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
                // linear index of mirror cell
                idx_y_m = idx_y;
                const int idx_m = idx_x_m + N_x_ghost * (idx_y_m + N_y_ghost * idx_z_m);

                quadrants[idx].N = quadrants[idx_m].N;
                                quadrants[idx].N_pad = quadrants[idx_m].N_pad;

                // associate the index of the ghost quadrant to the real index of the mirror
                quadrants[idx].idx = quadrants[idx_m].idx;

                quadrants[idx].f_x = quadrants[idx_m].f_x;
                quadrants[idx].f_y = quadrants[idx_m].f_y;
                quadrants[idx].f_z = quadrants[idx_m].f_z;

                quadrants[idx].y = quadrants[idx_m].y;

                // copy-shift the y and z coordinats of the mirror cell by +L_x and +L_z
                MOLEC_MALLOC(quadrants[idx].x, quadrants[idx_m].N_pad * sizeof(float));
                MOLEC_MALLOC(quadrants[idx].z, quadrants[idx_m].N_pad * sizeof(float));
                for(int k = 0; k < quadrants[idx_m].N_pad; ++k)
                {
                    quadrants[idx].x[k] = quadrants[idx_m].x[k] + molec_parameter->L_x;
                    quadrants[idx].z[k] = quadrants[idx_m].z[k] + molec_parameter->L_z;
                }
            }
        }

        // edges parallel to z direction
        {
            // 9
            // idx_y = 0 && idy_x = 0
            idx_y = 0;
            idx_x = 0;
            idx_y_m = cellList_parameter.N_y;
            idx_x_m = cellList_parameter.N_x;
            for(idx_z = 1; idx_z <= cellList_parameter.N_z; ++idx_z)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
                // linear index of mirror cell
                idx_z_m = idx_z;
                const int idx_m = idx_x_m + N_x_ghost * (idx_y_m + N_y_ghost * idx_z_m);

                quadrants[idx].N = quadrants[idx_m].N;
                                quadrants[idx].N_pad = quadrants[idx_m].N_pad;

                // associate the index of the ghost quadrant to the real index of the mirror
                quadrants[idx].idx = quadrants[idx_m].idx;

                quadrants[idx].f_x = quadrants[idx_m].f_x;
                quadrants[idx].f_y = quadrants[idx_m].f_y;
                quadrants[idx].f_z = quadrants[idx_m].f_z;

                quadrants[idx].z = quadrants[idx_m].z;

                // copy-shift the x and y coordinats of the mirror cell by -L_x and -L_y
                MOLEC_MALLOC(quadrants[idx].x, quadrants[idx_m].N_pad * sizeof(float));
                MOLEC_MALLOC(quadrants[idx].y, quadrants[idx_m].N_pad * sizeof(float));
                for(int k = 0; k < quadrants[idx_m].N_pad; ++k)
                {
                    quadrants[idx].x[k] = quadrants[idx_m].x[k] - molec_parameter->L_x;
                    quadrants[idx].y[k] = quadrants[idx_m].y[k] - molec_parameter->L_y;
                }
            }

            // 10
            // idx_y = 0 && idy_x = cellList_parameter.N_x + 1
            idx_y = 0;
            idx_x = cellList_parameter.N_x + 1;
            idx_y_m = cellList_parameter.N_y;
            idx_x_m = 1;
            for(idx_z = 1; idx_z <= cellList_parameter.N_z; ++idx_z)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
                // linear index of mirror cell
                idx_z_m = idx_z;
                const int idx_m = idx_x_m + N_x_ghost * (idx_y_m + N_y_ghost * idx_z_m);

                quadrants[idx].N = quadrants[idx_m].N;
                                quadrants[idx].N_pad = quadrants[idx_m].N_pad;

                // associate the index of the ghost quadrant to the real index of the mirror
                quadrants[idx].idx = quadrants[idx_m].idx;

                quadrants[idx].f_x = quadrants[idx_m].f_x;
                quadrants[idx].f_y = quadrants[idx_m].f_y;
                quadrants[idx].f_z = quadrants[idx_m].f_z;

                quadrants[idx].z = quadrants[idx_m].z;

                // copy-shift the x and y coordinats of the mirror cell by +L_x and -L_y
                MOLEC_MALLOC(quadrants[idx].x, quadrants[idx_m].N_pad * sizeof(float));
                MOLEC_MALLOC(quadrants[idx].y, quadrants[idx_m].N_pad * sizeof(float));
                for(int k = 0; k < quadrants[idx_m].N_pad; ++k)
                {
                    quadrants[idx].x[k] = quadrants[idx_m].x[k] + molec_parameter->L_x;
                    quadrants[idx].y[k] = quadrants[idx_m].y[k] - molec_parameter->L_y;
                }
            }

            // 11
            // idx_y = cellList_parameter.N_y + 1 && idy_x = 0
            idx_y = cellList_parameter.N_y + 1;
            idx_x = 0;
            idx_y_m = 1;
            idx_x_m = cellList_parameter.N_x;
            for(idx_z = 1; idx_z <= cellList_parameter.N_z; ++idx_z)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
                // linear index of mirror cell
                idx_z_m = idx_z;
                const int idx_m = idx_x_m + N_x_ghost * (idx_y_m + N_y_ghost * idx_z_m);

                quadrants[idx].N = quadrants[idx_m].N;
                                quadrants[idx].N_pad = quadrants[idx_m].N_pad;

                // associate the index of the ghost quadrant to the real index of the mirror
                quadrants[idx].idx = quadrants[idx_m].idx;

                quadrants[idx].f_x = quadrants[idx_m].f_x;
                quadrants[idx].f_y = quadrants[idx_m].f_y;
                quadrants[idx].f_z = quadrants[idx_m].f_z;

                quadrants[idx].z = quadrants[idx_m].z;

                // copy-shift the x and y coordinats of the mirror cell by -L_x and +L_y
                MOLEC_MALLOC(quadrants[idx].x, quadrants[idx_m].N_pad * sizeof(float));
                MOLEC_MALLOC(quadrants[idx].y, quadrants[idx_m].N_pad * sizeof(float));
                for(int k = 0; k < quadrants[idx_m].N_pad; ++k)
                {
                    quadrants[idx].x[k] = quadrants[idx_m].x[k] - molec_parameter->L_x;
                    quadrants[idx].y[k] = quadrants[idx_m].y[k] + molec_parameter->L_y;
                }
            }

            // 12
            // idx_y = cellList_parameter.N_y + 1 && idy_x = cellList_parameter.N_x + 1
            idx_y = cellList_parameter.N_y + 1;
            idx_x = cellList_parameter.N_x + 1;
            idx_y_m = 1;
            idx_x_m = 1;
            for(idx_z = 1; idx_z <= cellList_parameter.N_z; ++idx_z)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
                // linear index of mirror cell
                idx_z_m = idx_z;
                const int idx_m = idx_x_m + N_x_ghost * (idx_y_m + N_y_ghost * idx_z_m);

                quadrants[idx].N = quadrants[idx_m].N;
                                quadrants[idx].N_pad = quadrants[idx_m].N_pad;

                // associate the index of the ghost quadrant to the real index of the mirror
                quadrants[idx].idx = quadrants[idx_m].idx;

                quadrants[idx].f_x = quadrants[idx_m].f_x;
                quadrants[idx].f_y = quadrants[idx_m].f_y;
                quadrants[idx].f_z = quadrants[idx_m].f_z;

                quadrants[idx].z = quadrants[idx_m].z;

                // copy-shift the x and y coordinats of the mirror cell by +L_x and +L_y
                MOLEC_MALLOC(quadrants[idx].x, quadrants[idx_m].N_pad * sizeof(float));
                MOLEC_MALLOC(quadrants[idx].y, quadrants[idx_m].N_pad * sizeof(float));
                for(int k = 0; k < quadrants[idx_m].N_pad; ++k)
                {
                    quadrants[idx].x[k] = quadrants[idx_m].x[k] + molec_parameter->L_x;
                    quadrants[idx].y[k] = quadrants[idx_m].y[k] + molec_parameter->L_y;
                }
            }
        }
    }

    // corners ghost quadrants
    {
        /*             _________________________
                      /7_____________________ 8/|
                     / / ___________________/ / |
                    / / /| |               / /  |
                   / / / | |              / / . |
                  / / /| | |             / / /| |
                 / / / | | |            / / / | |      ^
                / / /  | | |           / / /| | |      | Z direction
               / /_/___|_|_|__________/ / / | | |      |
              /5______________________6/ /  | | |
              | ______________________ | |  | | |
              | | |    | | |_________| | |__| | |
              | | |    |3|___________| | |___4| |
              | | |   / / ___________| | |_  / /
              | | |  / / /           | | |/ / /
              | | | / / /            | | | / /
              | | |/ / /             | | |/ /
              | | | / /              | | ' /       ^
              | | |/_/_______________| |  /       / Y dirrection
              | |____________________| | /       /
              |1______________________2|/


                  --> X direction

        */

        // 1
        // idx_x = 0 && idx_y = 0 && idx_z = 0
        idx_x = 0;
        idx_y = 0;
        idx_z = 0;
        idx_x_m = cellList_parameter.N_x;
        idx_y_m = cellList_parameter.N_y;
        idx_z_m = cellList_parameter.N_z;
        {
            // linear index of ghost cell
            const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
            // linear index of mirror cell
            const int idx_m = idx_x_m + N_x_ghost * (idx_y_m + N_y_ghost * idx_z_m);

            quadrants[idx].N = quadrants[idx_m].N;
                            quadrants[idx].N_pad = quadrants[idx_m].N_pad;

            // associate the index of the ghost quadrant to the real index of the mirror
            quadrants[idx].idx = quadrants[idx_m].idx;

            quadrants[idx].f_x = quadrants[idx_m].f_x;
            quadrants[idx].f_y = quadrants[idx_m].f_y;
            quadrants[idx].f_z = quadrants[idx_m].f_z;

            // copy-shift the x, y and z coordinats of the mirror cell by -L_x, -L_y and -L_z
            MOLEC_MALLOC(quadrants[idx].x, quadrants[idx_m].N_pad * sizeof(float));
            MOLEC_MALLOC(quadrants[idx].y, quadrants[idx_m].N_pad * sizeof(float));
            MOLEC_MALLOC(quadrants[idx].z, quadrants[idx_m].N_pad * sizeof(float));
            for(int k = 0; k < quadrants[idx_m].N_pad; ++k)
            {
                quadrants[idx].x[k] = quadrants[idx_m].x[k] - molec_parameter->L_x;
                quadrants[idx].y[k] = quadrants[idx_m].y[k] - molec_parameter->L_y;
                quadrants[idx].z[k] = quadrants[idx_m].z[k] - molec_parameter->L_z;
            }
        }

        // 2
        // idx_x = cellList_parameter.N_x + 1 && idx_y = 0 && idx_z = 0
        idx_x = cellList_parameter.N_x + 1;
        idx_y = 0;
        idx_z = 0;
        idx_x_m = 1;
        idx_y_m = cellList_parameter.N_y;
        idx_z_m = cellList_parameter.N_z;
        {
            // linear index of ghost cell
            const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
            // linear index of mirror cell
            const int idx_m = idx_x_m + N_x_ghost * (idx_y_m + N_y_ghost * idx_z_m);

            quadrants[idx].N = quadrants[idx_m].N;
                            quadrants[idx].N_pad = quadrants[idx_m].N_pad;

            // associate the index of the ghost quadrant to the real index of the mirror
            quadrants[idx].idx = quadrants[idx_m].idx;

            quadrants[idx].f_x = quadrants[idx_m].f_x;
            quadrants[idx].f_y = quadrants[idx_m].f_y;
            quadrants[idx].f_z = quadrants[idx_m].f_z;

            // copy-shift the x, y and z coordinats of the mirror cell by +L_x, -L_y and -L_z
            MOLEC_MALLOC(quadrants[idx].x, quadrants[idx_m].N_pad * sizeof(float));
            MOLEC_MALLOC(quadrants[idx].y, quadrants[idx_m].N_pad * sizeof(float));
            MOLEC_MALLOC(quadrants[idx].z, quadrants[idx_m].N_pad * sizeof(float));
            for(int k = 0; k < quadrants[idx_m].N; ++k)
            {
                quadrants[idx].x[k] = quadrants[idx_m].x[k] + molec_parameter->L_x;
                quadrants[idx].y[k] = quadrants[idx_m].y[k] - molec_parameter->L_y;
                quadrants[idx].z[k] = quadrants[idx_m].z[k] - molec_parameter->L_z;
            }
        }

        // 3
        // idx_x = 0 && idx_y = cellList_parameter.N_y + 1 && idx_z = 0
        idx_x = 0;
        idx_y = cellList_parameter.N_y + 1;
        idx_z = 0;
        idx_x_m = cellList_parameter.N_x;
        idx_y_m = 1;
        idx_z_m = cellList_parameter.N_z;
        {
            // linear index of ghost cell
            const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
            // linear index of mirror cell
            const int idx_m = idx_x_m + N_x_ghost * (idx_y_m + N_y_ghost * idx_z_m);

            quadrants[idx].N = quadrants[idx_m].N;
                            quadrants[idx].N_pad = quadrants[idx_m].N_pad;

            // associate the index of the ghost quadrant to the real index of the mirror
            quadrants[idx].idx = quadrants[idx_m].idx;

            quadrants[idx].f_x = quadrants[idx_m].f_x;
            quadrants[idx].f_y = quadrants[idx_m].f_y;
            quadrants[idx].f_z = quadrants[idx_m].f_z;

            // copy-shift the x, y and z coordinats of the mirror cell by -L_x, +L_y and -L_z
            MOLEC_MALLOC(quadrants[idx].x, quadrants[idx_m].N_pad * sizeof(float));
            MOLEC_MALLOC(quadrants[idx].y, quadrants[idx_m].N_pad * sizeof(float));
            MOLEC_MALLOC(quadrants[idx].z, quadrants[idx_m].N_pad * sizeof(float));
            for(int k = 0; k < quadrants[idx_m].N; ++k)
            {
                quadrants[idx].x[k] = quadrants[idx_m].x[k] - molec_parameter->L_x;
                quadrants[idx].y[k] = quadrants[idx_m].y[k] + molec_parameter->L_y;
                quadrants[idx].z[k] = quadrants[idx_m].z[k] - molec_parameter->L_z;
            }
        }

        // 4
        // idx_x = cellList_parameter.N_x + 1 && idx_y = cellList_parameter.N_y + 1 && idx_z = 0
        idx_x = cellList_parameter.N_x + 1;
        idx_y = cellList_parameter.N_y + 1;
        idx_z = 0;
        idx_x_m = 1;
        idx_y_m = 1;
        idx_z_m = cellList_parameter.N_z;
        {
            // linear index of ghost cell
            const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
            // linear index of mirror cell
            const int idx_m = idx_x_m + N_x_ghost * (idx_y_m + N_y_ghost * idx_z_m);

            quadrants[idx].N = quadrants[idx_m].N;
                            quadrants[idx].N_pad = quadrants[idx_m].N_pad;

            // associate the index of the ghost quadrant to the real index of the mirror
            quadrants[idx].idx = quadrants[idx_m].idx;

            quadrants[idx].f_x = quadrants[idx_m].f_x;
            quadrants[idx].f_y = quadrants[idx_m].f_y;
            quadrants[idx].f_z = quadrants[idx_m].f_z;

            // copy-shift the x, y and z coordinats of the mirror cell by +L_x, +L_y and -L_z
            MOLEC_MALLOC(quadrants[idx].x, quadrants[idx_m].N_pad * sizeof(float));
            MOLEC_MALLOC(quadrants[idx].y, quadrants[idx_m].N_pad * sizeof(float));
            MOLEC_MALLOC(quadrants[idx].z, quadrants[idx_m].N_pad * sizeof(float));
            for(int k = 0; k < quadrants[idx_m].N_pad; ++k)
            {
                quadrants[idx].x[k] = quadrants[idx_m].x[k] + molec_parameter->L_x;
                quadrants[idx].y[k] = quadrants[idx_m].y[k] + molec_parameter->L_y;
                quadrants[idx].z[k] = quadrants[idx_m].z[k] - molec_parameter->L_z;
            }
        }

        // 5
        // idx_x = 0 && idx_y = 0 && idx_z = cellList_parameter.N_z + 1
        idx_x = 0;
        idx_y = 0;
        idx_z = cellList_parameter.N_z + 1;
        idx_x_m = cellList_parameter.N_x;
        idx_y_m = cellList_parameter.N_y;
        idx_z_m = 1;
        {
            // linear index of ghost cell
            const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
            // linear index of mirror cell
            const int idx_m = idx_x_m + N_x_ghost * (idx_y_m + N_y_ghost * idx_z_m);

            quadrants[idx].N = quadrants[idx_m].N;
                            quadrants[idx].N_pad = quadrants[idx_m].N_pad;

            // associate the index of the ghost quadrant to the real index of the mirror
            quadrants[idx].idx = quadrants[idx_m].idx;

            quadrants[idx].f_x = quadrants[idx_m].f_x;
            quadrants[idx].f_y = quadrants[idx_m].f_y;
            quadrants[idx].f_z = quadrants[idx_m].f_z;

            // copy-shift the x, y and z coordinats of the mirror cell by -L_x, -L_y and +L_z
            MOLEC_MALLOC(quadrants[idx].x, quadrants[idx_m].N_pad * sizeof(float));
            MOLEC_MALLOC(quadrants[idx].y, quadrants[idx_m].N_pad * sizeof(float));
            MOLEC_MALLOC(quadrants[idx].z, quadrants[idx_m].N_pad * sizeof(float));
            for(int k = 0; k < quadrants[idx_m].N_pad; ++k)
            {
                quadrants[idx].x[k] = quadrants[idx_m].x[k] - molec_parameter->L_x;
                quadrants[idx].y[k] = quadrants[idx_m].y[k] - molec_parameter->L_y;
                quadrants[idx].z[k] = quadrants[idx_m].z[k] + molec_parameter->L_z;
            }
        }

        // 6
        // idx_x = cellList_parameter.N_x + 1 && idx_y = 0 && idx_z = cellList_parameter.N_z + 1
        idx_x = cellList_parameter.N_x + 1;
        idx_y = 0;
        idx_z = cellList_parameter.N_z + 1;
        idx_x_m = 1;
        idx_y_m = cellList_parameter.N_y;
        idx_z_m = 1;
        {
            // linear index of ghost cell
            const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
            // linear index of mirror cell
            const int idx_m = idx_x_m + N_x_ghost * (idx_y_m + N_y_ghost * idx_z_m);

            quadrants[idx].N = quadrants[idx_m].N;
                            quadrants[idx].N_pad = quadrants[idx_m].N_pad;

            // associate the index of the ghost quadrant to the real index of the mirror
            quadrants[idx].idx = quadrants[idx_m].idx;

            quadrants[idx].f_x = quadrants[idx_m].f_x;
            quadrants[idx].f_y = quadrants[idx_m].f_y;
            quadrants[idx].f_z = quadrants[idx_m].f_z;

            // copy-shift the x, y and z coordinats of the mirror cell by +L_x, -L_y and +L_z
            MOLEC_MALLOC(quadrants[idx].x, quadrants[idx_m].N_pad * sizeof(float));
            MOLEC_MALLOC(quadrants[idx].y, quadrants[idx_m].N_pad * sizeof(float));
            MOLEC_MALLOC(quadrants[idx].z, quadrants[idx_m].N_pad * sizeof(float));
            for(int k = 0; k < quadrants[idx_m].N_pad; ++k)
            {
                quadrants[idx].x[k] = quadrants[idx_m].x[k] + molec_parameter->L_x;
                quadrants[idx].y[k] = quadrants[idx_m].y[k] - molec_parameter->L_y;
                quadrants[idx].z[k] = quadrants[idx_m].z[k] + molec_parameter->L_z;
            }
        }

        // 7
        // idx_x = 0 && idx_y = cellList_parameter.N_y + 1 && idx_z = cellList_parameter.N_z + 1
        idx_x = 0;
        idx_y = cellList_parameter.N_y + 1;
        idx_z = cellList_parameter.N_z + 1;
        idx_x_m = cellList_parameter.N_x;
        idx_y_m = 1;
        idx_z_m = 1;
        {
            // linear index of ghost cell
            const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
            // linear index of mirror cell
            const int idx_m = idx_x_m + N_x_ghost * (idx_y_m + N_y_ghost * idx_z_m);

            quadrants[idx].N = quadrants[idx_m].N;
                            quadrants[idx].N_pad = quadrants[idx_m].N_pad;

            // associate the index of the ghost quadrant to the real index of the mirror
            quadrants[idx].idx = quadrants[idx_m].idx;

            quadrants[idx].f_x = quadrants[idx_m].f_x;
            quadrants[idx].f_y = quadrants[idx_m].f_y;
            quadrants[idx].f_z = quadrants[idx_m].f_z;

            // copy-shift the x, y and z coordinats of the mirror cell by -L_x, +L_y and +L_z
            MOLEC_MALLOC(quadrants[idx].x, quadrants[idx_m].N_pad * sizeof(float));
            MOLEC_MALLOC(quadrants[idx].y, quadrants[idx_m].N_pad * sizeof(float));
            MOLEC_MALLOC(quadrants[idx].z, quadrants[idx_m].N_pad * sizeof(float));
            for(int k = 0; k < quadrants[idx_m].N_pad; ++k)
            {
                quadrants[idx].x[k] = quadrants[idx_m].x[k] - molec_parameter->L_x;
                quadrants[idx].y[k] = quadrants[idx_m].y[k] + molec_parameter->L_y;
                quadrants[idx].z[k] = quadrants[idx_m].z[k] + molec_parameter->L_z;
            }
        }

        // 8
        // idx_x = cellList_parameter.N_x + 1 && idx_y = cellList_parameter.N_y + 1 && idx_z =
        // cellList_parameter.N_z + 1
        idx_x = cellList_parameter.N_x + 1;
        idx_y = cellList_parameter.N_y + 1;
        idx_z = cellList_parameter.N_z + 1;
        idx_x_m = 1;
        idx_y_m = 1;
        idx_z_m = 1;
        {
            // linear index of ghost cell
            const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
            // linear index of mirror cell
            const int idx_m = idx_x_m + N_x_ghost * (idx_y_m + N_y_ghost * idx_z_m);

            quadrants[idx].N = quadrants[idx_m].N;
                            quadrants[idx].N_pad = quadrants[idx_m].N_pad;

            // associate the index of the ghost quadrant to the real index of the mirror
            quadrants[idx].idx = quadrants[idx_m].idx;

            quadrants[idx].f_x = quadrants[idx_m].f_x;
            quadrants[idx].f_y = quadrants[idx_m].f_y;
            quadrants[idx].f_z = quadrants[idx_m].f_z;

            // copy-shift the x, y and z coordinats of the mirror cell by +L_x, +L_y and +L_z
            MOLEC_MALLOC(quadrants[idx].x, quadrants[idx_m].N_pad * sizeof(float));
            MOLEC_MALLOC(quadrants[idx].y, quadrants[idx_m].N_pad * sizeof(float));
            MOLEC_MALLOC(quadrants[idx].z, quadrants[idx_m].N_pad * sizeof(float));
            for(int k = 0; k < quadrants[idx_m].N_pad; ++k)
            {
                quadrants[idx].x[k] = quadrants[idx_m].x[k] + molec_parameter->L_x;
                quadrants[idx].y[k] = quadrants[idx_m].y[k] + molec_parameter->L_y;
                quadrants[idx].z[k] = quadrants[idx_m].z[k] + molec_parameter->L_z;
            }
        }
    }

    free(cell_idx);
    free(current_particle_number_of_quadrant);

    return quadrants;
}


void molec_quadrants_finalize_ghost(molec_Quadrant_t* quadrants,
                                    molec_CellList_Parameter_t cellList_parameter,
                                    molec_Simulation_SOA_t* sim)
{
    // Number of cells in each direction considering also ghost cells
    const int N_x_ghost = cellList_parameter.N_x + 2;
    const int N_y_ghost = cellList_parameter.N_y + 2;

    // copy data from quadrants to 1D arrays looping over internal cells
    int n_1D = 0;
    for(int idx_z = 1; idx_z <= cellList_parameter.N_z; ++idx_z)
        for(int idx_y = 1; idx_y <= cellList_parameter.N_y; ++idx_y)
            for(int idx_x = 1; idx_x <= cellList_parameter.N_x; ++idx_x)
            {
                // linear index of cell
                int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);

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

    // free memory of internal quadrants
    for(int idx_z = 1; idx_z <= cellList_parameter.N_z; ++idx_z)
        for(int idx_y = 1; idx_y <= cellList_parameter.N_y; ++idx_y)
            for(int idx_x = 1; idx_x <= cellList_parameter.N_x; ++idx_x)
            {
                // linear index of cell
                int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);

                MOLEC_FREE(quadrants[idx].x);
                MOLEC_FREE(quadrants[idx].y);
                MOLEC_FREE(quadrants[idx].z);

                MOLEC_FREE(quadrants[idx].v_x);
                MOLEC_FREE(quadrants[idx].v_y);
                MOLEC_FREE(quadrants[idx].v_z);

                MOLEC_FREE(quadrants[idx].f_x);
                MOLEC_FREE(quadrants[idx].f_y);
                MOLEC_FREE(quadrants[idx].f_z);
            }

    // free memory of ghost quadrants
    int idx_x, idx_y, idx_z;
    // internal 'flat' ghost quadrants
    {
        // idx_z = 0
        idx_z = 0;
        for(idx_y = 1; idx_y <= cellList_parameter.N_y; ++idx_y)
            for(idx_x = 1; idx_x <= cellList_parameter.N_x; ++idx_x)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);

                MOLEC_FREE(quadrants[idx].z);
            }

        // idx_z = cellList_parameter.N_z + 1
        idx_z = cellList_parameter.N_z + 1;
        for(idx_y = 1; idx_y <= cellList_parameter.N_y; ++idx_y)
            for(idx_x = 1; idx_x <= cellList_parameter.N_x; ++idx_x)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);

                MOLEC_FREE(quadrants[idx].z);
            }

        // idx_y = 0
        idx_y = 0;
        for(idx_z = 1; idx_z <= cellList_parameter.N_z; ++idx_z)
            for(idx_x = 1; idx_x <= cellList_parameter.N_x; ++idx_x)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
                MOLEC_FREE(quadrants[idx].y);
            }

        // idx_y = cellList_parameter.N_y + 1
        idx_y = cellList_parameter.N_y + 1;
        for(idx_z = 1; idx_z <= cellList_parameter.N_z; ++idx_z)
            for(idx_x = 1; idx_x <= cellList_parameter.N_x; ++idx_x)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);

                MOLEC_FREE(quadrants[idx].y);
            }

        // idx_x = 0
        idx_x = 0;
        for(idx_z = 1; idx_z <= cellList_parameter.N_z; ++idx_z)
            for(idx_y = 1; idx_y <= cellList_parameter.N_y; ++idx_y)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);

                MOLEC_FREE(quadrants[idx].x);
            }

        // idx_x = cellList_parameter.N_x + 1
        idx_x = cellList_parameter.N_x + 1;
        for(idx_z = 1; idx_z <= cellList_parameter.N_z; ++idx_z)
            for(idx_y = 1; idx_y <= cellList_parameter.N_y; ++idx_y)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);

                MOLEC_FREE(quadrants[idx].x);
            }
    }

    // edge ghost quadrants
    {
        // edges parallel to x direction
        {
            // idx_z = 0 && idy_y = 0
            idx_z = 0;
            idx_y = 0;
            for(idx_x = 1; idx_x <= cellList_parameter.N_x; ++idx_x)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);

                MOLEC_FREE(quadrants[idx].y);
                MOLEC_FREE(quadrants[idx].z);
            }

            // idx_z = 0 && idy_y = cellList_parameter.N_y + 1
            idx_z = 0;
            idx_y = cellList_parameter.N_y + 1;
            for(idx_x = 1; idx_x <= cellList_parameter.N_x; ++idx_x)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);

                MOLEC_FREE(quadrants[idx].y);
                MOLEC_FREE(quadrants[idx].z);
            }

            // idx_z = cellList_parameter.N_z + 1 && idy_y = 0
            idx_z = cellList_parameter.N_z + 1;
            idx_y = 0;
            for(idx_x = 1; idx_x <= cellList_parameter.N_x; ++idx_x)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);

                MOLEC_FREE(quadrants[idx].y);
                MOLEC_FREE(quadrants[idx].z);
            }

            // idx_z = cellList_parameter.N_z + 1 && idy_y = cellList_parameter.N_y  + 1
            idx_z = cellList_parameter.N_z + 1;
            idx_y = cellList_parameter.N_y + 1;
            for(idx_x = 1; idx_x <= cellList_parameter.N_x; ++idx_x)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);

                MOLEC_FREE(quadrants[idx].y);
                MOLEC_FREE(quadrants[idx].z);
            }
        }

        // edges parallel to y direction
        {
            // idx_z = 0 && idy_x = 0
            idx_z = 0;
            idx_x = 0;
            for(idx_y = 1; idx_y <= cellList_parameter.N_y; ++idx_y)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);

                MOLEC_FREE(quadrants[idx].x);
                MOLEC_FREE(quadrants[idx].z);
            }

            // idx_z = 0 && idy_x = cellList_parameter.N_x + 1
            idx_z = 0;
            idx_x = cellList_parameter.N_x + 1;
            for(idx_y = 1; idx_y <= cellList_parameter.N_y; ++idx_y)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);

                MOLEC_FREE(quadrants[idx].x);
                MOLEC_FREE(quadrants[idx].z);
            }

            // idx_z = cellList_parameter.N_z + 1 && idy_x = 0
            idx_z = cellList_parameter.N_z + 1;
            idx_x = 0;
            for(idx_y = 1; idx_y <= cellList_parameter.N_y; ++idx_y)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);

                MOLEC_FREE(quadrants[idx].x);
                MOLEC_FREE(quadrants[idx].z);
            }

            // idx_z = cellList_parameter.N_z + 1 && idy_x = cellList_parameter.N_x + 1
            idx_z = cellList_parameter.N_z + 1;
            idx_x = cellList_parameter.N_x + 1;
            for(idx_y = 1; idx_y <= cellList_parameter.N_y; ++idx_y)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);

                MOLEC_FREE(quadrants[idx].x);
                MOLEC_FREE(quadrants[idx].z);
            }
        }

        // edges parallel to z direction
        {
            // idx_y = 0 && idy_x = 0
            idx_y = 0;
            idx_x = 0;
            for(idx_z = 1; idx_z <= cellList_parameter.N_z; ++idx_z)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);

                MOLEC_FREE(quadrants[idx].x);
                MOLEC_FREE(quadrants[idx].y);
            }

            // idx_y = 0 && idy_x = cellList_parameter.N_x + 1
            idx_y = 0;
            idx_x = cellList_parameter.N_x + 1;
            for(idx_z = 1; idx_z <= cellList_parameter.N_z; ++idx_z)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);

                MOLEC_FREE(quadrants[idx].x);
                MOLEC_FREE(quadrants[idx].y);
            }

            // idx_y = cellList_parameter.N_y + 1 && idy_x = 0
            idx_y = cellList_parameter.N_y + 1;
            idx_x = 0;
            for(idx_z = 1; idx_z <= cellList_parameter.N_z; ++idx_z)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);

                MOLEC_FREE(quadrants[idx].x);
                MOLEC_FREE(quadrants[idx].y);
            }

            // idx_y = cellList_parameter.N_y + 1 && idy_x = cellList_parameter.N_x + 1
            idx_y = cellList_parameter.N_y + 1;
            idx_x = cellList_parameter.N_x + 1;
            for(idx_z = 1; idx_z <= cellList_parameter.N_z; ++idx_z)
            {
                // linear index of ghost cell
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);

                MOLEC_FREE(quadrants[idx].x);
                MOLEC_FREE(quadrants[idx].y);
            }
        }
    }

    // corners ghost quadrants
    {

        // idx_x = 0 && idx_y = 0 && idx_z = 0
        idx_x = 0;
        idx_y = 0;
        idx_z = 0;
        {
            // linear index of ghost cell
            const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);

            MOLEC_FREE(quadrants[idx].x);
            MOLEC_FREE(quadrants[idx].y);
            MOLEC_FREE(quadrants[idx].z);
        }

        // idx_x = cellList_parameter.N_x + 1 && idx_y = 0 && idx_z = 0
        idx_x = cellList_parameter.N_x + 1;
        idx_y = 0;
        idx_z = 0;
        {
            // linear index of ghost cell
            const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);

            MOLEC_FREE(quadrants[idx].x);
            MOLEC_FREE(quadrants[idx].y);
            MOLEC_FREE(quadrants[idx].z);
        }

        // idx_x = 0 && idx_y = cellList_parameter.N_y + 1 && idx_z = 0
        idx_x = 0;
        idx_y = cellList_parameter.N_y + 1;
        idx_z = 0;
        {
            // linear index of ghost cell
            const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);

            MOLEC_FREE(quadrants[idx].x);
            MOLEC_FREE(quadrants[idx].y);
            MOLEC_FREE(quadrants[idx].z);
        }

        // idx_x = cellList_parameter.N_x + 1 && idx_y = cellList_parameter.N_y + 1 && idx_z = 0
        idx_x = cellList_parameter.N_x + 1;
        idx_y = cellList_parameter.N_y + 1;
        idx_z = 0;
        {
            // linear index of ghost cell
            const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);

            MOLEC_FREE(quadrants[idx].x);
            MOLEC_FREE(quadrants[idx].y);
            MOLEC_FREE(quadrants[idx].z);
        }

        // idx_x = 0 && idx_y = 0 && idx_z = cellList_parameter.N_z + 1
        idx_x = 0;
        idx_y = 0;
        idx_z = cellList_parameter.N_z + 1;
        {
            // linear index of ghost cell
            const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);

            MOLEC_FREE(quadrants[idx].x);
            MOLEC_FREE(quadrants[idx].y);
            MOLEC_FREE(quadrants[idx].z);
        }

        // idx_x = cellList_parameter.N_x + 1 && idx_y = 0 && idx_z = cellList_parameter.N_z + 1
        idx_x = cellList_parameter.N_x + 1;
        idx_y = 0;
        idx_z = cellList_parameter.N_z + 1;
        {
            // linear index of ghost cell
            const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);

            MOLEC_FREE(quadrants[idx].x);
            MOLEC_FREE(quadrants[idx].y);
            MOLEC_FREE(quadrants[idx].z);
        }

        // idx_x = 0 && idx_y = cellList_parameter.N_y + 1 && idx_z = cellList_parameter.N_z + 1
        idx_x = 0;
        idx_y = cellList_parameter.N_y + 1;
        idx_z = cellList_parameter.N_z + 1;
        {
            // linear index of ghost cell
            const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);

            MOLEC_FREE(quadrants[idx].x);
            MOLEC_FREE(quadrants[idx].y);
            MOLEC_FREE(quadrants[idx].z);
        }

        // idx_x = cellList_parameter.N_x + 1 && idx_y = cellList_parameter.N_y + 1 && idx_z =
        // cellList_parameter.N_z + 1
        idx_x = cellList_parameter.N_x + 1;
        idx_y = cellList_parameter.N_y + 1;
        idx_z = cellList_parameter.N_z + 1;
        {
            // linear index of ghost cell
            const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);

            MOLEC_FREE(quadrants[idx].x);
            MOLEC_FREE(quadrants[idx].y);
            MOLEC_FREE(quadrants[idx].z);
        }
    }

    free(quadrants);
}
