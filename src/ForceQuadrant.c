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
#include <molec/Timer.h>


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

    float Epot_ = 0;

    // Build the neighbor_cell array only if not initialized before
    if(neighbor_cells == NULL)
    {
        MOLEC_MALLOC(neighbor_cells, cellList_parameter.N * sizeof(int*));
        for(int i = 0; i < cellList_parameter.N; ++i)
            MOLEC_MALLOC(neighbor_cells[i], 27 * sizeof(int));

        molec_build_cell_neighbors(neighbor_cells, cellList_parameter);
    }

    molec_Quadrant_t* quadrants = molec_quadrant_init(N, cellList_parameter, sim);


    // loop over the quadrants
    for(int idx = 0; idx < cellList_parameter.N; ++idx)
    {
        // load the quadrant associated to cell idx
        molec_Quadrant_t q_idx = quadrants[idx];
        int N_idx = q_idx.N;

        for(int neighbor_cell = 0; neighbor_cell < 27; ++neighbor_cell)
        {
            int n_idx = neighbor_cells[idx][neighbor_cell];

            if(idx > n_idx)
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
            else if(idx == n_idx)
            {
                for(int i = 0; i < N_idx; ++i)
                {
                    float xi = q_idx.x[i];
                    float yi = q_idx.y[i];
                    float zi = q_idx.z[i];

                    float f_xi = q_idx.f_x[i];
                    float f_yi = q_idx.f_y[i];
                    float f_zi = q_idx.f_z[i];

                    for(int j = i + 1; j < N_idx; ++j)
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

    molec_quadrants_finalize(quadrants, cellList_parameter, sim);

    *Epot = Epot_;
}
