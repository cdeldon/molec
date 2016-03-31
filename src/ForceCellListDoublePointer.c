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
#include <molec/Sort.h>
#include <math.h>

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

/**
 * Calculate the positive modulo between two integers, used for periodic BC
 */
MOLEC_INLINE int mod(int b, int m)
{
    return (b % m + m) % m;
}

void molec_force_cellList_dummy(molec_Simulation_SOA_t* sim, Real* Epot, const int N)
{
    assert(molec_parameter);

    molec_uint64_t num_potential_interactions = 0;
    molec_uint32_t num_effective_interactions = 0;

    molec_CellList_Parameter_t cellList_parameter = molec_parameter->cellList;

    // Local aliases
    const Real* x = sim->x;
    const Real* y = sim->y;
    const Real* z = sim->z;
    Real* f_x = sim->f_x;
    Real* f_y = sim->f_y;
    Real* f_z = sim->f_z;

    Real Epot_ = 0;

    // for each particle compute cell index
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
    }

    // count number of particles in each cell
    int* particles_in_cell_idx = malloc(sizeof(int) * cellList_parameter.N);
    for(int idx = 0; idx < cellList_parameter.N; ++idx)
        particles_in_cell_idx[idx] = 0;

    for(int i = 0; i < N; ++i)
        particles_in_cell_idx[c_idx[i]] += 1;

    // count max number of particles per cell
    int max_particles_per_cell = 0;
    for (int idx = 0; idx < cellList_parameter.N; ++idx)
        max_particles_per_cell = fmax(max_particles_per_cell, particles_in_cell_idx[idx]);

    // generate cell list, cellList[idx][k] is the k-th particle of cell idx
    int **cellList = malloc(sizeof(int*) * cellList_parameter.N);
    for(int idx = 0; idx < cellList_parameter.N; ++idx)
        cellList[idx] = malloc(sizeof(int) * max_particles_per_cell);

    // set of counters for next for loop
    int *cellCounter = malloc(sizeof(int) * cellList_parameter.N);
    for(int idx = 0; idx < cellList_parameter.N; ++idx)
        cellCounter[idx] = 0;

    for(int i = 0; i < N; ++i)
    {
        int idx = c_idx[i];
        cellList[idx][cellCounter[idx]++] = i;
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
                const int idx
                    = idx_x + cellList_parameter.N_x * (idx_y + cellList_parameter.N_y * idx_z);

                // loop over neighbour cells
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

                            // iterate over particles in cell idx starting with k = 0;
                            for(int k_idx = 0; k_idx < particles_in_cell_idx[idx]; ++k_idx)
                            {
                                int i = cellList[idx][k_idx];

                                // scan particles in cell n_idx
                                for(int k_n_idx = 0; k_n_idx < particles_in_cell_idx[n_idx]; ++k_n_idx)
                                {
                                    int j = cellList[n_idx][k_n_idx];
                                    // avoid double counting of interactions
                                    if(i < j)
                                    {
                                        // count number of interactions
                                        if(MOLEC_CELLLIST_COUNT_INTERACTION)
                                            ++num_potential_interactions;

                                        const Real xij = dist(x[i], x[j], molec_parameter->L);
                                        const Real yij = dist(y[i], y[j], molec_parameter->L);
                                        const Real zij = dist(z[i], z[j], molec_parameter->L);

                                        const Real r2 = xij * xij + yij * yij + zij * zij;

                                        if(r2 < molec_parameter->Rcut2)
                                        {
                                            // count effective number of interactions
                                            if(MOLEC_CELLLIST_COUNT_INTERACTION)
                                                ++num_effective_interactions;

                                            // V(s) = 4 * eps * (s^12 - s^6) with  s = sig/r
                                            const Real s2 = (molec_parameter->sigLJ * molec_parameter->sigLJ) / r2;
                                            const Real s6 = s2 * s2 * s2;

                                            Epot_ += 4 * molec_parameter->epsLJ * (s6 * s6 - s6);

                                            const Real fr = 24 * molec_parameter->epsLJ / r2 * (2 * s6 * s6 - s6);

                                            sim->f_x[i] += fr * xij;
                                            sim->f_y[i] += fr * yij;
                                            sim->f_z[i] += fr * zij;

                                            sim->f_x[j] -= fr * xij;
                                            sim->f_y[j] -= fr * yij;
                                            sim->f_z[j] -= fr * zij;
                                        }
                                    }
                                }

                            } // finished particles in cell idx

                        } // end loop over neighbur cells n_idx

            } // end loop over cells idx

    // free memory
    free(c_idx);
    free(particles_in_cell_idx);
    for(int idx = 0; idx < cellList_parameter.N; ++idx)
        free(cellList[idx]);
    free(cellList);
    free(cellCounter);

    *Epot = Epot_;

    // print out percentage of effective interactions
    if(MOLEC_CELLLIST_COUNT_INTERACTION)
        printf("\tPercentage of failed potential interactions: %3.2f\n",
               1. - ((double) num_effective_interactions) / ((double) num_potential_interactions));
}

