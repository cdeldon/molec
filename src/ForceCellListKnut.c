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

/**
 * Calculate the positive modulo between two integers, used for periodic BC
 */
MOLEC_INLINE int mod(int b, int m)
{
    return (b % m + m) % m;
}

void molec_force_cellList_knut(molec_Simulation_SOA_t* sim, float* Epot, const int N)
{
    assert(molec_parameter);
    const float sigLJ = molec_parameter->sigLJ;
    const float epsLJ = molec_parameter->epsLJ;
    const float L = molec_parameter->L;
    const float Rcut2 = molec_parameter->Rcut2;

    molec_uint64_t num_potential_interactions = 0;
    molec_uint32_t num_effective_interactions = 0;

    molec_CellList_Parameter_t cellList_parameter = molec_parameter->cellList;

    // Local aliases
    const float* x = sim->x;
    const float* y = sim->y;
    const float* z = sim->z;
    float* f_x = sim->f_x;
    float* f_y = sim->f_y;
    float* f_z = sim->f_z;

    float Epot_ = 0;

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

    // cell list construction
    int *head, *lscl;
    MOLEC_MALLOC(head, sizeof(int) * cellList_parameter.N);
    MOLEC_MALLOC(lscl, sizeof(int) * N);

    // fill head with '-1', indicating that the cell is empty
    for(int c = 0; c < cellList_parameter.N; ++c)
        head[c] = -1;

    // generate cell list, lscl[i] contains index of next particle inside
    // the same cell, if lscs[i] == -1, then 'i' was the last particle of the cell
    for(int i = 0; i < N; ++i)
    {
        lscl[i] = head[c_idx[i]];
        head[c_idx[i]] = i;
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

                            // iterate over particles in cell idx starting at the head
                            int i = head[idx];
                            while(i != -1)
                            {
                                // local aliases for particle i in cell idx
                                const float xi = x[i];
                                const float yi = y[i];
                                const float zi = z[i];

                                float f_xi = f_x[i];
                                float f_yi = f_y[i];
                                float f_zi = f_z[i];

                                // scan particles in cell n_idx
                                int j = head[n_idx];
                                while(j != -1)
                                {
                                    // avoid double counting of interactions
                                    if(i < j)
                                    {
                                        // count number of interactions
                                        if(MOLEC_CELLLIST_COUNT_INTERACTION)
                                            ++num_potential_interactions;

                                        const float xij = dist(xi, x[j], L);
                                        const float yij = dist(yi, y[j], L);
                                        const float zij = dist(zi, z[j], L);

                                        const float r2 = xij * xij + yij * yij + zij * zij;

                                        if(r2 < Rcut2)
                                        {
                                            // count effective number of interactions
                                            if(MOLEC_CELLLIST_COUNT_INTERACTION)
                                                ++num_effective_interactions;

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

                                    // next particle inside cell n_idx
                                    j = lscl[j];
                                }

                                f_x[i] = f_xi;
                                f_y[i] = f_yi;
                                f_z[i] = f_zi;

                                // next particle inside cell idx
                                i = lscl[i];

                            } // finished particles in cell idx

                        } // end loop over neighbur cells n_idx

            } // end loop over cells idx

    // free memory
    MOLEC_FREE(c_idx);
    MOLEC_FREE(head);
    MOLEC_FREE(lscl);

    *Epot = Epot_;

    // print out percentage of effective interactions
    if(MOLEC_CELLLIST_COUNT_INTERACTION)
        printf("\tPercentage of failed potential interactions: %3.2f\n",
               1. - ((double) num_effective_interactions) / ((double) num_potential_interactions));
}
