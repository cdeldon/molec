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
#include <string.h>
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

void molec_force_gonnet(molec_Simulation_SOA_t* sim, Real* Epot, const int N)
{
    assert(molec_parameter);
    const Real sigLJ = molec_parameter->sigLJ;
    const Real epsLJ = molec_parameter->epsLJ;
    const Real L = molec_parameter->L;
    const Real Rcut2 = molec_parameter->Rcut2;

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


    //======== CELL LIST CONSTRUCTION ========//

    // compute cell index for each particle in the system
    int *c_idx;                     // c_idx[i] contains cell index idx of particle i
    int *n_particles_per_cell;      // n_particles_per_cell[idx] contains the number
                                    // of particles in cell idx
    int max_particles_per_cell = 0;

    MOLEC_MALLOC(c_idx, sizeof(int) * N);
    MOLEC_MALLOC(n_particles_per_cell, sizeof(int) * cellList_parameter.N);
    memset(n_particles_per_cell, 0, sizeof(int) * cellList_parameter.N);

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

        // increase the counter that stores the number of particles located in cell idx
        ++n_particles_per_cell[idx];
    }

    // get the maximum number of particles in the cells (in order to malloc and free
    // n_particles_per_cell only once
    max_particles_per_cell = n_particles_per_cell[0];
    for(int c = 1; c < cellList_parameter.N; ++c)
    {
        max_particles_per_cell = fmax(max_particles_per_cell, n_particles_per_cell[c]);
    }

    // datastructure representing the cell list
    int *head;                      // head[idx] contains the index of first particle in cell idx
                                    // if head[idx] == -1, then cell is empty
    int *lscl;                      // lscl[i] contains index of next particle in same cell as i
                                    // if lscl[i] == -1, then i was last particle in cell idx

    MOLEC_MALLOC(head, sizeof(int) * cellList_parameter.N);
    MOLEC_MALLOC(lscl, sizeof(int) * N);

    // arrays containing the indices od the particles in cell idx and n_idx respectively
    // note that their size is greater (or equal) than the number of particles in the
    // respective cell; this is done to avoid to malloc and free the arrays for every cell iteration
    int *particles_in_cell_idx;
    int *particles_in_cell_n_idx;

    MOLEC_MALLOC(particles_in_cell_idx, max_particles_per_cell * sizeof(int));
    MOLEC_MALLOC(particles_in_cell_n_idx, max_particles_per_cell * sizeof(int));

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


    //======== CELL ITERATION ========//

    // Reset forces
    memset(f_x, 0, N * sizeof(Real));
    memset(f_y, 0, N * sizeof(Real));
    memset(f_z, 0, N * sizeof(Real));

    // Loop over the cells
    for(int idx_z = 0; idx_z < cellList_parameter.N_z; ++idx_z)
        for(int idx_y = 0; idx_y < cellList_parameter.N_y; ++idx_y)
            for(int idx_x = 0; idx_x < cellList_parameter.N_x; ++idx_x)
            {
                // compute scalar cell index
                const int idx
                    = idx_x + cellList_parameter.N_x * (idx_y + cellList_parameter.N_y * idx_z);

                // extract the particles that lie in cell idx
                const int n_particles_in_cell_idx = n_particles_per_cell[idx];
                particles_in_cell_idx[0] = head[idx];

                if(particles_in_cell_idx[0] != -1)
                {
                    for(int p = 1; p < n_particles_in_cell_idx; ++p)
                    {
                        // next particle inside cell idx
                        particles_in_cell_idx[p] = lscl[particles_in_cell_idx[p-1]];
                    }
                }

                // loop over neighbour cells
                for(int d_z = -1; d_z <= 1; ++d_z)
                    for(int d_y = -1; d_y <= 1; ++d_y)
                        for(int d_x = -1; d_x <= 1; ++d_x)
                        {
                            // compute cell index considering periodic BC
                            int n_idx_z = mod(idx_z + d_z, cellList_parameter.N_z);
                            int n_idx_y = mod(idx_y + d_y, cellList_parameter.N_y);
                            int n_idx_x = mod(idx_x + d_x, cellList_parameter.N_x);

                            // linear index of neighbour cell
                            int n_idx = n_idx_x
                                        + cellList_parameter.N_x
                                              * (n_idx_y + cellList_parameter.N_y * n_idx_z);

                            // extract the particles that lie in cell n_idx
                            const int n_particles_in_cell_n_idx = n_particles_per_cell[n_idx];

                            particles_in_cell_n_idx[0] = head[n_idx];
                            if(particles_in_cell_n_idx[0] != -1)
                            {
                                for(int p = 1; p < n_particles_in_cell_n_idx; ++p)
                                {
                                    // next particle inside cell idx
                                    particles_in_cell_n_idx[p] = lscl[particles_in_cell_n_idx[p-1]];
                                }
                            }

                            // iterate over particles in cell idx
                            for(int pi = 0; pi < n_particles_in_cell_idx; ++pi)
                            {
                                int i = particles_in_cell_idx[pi];

                                // local aliases for particle i in cell idx
                                const Real xi = x[i];
                                const Real yi = y[i];
                                const Real zi = z[i];

                                Real f_xi = f_x[i];
                                Real f_yi = f_y[i];
                                Real f_zi = f_z[i];

                                // iterate over particles in cell n_idx
                                for(int pj = 0; pj < n_particles_in_cell_n_idx; ++pj)
                                {
                                    int j = particles_in_cell_n_idx[pj];

                                    // avoid double counting of interactions
                                    if(i < j)
                                    {
                                        // count number of interactions
                                        if(MOLEC_CELLLIST_COUNT_INTERACTION)
                                            ++num_potential_interactions;

                                        const Real xij = dist(xi, x[j], L);
                                        const Real yij = dist(yi, y[j], L);
                                        const Real zij = dist(zi, z[j], L);

                                        const Real r2 = xij * xij + yij * yij + zij * zij;

                                        if(r2 < Rcut2)
                                        {
                                            // count effective number of interactions
                                            if(MOLEC_CELLLIST_COUNT_INTERACTION)
                                                ++num_effective_interactions;

                                            // V(s) = 4 * eps * (s^12 - s^6) with  s = sig/r
                                            const Real s2 = (sigLJ * sigLJ) / r2;
                                            const Real s6 = s2 * s2 * s2;

                                            Epot_ += 4 * epsLJ * (s6 * s6 - s6);

                                            const Real fr = 24 * epsLJ / r2 * (2 * s6 * s6 - s6);

                                            f_xi += fr * xij;
                                            f_yi += fr * yij;
                                            f_zi += fr * zij;

                                            f_x[j] -= fr * xij;
                                            f_y[j] -= fr * yij;
                                            f_z[j] -= fr * zij;
                                        } // end if(r2<R2)

                                    } // end if (i<j)

                                } // finished particles in cell n_idx

                                f_x[i] = f_xi;
                                f_y[i] = f_yi;
                                f_z[i] = f_zi;

                            } // finished particles in cell idx

                        } // end loop over neighbor cells n_idx

            } // end loop over cells idx


    //======== FREE MEMORY ========//

    MOLEC_FREE(c_idx);
    MOLEC_FREE(n_particles_per_cell);

    MOLEC_FREE(head);
    MOLEC_FREE(lscl);

    MOLEC_FREE(particles_in_cell_idx);
    MOLEC_FREE(particles_in_cell_n_idx);

    *Epot = Epot_;

    // print out percentage of effective interactions
    if(MOLEC_CELLLIST_COUNT_INTERACTION)
        printf("\tPercentage of failed potential interactions: %3.2f\n",
               1. - ((double) num_effective_interactions) / ((double) num_potential_interactions));

}
