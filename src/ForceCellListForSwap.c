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

/**
 * @brief Static arrays needed to perform an 'on the fly' sorting
 *
 * The following arrays are used as temporary containers for the
 * newly sorted particles. At the end of each timestep, these temporary
 * arrays are swapped with the ones contained in @c molec_Simulation_SOA_t
 */
static float* x_copy;
static float* y_copy;
static float* z_copy;
static float* v_x_copy;
static float* v_y_copy;
static float* v_z_copy;

void molec_force_cellList_for_swap(molec_Simulation_SOA_t* sim, float* Epot, const int N)
{
    assert(molec_parameter);
    const float sigLJ = molec_parameter->sigLJ;
    const float epsLJ = molec_parameter->epsLJ;
    const float L = molec_parameter->L;
    const float Rcut2 = molec_parameter->Rcut2;

    molec_uint64_t num_potential_interactions = 0;
    molec_uint32_t num_effective_interactions = 0;

    molec_CellList_Parameter_t cellList_parameter = molec_parameter->cellList;

    // If it is the first iteration, we can malloc memory for the static swap arrays
    if(!x_copy)
    {
        MOLEC_MALLOC(x_copy, N * sizeof(float));
        MOLEC_MALLOC(y_copy, N * sizeof(float));
        MOLEC_MALLOC(z_copy, N * sizeof(float));
        MOLEC_MALLOC(v_x_copy, N * sizeof(float));
        MOLEC_MALLOC(v_y_copy, N * sizeof(float));
        MOLEC_MALLOC(v_z_copy, N * sizeof(float));
    }

    // Local aliases
    const float* x = sim->x;
    const float* y = sim->y;
    const float* z = sim->z;
    const float* v_x = sim->v_x;
    const float* v_y = sim->v_y;
    const float* v_z = sim->v_z;
    float* f_x = sim->f_x;
    float* f_y = sim->f_y;
    float* f_z = sim->f_z;

    float Epot_ = 0;

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
        c_idx[i] = idx_x + cellList_parameter.N_x * (idx_y + cellList_parameter.N_y * idx_z);

        // increase the counter that stores the number of particles located in cell idx
        ++n_particles_per_cell[c_idx[i]];
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
                                    // if lscl[i] == -1, then i was last particle of the cell

    MOLEC_MALLOC(head, sizeof(int) * cellList_parameter.N);
    MOLEC_MALLOC(lscl, sizeof(int) * N);

    // arrays containing the indices of the particles in cell idx and n_idx respectively
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
    for(int i = N-1; i >= 0; --i)
    {
        lscl[i] = head[c_idx[i]];
        // now 'i' is the 'head' of cell 'idx'
        head[c_idx[i]] = i;
    }


    //======== CELL ITERATION ========//

    // Reset forces
    memset(f_x, 0, N * sizeof(float));
    memset(f_y, 0, N * sizeof(float));
    memset(f_z, 0, N * sizeof(float));

    // stores the (temporal) accessed position of the particles
    int particle_access_order = 0;
    int *particle_order;
    MOLEC_MALLOC(particle_order, N * sizeof(int));

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

                // store the particle into the swap arrays for the next loop
                for(int p = 0; p < n_particles_in_cell_idx; ++p, ++particle_access_order)
                {
                    int i = particles_in_cell_idx[p];

                    x_copy[particle_access_order] = x[i];
                    y_copy[particle_access_order] = y[i];
                    z_copy[particle_access_order] = z[i];
                    v_x_copy[particle_access_order] = v_x[i];
                    v_y_copy[particle_access_order] = v_y[i];
                    v_z_copy[particle_access_order] = v_z[i];

                    particle_order[particle_access_order] = i;
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
                                const float xi = x[i];
                                const float yi = y[i];
                                const float zi = z[i];

                                float f_xi = f_x[i];
                                float f_yi = f_y[i];
                                float f_zi = f_z[i];

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
                                        } // end if(r2<R2)

                                    } // end if (i<j)

                                } // finished particles in cell n_idx

                                f_x[i] = f_xi;
                                f_y[i] = f_yi;
                                f_z[i] = f_zi;

                            } // finished particles in cell idx

                        } // end loop over neighbor cells n_idx

            } // end loop over cells idx


    // swap the temporary sorted arrays with the default ones
    float *t;

    t = sim->x;
    sim->x = x_copy;
    x_copy = t;

    t = sim->y;
    sim->y = y_copy;
    y_copy = t;

    t = sim->z;
    sim->z = z_copy;
    z_copy = t;

    t = sim->v_x;
    sim->v_x = v_x_copy;
    v_x_copy = t;

    t = sim->v_y;
    sim->v_y = v_y_copy;
    v_y_copy = t;

    t = sim->v_z;
    sim->v_z = v_z_copy;
    v_z_copy = t;

    // reorder the forces such that the new particle ordering is consistent
    // with the calculated forces

    float *f_x_temp, *f_y_temp, *f_z_temp;
    MOLEC_MALLOC(f_x_temp, N * sizeof(float));
    MOLEC_MALLOC(f_y_temp, N * sizeof(float));
    MOLEC_MALLOC(f_z_temp, N * sizeof(float));
    for(int i = 0; i < N; ++i)
    {
        f_x_temp[i] = f_x[particle_order[i]];
        f_y_temp[i] = f_y[particle_order[i]];
        f_z_temp[i] = f_z[particle_order[i]];
    }
    memcpy(f_x, f_x_temp, N * sizeof(float));
    memcpy(f_y, f_y_temp, N * sizeof(float));
    memcpy(f_z, f_z_temp, N * sizeof(float));

    MOLEC_FREE(f_x_temp);
    MOLEC_FREE(f_y_temp);
    MOLEC_FREE(f_z_temp);


    //======== FREE MEMORY ========//

    MOLEC_FREE(c_idx);
    MOLEC_FREE(n_particles_per_cell);

    MOLEC_FREE(head);
    MOLEC_FREE(lscl);

    MOLEC_FREE(particles_in_cell_idx);
    MOLEC_FREE(particles_in_cell_n_idx);

    MOLEC_FREE(particle_order);

    *Epot = Epot_;

    // print out percentage of effective interactions
    if(MOLEC_CELLLIST_COUNT_INTERACTION)
        printf("\tPercentage of failed potential interactions: %3.2f\n",
               1. - ((double) num_effective_interactions) / ((double) num_potential_interactions));

}


