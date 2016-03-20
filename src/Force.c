/*                   _
 *   _ __ ___   ___ | | ___  ___
 *  | '_ ` _ \ / _ \| |/ _ \/ __|
 *  | | | | | | (_) | |  __/ (__
 *  |_| |_| |_|\___/|_|\___|\___| - Molecular Dynamics Framework
 *
 *  Copyright (C) 2016  Carlo Del Don  (deldonc@student.ethz.ch)
 *                      Michel Breyer  (mbreyer@student.ethz.ch)
 *                      Florian Frei   (frofrei@student.ethz.ch)
 *                      Fabian Thuring (thfabian@student.ethz.ch)
 *
 *  This file is distributed under the MIT Open Source License.
 *  See LICENSE.txt for details.
 */

#include <molec/Force.h>
#include <molec/Parameter.h>

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

void molec_force_N2_refrence(molec_Simulation_SOA_t* sim, Real* Epot, const int N)
{
    assert(molec_parameter);
    const Real sigLJ = molec_parameter->sigLJ;
    const Real epsLJ = molec_parameter->epsLJ;
    const Real L = molec_parameter->L;
    const Real Rcut2 = molec_parameter->Rcut2;

    // Local aliases
    const Real* x = sim->x;
    const Real* y = sim->y;
    const Real* z = sim->z;
    Real* f_x = sim->f_x;
    Real* f_y = sim->f_y;
    Real* f_z = sim->f_z;

    Real Epot_ = 0;

    // Reset forces
    for(int i = 0; i < N; ++i)
        f_x[i] = f_y[i] = f_z[i] = 0.0;

    for(int i = 0; i < N; ++i)
    {
        const Real xi = x[i];
        const Real yi = y[i];
        const Real zi = z[i];

        Real f_xi = f_x[i];
        Real f_yi = f_y[i];
        Real f_zi = f_z[i];

        for(int j = i + 1; j < N; ++j)
        {
            const Real xij = dist(xi, x[j], L);
            const Real yij = dist(yi, y[j], L);
            const Real zij = dist(zi, z[j], L);

            const Real r2 = xij * xij + yij * yij + zij * zij;

            if(r2 < Rcut2)
            {
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
            }
        }

        f_x[i] = f_xi;
        f_y[i] = f_yi;
        f_z[i] = f_zi;
    }

    *Epot = Epot_;
}

void molec_force_cellList(molec_Simulation_SOA_t* sim, Real* Epot, const int N)
{
    assert(molec_parameter);
    const Real sigLJ = molec_parameter->sigLJ;
    const Real epsLJ = molec_parameter->epsLJ;
    const Real L = molec_parameter->L;
    const Real Rcut2 = molec_parameter->Rcut2;

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
    int* c_idx;
    MOLEC_MALLOC(c_idx, sizeof(int) * N);

    /* note that the cell list will be traversed in the following form:
     *    for(z=0;z<cellList_parameter.N_z;++z)
     *      for(y=0;y<cellList_parameter.N_y;++y)
     *          for(x=0;x<cellList_parameter.N_x;++x)
     *
     * the fastet running index is x, while the slowest is z
     */
    for(int i = 0; i < N; ++i)
    {
        // linear one dimensional index of cell associated to i-th particle
        int idx_x = x[i]/cellList_parameter.c_x;
        int idx_y = y[i]/cellList_parameter.c_y;
        int idx_z = z[i]/cellList_parameter.c_z;

        // linear index of cell
        int idx = idx_x + cellList_parameter.N_x*(
                      idx_y + cellList_parameter.N_y*idx_z);
        c_idx[i] = idx;
    }

    // cell list construction
    int* head, * lscl;
    MOLEC_MALLOC(head, sizeof(int)*cellList_parameter.N);
    MOLEC_MALLOC(lscl, sizeof(int)*N);

    // fill head with '-1', indicating that the cell is empty
    for(int c = 0; c < cellList_parameter.N; ++c)
        head[c] = -1;

    for(int i = 0; i < N; ++i)
    {
        lscl[i] = head[c_idx[i]];
        head[c_idx[i]] = i;
    }

    // Reset forces
    for(int i = 0; i < N; ++i)
        f_x[i] = f_y[i] = f_z[i] = 0.0;

    // Loop over the cells
    for (int idx_z = 0; idx_z < cellList_parameter.N_z; ++idx_z)
    for (int idx_y = 0; idx_y < cellList_parameter.N_y; ++idx_y)
    for (int idx_x = 0; idx_x < cellList_parameter.N_x; ++idx_x)
    {
        // compute scalar cell index
        const int idx = idx_x + cellList_parameter.N_x*(
                            idx_y + cellList_parameter.N_y*idx_z);

        // loop over neighbour cells
        for (int n_idx_z = idx_z-1; n_idx_z <= idx_z + 1; ++n_idx_z)
        for (int n_idx_y = idx_y-1; n_idx_y <= idx_y + 1; ++n_idx_y)
        for (int n_idx_x = idx_x-1; n_idx_x <= idx_x + 1; ++n_idx_x)
        {
            int shift_x = 0;
            int shift_y = 0;
            int shift_z = 0;
            // Periodic boundary condition by shifting coordinates
            if(n_idx_z < 0)
                shift_z = 1;
            else if(n_idx_z == cellList_parameter.N_z)
                shift_z = -1;

            if(n_idx_y < 0)
                shift_y = 1;
            else if(n_idx_y == cellList_parameter.N_y)
                shift_y = -1;

            if(n_idx_x < 0)
                shift_x = 1;
            else if(n_idx_x == cellList_parameter.N_x)
                shift_x = -1;

            // FIXME -> can compute scalar index of neighbour cell by simple
            // computing with modulos!!
            int n_idx = (n_idx_x + shift_x) + cellList_parameter.N_x*(
                            (n_idx_y + shift_y) + cellList_parameter.N_y*(
                                n_idx_z + shift_z));

            // iterate over particles in cell idx starting at the head
            int i = head[idx];
            while(i != -1)
            {
                // local aliases for particle i in cell idx
                const Real xi = x[i];
                const Real yi = y[i];
                const Real zi = z[i];

                Real f_xi = f_x[i];
                Real f_yi = f_y[i];
                Real f_zi = f_z[i];

                // scan particles in cell n_idx
                int j = head[n_idx];
                while(j != -1)
                {
                    // avoid double counting of interactions
                    if (i < j)
                    {
                        const Real xij = dist(xi, x[j], L);
                        const Real yij = dist(yi, y[j], L);
                        const Real zij = dist(zi, z[j], L);

                        const Real r2 = xij * xij + yij * yij + zij * zij;

                        if(r2 < Rcut2)
                        {
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
}

