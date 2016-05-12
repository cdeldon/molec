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

    const float* x = sim->x;
    const float* y = sim->y;
    const float* z = sim->z;
    float* f_x = sim->f_x;
    float* f_y = sim->f_y;
    float* f_z = sim->f_z;

    float Epot_ = 0;

    molec_Quadrant_t* quadrants = malloc(cellList_parameter.N * sizeof(molec_Quadrant_t));

    for (int i = 0; i < cellList_parameter.N; ++i)
        quadrants[i].N = 0;


    // for each particle compute cell index and count number of particles inside each cell
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
        quadrants[idx].N += 1;
    }

        for (int i = 0; i < cellList_parameter.N; ++i)
        {
            MOLEC_MALLOC(quadrants[i].x, quadrants[i].N * sizeof(float));
            MOLEC_MALLOC(quadrants[i].y, quadrants[i].N * sizeof(float));
            MOLEC_MALLOC(quadrants[i].z, quadrants[i].N * sizeof(float));

            MOLEC_MALLOC(quadrants[i].f_x, quadrants[i].N * sizeof(float));
            MOLEC_MALLOC(quadrants[i].f_y, quadrants[i].N * sizeof(float));
            MOLEC_MALLOC(quadrants[i].f_z, quadrants[i].N * sizeof(float));
        }



        for (int i = 0; i < cellList_parameter.N; ++i)
        {
            MOLEC_FREE(quadrants[i].x);
            MOLEC_FREE(quadrants[i].y);
            MOLEC_FREE(quadrants[i].z);
            MOLEC_FREE(quadrants[i].f_x);
            MOLEC_FREE(quadrants[i].f_y);
            MOLEC_FREE(quadrants[i].f_z);
        }

        free(quadrants);
}
