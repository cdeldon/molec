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
#include <molec/CellVector.h>
#include <string.h>
#include <math.h>

MOLEC_INLINE void get_neighbors(int quadrant_number, int* working_mask)
{
    // TODO
    working_mask[quadrant_number] = 1;
}


void calculate_distance(float* current_quadrant, float* working_set, float* distant_vector)
{
    // TODO
}

void set_update_mask(float* distant_vector, float Rcut2)
{
    // TODO
}
void calculate_forces(float* distant_vector, float* force_vector, int* update_mask)
{
    // TODO
}
void update_forces(float* force_vector, float* f_x, float* f_y, float* f_z, float* Epot)
{
    // TODO
}

void molec_force_quadrant(molec_Simulation_SOA_t* sim, float* Epot, const int N)
{
    molec_CellList_Parameter_t cellList_parameter = molec_parameter->cellList;

    //const float sigLJ = molec_parameter->sigLJ;
    //const float epsLJ = molec_parameter->epsLJ;
    const float Rcut2 = molec_parameter->Rcut2;
    float* x = sim->x;
    float* y = sim->y;
    float* z = sim->z;
    float* distance_vector;
    float* force_vector;
    MOLEC_MALLOC(distance_vector, N * sizeof(float));
    MOLEC_MALLOC(force_vector, 3 * N * sizeof(float));

    // mask of quadrant neighborhood
    int working_mask[cellList_parameter.N];
    // update mask
    int update_mask[N];

    float Epot_ = 0;

    // saves all atoms in working quadrant set
    float* working_set;
    MOLEC_MALLOC(working_set, N * 3 * sizeof(float));
    // saves all atoms in of current quadrant
    float* current_quadrant;
    MOLEC_MALLOC(current_quadrant, N * 3 * sizeof(float))

    int quadrant_number;
    // iterate over all quadrants
    for(quadrant_number = 0; quadrant_number < cellList_parameter.N; ++quadrant_number)
    {
        // construct working set for this quadrant
        int atom_number;
        for(atom_number = 0; atom_number < N; ++atom_number)
        {
            // get corresponding quadrant working set
            get_neighbors(quadrant_number, working_mask); // TODO

            // copy working set
            int k;     // atom number
            int j = 0; // working quadrant index
            int i = 0; // workint current quadrant index
            for(k = 0; k < N; ++k)
            {
                int idx = x[k] / cellList_parameter.N_x;
                int idy = y[k] / cellList_parameter.N_y;
                int idz = z[k] / cellList_parameter.N_z;
                int cellID = idx + cellList_parameter.N_x * (idy + cellList_parameter.N_y * idz);
                if(working_mask[cellID])
                {
                    working_set[j] = x[k];
                    working_set[j + 1] = y[k];
                    working_set[j + 2] = z[k];
                    j += 3;
                }
                if(cellID == quadrant_number)
                {
                    current_quadrant[i] = x[k];
                    current_quadrant[i + 1] = y[k];
                    current_quadrant[i + 2] = z[k];
                    i += 3;
                }
            }
        }

        // calculate the distant matrix
        calculate_distance(current_quadrant, working_set, distance_vector); // TODO

        // get the mask vector for rcut
        set_update_mask(distance_vector, Rcut2); // TODO

        // calculate force for vector
        calculate_forces(distance_vector, force_vector, update_mask); // TODO

        // update forces
        update_forces(force_vector, sim->f_x, sim->f_y, sim->f_z, &Epot_); // TODO
    }


    *Epot = Epot_;
}
