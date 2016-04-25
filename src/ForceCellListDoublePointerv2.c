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

static int neighbors[27];
molec_CellList_Parameter_t cellList_parameter;
static int average_atoms_per_cells;
static int increased_average;
static int** molec_cellList;

/**
 * Calculate the positive modulo between two integers, used for periodic BC
 */
MOLEC_INLINE int mod(int b, int m)
{
    return (b % m + m) % m;
}

/**
  * Transform cellNumber to respective 1D indices
  */
MOLEC_INLINE void molec_1to3trans(int cellNumber, int* idx_x, int* idx_y, int* idx_z)
{
    *idx_z = cellNumber / (cellList_parameter.N_x * cellList_parameter.N_y);
    *idx_y = cellNumber / cellList_parameter.N_x - cellList_parameter.N_y*(*idx_z);
    *idx_x = cellNumber - cellList_parameter.N_x*(*idx_y + cellList_parameter.N_y * *idx_z);
}

/**
  * Transform 1D indices to respective cellNumber
  */
MOLEC_INLINE void molec_3to1trans(int idx_x, int idx_y, int idx_z, int* cellNumber)
{
    *cellNumber = idx_x + cellList_parameter.N_x*(idx_y + cellList_parameter.N_y * idx_z);
}

/**
  * Calculate neighbors
  */
MOLEC_INLINE void molec_get_neighbors(int cellNumber, int* neighbors)
{
    int idx_x, idx_y, idx_z;
    int index = 0;
    molec_1to3trans(cellNumber, &idx_x, &idx_y, &idx_z);
    for(int d_z = -1; d_z <= 1; ++d_z)
        for(int d_y = -1; d_y <= 1; ++d_y)
            for(int d_x = -1; d_x <= 1; ++d_x)
            {
                int n_idx_z = mod(idx_z + d_z, cellList_parameter.N_z);
                int n_idx_y = mod(idx_y + d_y, cellList_parameter.N_y);
                int n_idx_x = mod(idx_x + d_x, cellList_parameter.N_x);
                molec_3to1trans(n_idx_x, n_idx_y, n_idx_z, &neighbors[index]);
                ++index;
            }
}

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

void molec_force_cellList_double_pointer_v2(molec_Simulation_SOA_t* sim, float* Epot, const int N)
{
    assert(molec_parameter);
    // get molec parameter
    cellList_parameter = molec_parameter->cellList;
    // get average
    average_atoms_per_cells = molec_parameter->N / cellList_parameter.N;
    // increase expectation of atoms in cell
    increased_average = 2 * average_atoms_per_cells;

    int index;
    // allocation of space for cellList
    if(molec_cellList == NULL)
    {
        MOLEC_MALLOC(molec_cellList, cellList_parameter.N * sizeof(int*));
        for(index = 0; index < cellList_parameter.N; ++index)
        {
            MOLEC_MALLOC(molec_cellList[index], increased_average * sizeof(int));
        }
    }

    const float sigLJ = molec_parameter->sigLJ;
    const float epsLJ = molec_parameter->epsLJ;
    const float L = molec_parameter->L;
    const float Rcut2 = molec_parameter->Rcut2;

    // Local aliases
    const float* x = sim->x;
    const float* y = sim->y;
    const float* z = sim->z;
    float* f_x = sim->f_x;
    float* f_y = sim->f_y;
    float* f_z = sim->f_z;

    float Epot_ = 0;

    // number of atoms in the cell[cellnumber]
    int* number_of_atoms_in_cell;
    MOLEC_MALLOC(number_of_atoms_in_cell, sizeof(int) * cellList_parameter.N);
    // set number to zero in each cell
    memset(number_of_atoms_in_cell, 0, sizeof(int) * cellList_parameter.N);

    int* cell_index_of_atom; // storing corresponding cellnumber
    MOLEC_MALLOC(cell_index_of_atom, sizeof(int) * molec_parameter->N);


    for(int i = 0; i < N; ++i) // structure the celllist
    {
        // linear one dimensional index of cell associated to i-th particle
        int cellnumber_in_x_direction = x[i] / cellList_parameter.c_x;
        int cellnumber_in_y_direction = y[i] / cellList_parameter.c_y;
        int cellnumber_in_z_direction = z[i] / cellList_parameter.c_z;

        // set forces to zero
        f_x[i] = f_y[i] = f_z[i] = 0.0;

        // linear index of cell
        int cellNmbr;
        molec_3to1trans(cellnumber_in_x_direction, cellnumber_in_y_direction, cellnumber_in_z_direction, &cellNmbr);
        cell_index_of_atom[i] = cellNmbr;

        molec_cellList[cellNmbr][number_of_atoms_in_cell[cellNmbr]] = i; // segfault?
        number_of_atoms_in_cell[cellNmbr] += 1;
        assert(number_of_atoms_in_cell[cellNmbr] <= increased_average);
    }

    int cell_number_one;
    // loop over all cells
    for(cell_number_one = 0; cell_number_one < cellList_parameter.N; ++cell_number_one)
    {

        // get neighbors
        molec_get_neighbors(cell_number_one, &neighbors[0]);

        // get my index list
        int* index_array_of_atoms_in_cell_one = molec_cellList[cell_number_one];

        int neighbor_number;
        // loop over all neighbors
        for(neighbor_number = 0; neighbor_number < 27; ++neighbor_number)
        {
            int cell_number_two = neighbors[neighbor_number];

            // get corresponding list
            int* index_array_of_atoms_in_cell_two = molec_cellList[cell_number_two];

            // case one for different cells
            if(cell_number_one < cell_number_two)
            {
                // calculate interaction with loop over the two lists
                int k, i;
                for(k = 0; k < number_of_atoms_in_cell[cell_number_one]; ++k) // first list
                {

                    // index one
                    int index_of_atom_one = index_array_of_atoms_in_cell_one[k];

                    const float xi = x[index_of_atom_one];
                    const float yi = y[index_of_atom_one];
                    const float zi = z[index_of_atom_one];

                    float f_xi = f_x[index_of_atom_one];
                    float f_yi = f_y[index_of_atom_one];
                    float f_zi = f_z[index_of_atom_one];

                    for(i = 0; i < number_of_atoms_in_cell[cell_number_two]; ++i) // second list
                    {

                        // index two
                        int index_of_atom_two = index_array_of_atoms_in_cell_two[i];

                        const float xij = dist(xi, x[index_of_atom_two], L);
                        const float yij = dist(yi, y[index_of_atom_two], L);
                        const float zij = dist(zi, z[index_of_atom_two], L);

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

                            f_x[index_of_atom_two] -= fr * xij;
                            f_y[index_of_atom_two] -= fr * yij;
                            f_z[index_of_atom_two] -= fr * zij;
                        }
                    } // end second list

                    f_x[index_of_atom_one] = f_xi;
                    f_y[index_of_atom_one] = f_yi;
                    f_z[index_of_atom_one] = f_zi;

                } // end first list
            }
            else if(cell_number_one == cell_number_two) // case two for the same cell
            {
                // calculate interaction with loop over the two lists
                int k, i;
                for(k = 0; k < number_of_atoms_in_cell[cell_number_one]; ++k) // first list
                {
                    // index one
                    int index_of_atom_one = index_array_of_atoms_in_cell_one[k];

                    const float xi = x[index_of_atom_one];
                    const float yi = y[index_of_atom_one];
                    const float zi = z[index_of_atom_one];

                    float f_xi = f_x[index_of_atom_one];
                    float f_yi = f_y[index_of_atom_one];
                    float f_zi = f_z[index_of_atom_one];

                    for(i = k + 1; i < number_of_atoms_in_cell[cell_number_two]; ++i) // second list
                    {
                        // index two
                        int index_of_atom_two = index_array_of_atoms_in_cell_two[i];

                        const float xij = dist(xi, x[index_of_atom_two], L);
                        const float yij = dist(yi, y[index_of_atom_two], L);
                        const float zij = dist(zi, z[index_of_atom_two], L);

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

                            f_x[index_of_atom_two] -= fr * xij;
                            f_y[index_of_atom_two] -= fr * yij;
                            f_z[index_of_atom_two] -= fr * zij;

                        } // end if

                    } // end second list

                    f_x[index_of_atom_one] = f_xi;
                    f_y[index_of_atom_one] = f_yi;
                    f_z[index_of_atom_one] = f_zi;

                } // end first list

            } // end else

        } // end neighbors

    } // end loop cells

    *Epot = Epot_;
}
