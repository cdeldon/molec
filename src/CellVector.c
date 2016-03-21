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
#include "molec/Cell_norm_vec.h"

MOLEC_INLINE Real norm(Real[3] r)
{
    return sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
}


MOLEC_INLINE Real[3] molec_cell_vector(Real* point, int cell_nr_1, int cell_nr_2)
{
    int Nx = cellList_parameter.N_x;
    int Ny = cellList_parameter.N_y;

    int nr1 = Nx * Ny;

    int idz1 = cell_nr_1 / nr1;
    int idz2 = cell_nr_2 / nr1;

    int idy1 = cell_nr_1 / Nx - Ny / Nx * idz1;
    int idy2 = cell_nr_2 / Nx - Ny / Nx * idz2;

    int idx1 = cell_nr_1 - Nx * (idy1 + Ny * idz1);
    int idx2 = cell_nr_2 - Nx * (idy2 + Ny * idz2);

    return molec_cell_vector(idx1, idy1, idz1, idx2, idy2, idz2);
}

MOLEC_INLINE Real[3] molec_cell_vector(
    int idx_x, int idx_y, int idx_z, int n_idx_x, int n_idx_y, int n_idx_z)
{
    Real[3] norm_vec;

    norm_vec[0] = (n_idx_x - idx_x); // x direction
    norm_vec[1] = (n_idx_y - idx_y); // y direction
    norm_vec[2] = (n_idx_z - idx_z); // z direction

    Real nmm = norm(norm_vec); // normalize

    norm_vec[0] = norm_vec[0] / nmm;
    norm_vec[1] = norm_vec[1] / nmm;
    norm_vec[2] = norm_vec[2] / nmm;

    return norm_vec;
}

