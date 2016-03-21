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


Real* molec_cell_vector(int cell_nr_1, int cell_nr_2)
{
    Real* norm_vec MOLEC_MALLOC(norm_vec, 3 * sizeof(Real));
    // idx_x + cellList_parameter.N_x * (idx_y + cellList_parameter.N_y * idx_z);
    int idz1 = cell_nr_1 / (cellList_parameter.N_x * cellList_parameter.N_y);
    int idy1 = cell_nr_1 / cellList_parameter.N_x
               - cellList_parameter.N_y / cellList_parameter.N_x * idz1;
    int idx1 = cell_nr_1 - cellList_parameter.N_x * (idy1 + cellList_parameter.N_y * idz1);

    int idz2 = cell_nr_2 / (cellList_parameter.N_x * cellList_parameter.N_y);
    int idy2 = cell_nr_2 / cellList_parameter.N_x
               - cellList_parameter.N_y / cellList_parameter.N_x * idz2;
    int idx2 = cell_nr_2 - cellList_parameter.N_x * (idy2 + cellList_parameter.N_y * idz2);
    return get_normal_vec(idx1, idy1, idz1, idx2, idy2, idz2);
}

Real* molec_cell_vector(int idx_x, int idx_y, int idx_z, int n_idx_x, int n_idx_y, int n_idx_z)
{
    Real* norm_vec MOLEC_MALLOC(norm_vec, 3 * sizeof(Real));

    Real* md_idx_x = idx_x * c_x + 0.5 * c_x;
    Real* md_idx_y = idx_y * c_y + 0.5 * c_y;
    Real* md_idx_z = idx_z * c_z + 0.5 * c_z;

    Real* md_nidx_x = n_idx_x * c_x + 0.5 * c_x;
    Real* md_nidx_y = n_idx_y * c_y + 0.5 * c_y;
    Real* md_nidx_z = n_idx_z * c_z + 0.5 * c_z;

    norm_vec[0] = abs(md_idx_x - md_nidx_x);
    norm_vec[1] = abs(md_idx_x - md_nidx_x);
    norm_vec[2] = abs(md_idx_x - md_nidx_x);

    Real nmm
        = sqrt(norm_vec[0] * norm_vec[0] + norm_vec[1] * norm_vec[1] + norm_vec[2] * norm_vec[2]);

    norm_vec[0] = norm_vec[0] / nmm;
    norm_vec[1] = norm_vec[1] / nmm;
    norm_vec[2] = norm_vec[2] / nmm;

    return norm_vec;
}

