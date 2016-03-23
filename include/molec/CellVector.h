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

#ifndef MOLEC_CELL_NORM_H
#define MOLEC_CELL_NORM_H

#include <molec/CellListParam.h>

#define MOLEC_DIST_PARALLEL 1.0
#define MOLEC_DIST_DIAG_1   0.70710678118
#define MOLEC_DIST_DIAG_2   0.57735026919

struct molec_Cellvector
{
    Real Nx[3] = {MOLEC_DIST_PARALLEL, 0.0, 0.0};
    Real Ny[3] = {0.0, MOLEC_DIST_PARALLEL, 0.0};
    Real Nz[3] = {0.0, 0.0, MOLEC_DIST_PARALLEL};

    Real NxNy1[3] = {MOLEC_DIST_DIAG_1, MOLEC_DIST_DIAG_1, 0.0};
    Real NxNy2[3] = {-MOLEC_DIST_DIAG_1, MOLEC_DIST_DIAG_1, 0.0};

    Real NxNz1[3] = {MOLEC_DIST_DIAG_1, 0.0, MOLEC_DIST_DIAG_1};
    Real NxNz2[3] = {-MOLEC_DIST_DIAG_1, 0.0, MOLEC_DIST_DIAG_1};

    Real NyNz1[3] = {0.0, MOLEC_DIST_DIAG_1, MOLEC_DIST_DIAG_1};
    Real NyNz2[3] = {0.0, -MOLEC_DIST_DIAG_1, MOLEC_DIST_DIAG_1};
	
	Real NxNyNz1[3]={MOLEC_DIST_DIAG_2, MOLEC_DIST_DIAG_2, MOLEC_DIST_DIAG_2};
	Real NxNyNz2[3]={-MOLEC_DIST_DIAG_2, MOLEC_DIST_DIAG_2, MOLEC_DIST_DIAG_2};
	Real NxNyNz3[3]={MOLEC_DIST_DIAG_2, -MOLEC_DIST_DIAG_2, MOLEC_DIST_DIAG_2;
	Real NxNyNz4[3]={MOLEC_DIST_DIAG_2, MOLEC_DIST_DIAG_2, MOLEC_DIST_DIAG_2};
};


/**
 * @brief Return the normal vector between 2 cells
 *
 * Computes the normalizing vector of cell number 1 and number 2 which can be arbitrary between 0
 *and
 * @c molec_CellList_Parameter_t.N
 */
Real[3] molec_cell_vector(int cell_nr_1, int cell_nr_2);


/**
 * @brief Return the normal vector between 2 cells
 *
 * Computes the normalizing vector of cells and number with all indices
 */
Real[3] molec_cell_vector(int idx_x, int idx_y, int idx_z, int n_idx_x, int n_idx_y, int n_idx_z);


Real[3] molec_lookup_vec(int idx_x, int idx_y, int idx_z, int n_idx_x, int n_idx_y, int n_idx_z);
#endif

