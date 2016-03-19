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

#ifndef MOLEC_CELLLIST_H
#define MOLEC_CELLLIST_H

#include <molec/Common.h>

/**
 * @brief Parameters of the cell list
 *
 * Contains parameters related to the size and topology
 * of the cell list
 *
 */
typedef struct molec_CellList_Parameter
{
    /** Number of cells per dimension */
    int N_x, N_y, N_z;

    /** Size of one cell of the cell list */
    Real c_x, c_y, c_z;

} molec_CellList_Parameter_t;

/**
 * @brief Initializes the parameters of the cell list
 *
 * Computes the size and number of the cells in the cell list
 * from the parameters of the simulation, such as the bounding box
 * extent  @c molec_Parameter.L and cut off radius @c molec_Parameter.Rcut
 */
void molec_cell_init();

#endif

