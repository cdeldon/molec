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

#include "Unittest.h"
#include <molec/Quadrant.h>

int molec_NAtoms;
float molec_Rho;

TEST_CASE(molec_UnittestGhost)
{
    // Initialize simulation
    molec_NAtoms = 10000;
    molec_Rho = 2.5;
    molec_Simulation_SOA_t* sim = molec_setup_simulation_SOA();
    molec_CellList_Parameter_t cellList_parameter = molec_parameter->cellList;

    // size of internal cells
    const int N_x = cellList_parameter.N_x;
    const int N_y = cellList_parameter.N_y;

    // size of celllist with ghost
    const int N_ghost_x = N_x + 2;
    const int N_ghost_y = N_y + 2;

    molec_Quadrant_t * quadrants = molec_quadrant_init_ghost(molec_NAtoms, cellList_parameter, sim);

    // check if the force pointers of ghost cells point to the correct array
    float *f_x_ghost, *f_y_ghost, *f_z_ghost;
    float *f_x_mirror, *f_y_mirror, *f_z_mirror;

    // choose one ghost cell
    int idx_x = 0, idx_y = 3, idx_z = 3;
    int idx = idx_x + N_ghost_x*(idx_y + N_ghost_y * idx_z);

    // get mirror quadrant
    int idx_x_m = N_ghost_x - 2;
    int idx_y_m = idx_y;
    int idx_z_m = idx_z;
    int idx_m = idx_x_m + N_ghost_x*(idx_y_m + N_ghost_y * idx_z_m);

    f_x_ghost =  quadrants[idx].f_x;
    f_y_ghost =  quadrants[idx].f_y;
    f_z_ghost =  quadrants[idx].f_z;

    f_x_mirror = quadrants[idx_m].f_x;
    f_y_mirror = quadrants[idx_m].f_y;
    f_z_mirror = quadrants[idx_m].f_z;


    molec_teardown_simulation_SOA(sim);
}

