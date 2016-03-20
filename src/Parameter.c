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

#include <molec/Parameter.h>
#include <stdlib.h>
#include <math.h>

molec_Parameter_t* molec_parameter = NULL;

void molec_parameter_init(int N)
{
    molec_parameter = malloc(sizeof(molec_Parameter_t));

    // Set some default parameters
    molec_parameter->N = N;
    molec_parameter->dt = 0.005;
    molec_parameter->mass = 1.0;
    molec_parameter->Rcut = 2.5;
    molec_parameter->Rcut2 = 2.5 * 2.5;
    molec_parameter->scaling = 0.05;
    molec_parameter->L = 10.0;
    molec_parameter->epsLJ = 1.0;
    molec_parameter->sigLJ = 1.0;

    // Initialize the cell list associated with the defined bounding box,
    // and cut off radius
    molec_cell_init();
}

void molec_cell_init()
{
    const Real L = molec_parameter->L;
    const Real Rcut = molec_parameter->Rcut;

    // compute the number of cells per dimension
    molec_parameter->cellList.N_x = floor(L / Rcut);
    molec_parameter->cellList.N_y = floor(L / Rcut);
    molec_parameter->cellList.N_z = floor(L / Rcut);

    molec_parameter->cellList.N = molec_parameter->cellList.N_x * molec_parameter->cellList.N_y
                                  * molec_parameter->cellList.N_z;

    // compute the size of the cells
    molec_parameter->cellList.c_x = L / molec_parameter->cellList.N_x;
    molec_parameter->cellList.c_y = L / molec_parameter->cellList.N_y;
    molec_parameter->cellList.c_z = L / molec_parameter->cellList.N_z;
}

