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

#include <math.h>
#include <molec/Parameter.h>
#include <stdlib.h>

molec_Parameter_t* molec_parameter = NULL;

void molec_parameter_init(int N)
{
    if(molec_parameter)
    {
// We will not deallocate the pointer on Windows as it was allocated in another
// DLL-heap an therefore triggers an exception ...
#ifndef MOLEC_PLATFORM_WINDOWS
        // free(molec_parameter);
        molec_parameter = NULL;
#endif
    }

    molec_parameter = (molec_Parameter_t*) malloc(sizeof(molec_Parameter_t));

    // Set some default parameters
    molec_parameter->N = N;
    molec_parameter->Nstep = 100;
    molec_parameter->dt = 0.005;
    molec_parameter->mass = 1.0;
    molec_parameter->Rcut = 2.5;
    molec_parameter->Rcut2 = 2.5 * 2.5;
    molec_parameter->scaling = 0.05;
    molec_parameter->rho = 1.25;
    molec_parameter->epsLJ = 1.0;
    molec_parameter->sigLJ = 1.0;


    // Initialize the cell list associated with the defined bounding box,
    // and cut off radius
    molec_cell_init();
}

static int parameter_cell_compare(const void* a, const void* b)
{
    return (*(int*) a - *(int*) b);
}

void molec_cell_init()
{
    const int N = molec_parameter->N;
    const float rho = molec_parameter->rho;

    // compute volume of bounding box
    const float BB_vol = ((float) N) / rho;
    const float Rcut = molec_parameter->Rcut;

    // start proposing cell grid size, and check if resulting volume is big enough
    int not_big_enough = 1;
    // minimal cell-grid size
    int cellSize[3] = {3, 3, 3};

    while(not_big_enough)
    {
        // total number of cells
        int num_cells = cellSize[0] * cellSize[1] * cellSize[2];

        // volume achieved with current cells
        float vol = num_cells * Rcut * Rcut * Rcut;

        if(vol < BB_vol)
            not_big_enough = 1;
        else
            not_big_enough = 0;

        // increase size of the smallest dimension
        if(not_big_enough)
        {
            qsort(cellSize, 3, sizeof(int), parameter_cell_compare);
            cellSize[0] += 1;
        }
    }

    molec_parameter->cellList.N_x = cellSize[0];
    molec_parameter->cellList.N_y = cellSize[1];
    molec_parameter->cellList.N_z = cellSize[2];
    molec_parameter->cellList.N = cellSize[0] * cellSize[1] * cellSize[2];

    molec_parameter->cellList.c_x = Rcut;
    molec_parameter->cellList.c_y = Rcut;
    molec_parameter->cellList.c_z = Rcut;

    molec_parameter->L_x = Rcut * molec_parameter->cellList.N_x;
    molec_parameter->L_y = Rcut * molec_parameter->cellList.N_y;
    molec_parameter->L_z = Rcut * molec_parameter->cellList.N_z;

    // compute resulting density
    molec_parameter->rho = ((float) N) / (molec_parameter->cellList.N * pow(Rcut, 3));
}

void molec_print_parameters()
{
    printf("\n      ========= MOLEC - Simulation paramteters =========\n\n");
    printf("      Number of particles: \t\t%d\n", molec_parameter->N);
    printf("      Time step: \t\t\t%f\n", molec_parameter->dt);
    printf("      Particle density: \t\t%f\n", molec_parameter->rho);
    printf("      Bounding box: \t\t\t%2.1f x %2.1f x %2.1f\n", molec_parameter->L_x,
           molec_parameter->L_y, molec_parameter->L_z);
    printf("      Particle mass: \t\t\t%f\n", molec_parameter->mass);
    printf("      Cutoff Radius: \t\t\t%f\n", molec_parameter->Rcut);
    printf("      Lennard Jones:\n \t\t\tepsilon:\t%f\n\t\t\tsigma:\t\t%f\n",
           molec_parameter->epsLJ, molec_parameter->sigLJ);

    if(molec_parameter->cellList.N != 0)
    {
        printf("\n      Cell List:\n");
        printf("\t\tNumber of cells:\t%d x %d x %d\n", molec_parameter->cellList.N_x,
               molec_parameter->cellList.N_y, molec_parameter->cellList.N_z);
        printf("\t\tNTotal number of cells:\t%d\n", molec_parameter->cellList.N);
        printf("\t\tCell lenght:\t\t%2.1f x %2.1f x %2.1f\n", molec_parameter->cellList.c_x,
               molec_parameter->cellList.c_y, molec_parameter->cellList.c_z);
        printf("\t\t<#particles> per cell:\t%d",
               ((int) molec_parameter->N / molec_parameter->cellList.N));
    }
    printf("\n\n");
}
