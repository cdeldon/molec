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
    if(molec_parameter)
        MOLEC_FREE(molec_parameter);

    molec_parameter = malloc(sizeof(molec_Parameter_t));

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

    // Determine bounding box size dependind of N and rho
    // such that rho = N/(L*L*L)
    molec_parameter->L = pow(((Real) molec_parameter->N) / molec_parameter->rho,(1./3));

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

void molec_print_parameters()
{
    printf("\n\t========= MOLEC - Simulation paramteters =========\n\n");
    printf("\tNumber of particles: \t\t%d\n", molec_parameter->N);
    printf("\tNumber of steps: \t\t%d\n", molec_parameter->Nstep);
    printf("\tTime step: \t\t\t%f\n", molec_parameter->dt);
    printf("\tParticle density: \t\t%f\n", molec_parameter->rho);
    printf("\tBounding box: \t\t\t%2.1f x %2.1f x %2.1f\n",
           molec_parameter->L, molec_parameter->L, molec_parameter->L);
    printf("\tParticle mass: \t\t\t%f\n", molec_parameter->mass);
    printf("\tCutoff Radius: \t\t\t%f\n", molec_parameter->Rcut);
    printf("\tLennard Jones:\n \t\t\tepsilon:\t%f\n\t\t\tsigma:\t\t%f\n",
           molec_parameter->epsLJ, molec_parameter->sigLJ);

    printf("\n\tCell List:\n");
    printf("\t\tNumber of cells:\t%d x %d x %d\n",
          molec_parameter->cellList.N_x,molec_parameter->cellList.N_y,molec_parameter->cellList.N_z);
    printf("\t\tNTotal number of cells:\t%d\n",molec_parameter->cellList.N);
    printf("\t\tCell lenght:\t\t%2.1f x %2.1f x %2.1f\n",
           molec_parameter->cellList.c_x, molec_parameter->cellList.c_y, molec_parameter->cellList.c_z);
    printf("\t\t<#particles> per cell:\t%d", ((int) molec_parameter->N / molec_parameter->cellList.N));

    printf("\n\n");
}

