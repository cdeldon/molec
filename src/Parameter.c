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

#include <molec/Parameter.h>
#include <stdlib.h>

molec_Parameter_t* molec_parameter = NULL;

void molec_parameter_init(int N)
{
    molec_parameter = malloc(sizeof(molec_Parameter_t));

    // Set some default parameters
    molec_parameter->N = N;
    molec_parameter->dt = 0.005;
    molec_parameter->mass = 1.0;
    molec_parameter->Rcut = 2.5;
    molec_parameter->scaling = 0.05;
    molec_parameter->L = 10.0;
    molec_parameter->epsLJ = 1.0;
    molec_parameter->sigLJ = 1.0;
}

