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

#include <molec/Simulation.h>
#include <molec/Force.h>
#include <molec/LoadConfig.h>
#include <stdlib.h>

int main(int argc, const char* argv[])
{
    srand(42);

    molec_load_parameters(argc, argv, 1);

    // Run the simulation using the routing specified in the passed argument
    // to compute the interactions between particles
    molec_force_calculation force_calculation = &molec_force_cellList;
    molec_run_simulation(force_calculation);
}

