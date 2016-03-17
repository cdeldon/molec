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

#include <molec/Simulation.h>
#include <stdlib.h>

int main(int argc, const char* argv[])
{
    srand(42);

    int N = argc > 1 ? atoi(argv[1]) : 1000;
    int Nstep = argc > 2 ? atoi(argv[2]) : 10;

    molec_run_simulation_refrence(N, Nstep);
}

