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

#include <molec/Dump.h>

char *dump_file_name;
FILE *molec_dump_file;

void molec_dump_coordinates(molec_Simulation_SOA_t* sim, const int N)
{
    // print the molecule coordinates in the form
    // 'x1, y1, z1
    //  x2, y2, z2
    //  ...
    //  xN, yN, zN'
    for(int i = 0; i < N; ++i)
        fprintf(molec_dump_file,"%5.6f, %5.6f, %5.6f\n", sim->x[i], sim->y[i], sim->z[i]);
}
