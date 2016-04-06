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

#ifndef MOLEC_DUMP_H
#define MOLEC_DUMP_H

#include <molec/Common.h>
#include <molec/Simulation.h>

extern char* dump_file_name;
extern FILE *molec_dump_file;

/**
 * @brief Dumps particle coordinates to file
 *
 * Dumps particle coordinates to file specified in this
 * header file for each timestep.
 * The dumped file has the following form:
 * <N>
 * < x1  y1  z1>
 * < x2  y2  z2>
 *    .......
 * < xN  yN  zN>
 * < x1  y1  z1>
 * < x2  y2  z2>
 *    .......
 * < xN  yN  zN>
 *
 * @see https://en.wikipedia.org/wiki/XYZ_file_format
 */
void molec_dump_coordinates(molec_Simulation_SOA_t* sim, const int N);

#endif
