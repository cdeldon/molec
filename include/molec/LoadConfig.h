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

#ifndef MOLEC_LOADCONFIG_H
#define MOLEC_LOADCONFIG_H

#include <molec/Common.h>
#include <molec/Parameter.h>

#define MOLEC_FILENAME_MAX_LENGTH 128

typedef struct molec_Loader
{
    /** Filename of parameters file */
    char *filename;

} molec_Loader_t;

/**
 * Global access to the loader
 */
extern molec_Loader_t* molec_loader;

/**
 * @brief Loads the simulation parameters into the program from an external file
 *
 * If a valid external file is passed as argument to the executable, the program
 * will run the simulation using the parameters specified in that file
 *
 * @param filename  Path to configuration parameter
 * @param N         Desired number of particles
 * @param rho       Desired particle density
 */
void molec_load_parameters(const char* filename, int verbose, int N, float rho);

#endif

