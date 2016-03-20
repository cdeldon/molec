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

#include <molec/Common.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void molec_print_array(const Real* array, const int N)
{
    for(int i = 0; i < N; ++i)
        printf("%12.6f\n", array[i]);
}

void molec_error(const char* format, ...)
{
    va_list args;
    va_start(args, format);

    fputs("molec: error: ", stderr);
    vfprintf(stderr, format, args);

    fflush(stderr);
    va_end(args);
    exit(EXIT_FAILURE);
}

