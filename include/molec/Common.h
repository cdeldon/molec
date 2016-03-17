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

#ifndef MOLEC_COMMON_H
#define MOLEC_COMMON_H

#include <molec/Config.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef MOLEC_SINGLE_PRECISION
typedef float Real;
#else
typedef double Real;
#endif

#ifdef MOLEC_PLATFORM_WINDOWS
typedef unsigned __int64 molec_uint64_t;
typedef unsigned __int32 molec_uint32_t;
#else
typedef unsigned long long molec_uint64_t;
typedef unsigned int molec_uint32_t;
#endif

/**
 * @brief Print an error message to stderr and abort the execution
 *
 * @param format    Printf-like format string
 */
void molec_error(const char* format, ...);

/**
 * Print an array of lenght N to stdout
 */
void molec_print_array(const Real* array, const int N);

/**
 * @brief This macro wraps libc's malloc and checks for success
 *
 * @param ptr   Pointer to the memory block allocated by malloc
 * @param size  Size of the memory block, in bytes
 */
#define MOLEC_MALLOC(ptr, size)                                                                     \
    MOLEC_INTERNAL_MALLOC_LIBC((ptr), (size), #ptr, __FILE__, __LINE__)

/**
 * @brief This macro wraps libc's free and NULLs out the pointer
 *
 * @param ptr   Pointer to be freed
 */
#define MOLEC_FREE(ptr) MOLEC_INTERNAL_FREE_LIBC((ptr))

// Internal macros
#define MOLEC_INTERNAL_FREE_LIBC(ptr)                                                              \
    free((ptr));                                                                                   \
    (ptr) = NULL;

#define MOLEC_INTERNAL_MALLOC_LIBC(ptr, size, ptr_name, file, line)                                \
    if(!((ptr) = malloc((size))))                                                                  \
        molec_error("%s:%i: failed to allocate memory for '%s'", (file), (line), (ptr_name));

#endif

