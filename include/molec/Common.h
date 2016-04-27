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

#ifndef MOLEC_COMMON_H
#define MOLEC_COMMON_H


#include <assert.h>
#include <molec/Config.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef MOLEC_PLATFORM_WINDOWS
typedef unsigned __int64 molec_uint64_t;
typedef unsigned __int32 molec_uint32_t;
#else
typedef unsigned long long molec_uint64_t;
typedef unsigned int molec_uint32_t;
#endif

/**
 * @brief if enabled, verbose output
 */
extern int molec_verbose;

/**
 * @brief Print an error message to stderr and abort the execution
 *
 * @param format    Printf-like format string
 */
void molec_error(const char* format, ...);

/**
 * Print an array of lenght N to stdout
 */
void molec_print_array(const float* array, const int N);

/**
  * @brief Print progress bar
  *
  * @param i  steps already done
  * @param n  total number of step
  * @param r  bar resolution
  * @param w  bar width
  */
void molec_progress_bar(int i, int n, int r, int w);

/**
 * @brief This macro wraps libc's malloc and checks for success
 *
 * @param ptr   Pointer to the memory block allocated by malloc
 * @param size  Size of the memory block, in bytes
 */
#define MOLEC_ALIGN_32

#if defined(MOLEC_ALIGN_16)
#define MOLEC_MALLOC(ptr, size) MOLEC_INTERNAL_MALLOC_ALIGN_16((ptr), (size))
#elif defined(MOLEC_ALIGN_32)
#define MOLEC_MALLOC(ptr, size) MOLEC_INTERNAL_MALLOC_ALIGN_32((ptr), (size))
#else
#define MOLEC_MALLOC(ptr, size) MOLEC_INTERNAL_MALLOC_LIBC((ptr), (size), #ptr, __FILE__, __LINE__)
#endif

/**
 * @brief This macro wraps libc's free and NULLs out the pointer
 *
 * @param ptr   Pointer to be freed
 */
#if defined(MOLEC_ALIGN_16) || defined(MOLEC_ALIGN_32)
#define MOLEC_FREE(ptr) MOLEC_INTERNAL_FREE_ALIGN((ptr))
#else
#define MOLEC_FREE(ptr) MOLEC_INTERNAL_FREE_LIBC((ptr))
#endif

// Internal macros
#define MOLEC_INTERNAL_FREE_LIBC(ptr)                                                              \
    {                                                                                              \
        free((ptr));                                                                               \
        (ptr) = NULL;                                                                              \
    }

#define MOLEC_INTERNAL_MALLOC_LIBC(ptr, size, ptr_name, file, line)                                \
    if(!((ptr) = malloc((size))))                                                                  \
        molec_error("%s:%i: failed to allocate memory for '%s'", (file), (line), (ptr_name));

#ifndef MOLEC_PLATFORM_WINDOWS

// Allocate 16 Byte aligned memory
#define MOLEC_INTERNAL_FREE_ALIGN(ptr)                                                             \
    free((ptr));                                                                                   \
    (ptr) = NULL;

// Allocate 16 Byte aligned memory
#define MOLEC_INTERNAL_MALLOC_ALIGN_16(ptr, size) posix_memalign(&(ptr), 16lu, (size));

// Allocate 32 Byte aligned memory
#define MOLEC_INTERNAL_MALLOC_ALIGN_32(ptr, size) posix_memalign(&(ptr), 32lu, (size));

#else

#define MOLEC_INTERNAL_FREE_ALIGN(ptr) { _aligned_free((ptr)); (ptr) = NULL; }

// Allocate 16 Byte aligned memory
#define MOLEC_INTERNAL_MALLOC_ALIGN_16(ptr, size)                                                  \
    {                                                                                              \
        (ptr) = _aligned_malloc((size), 16lu);                                                     \
    }

// Allocate 32 Byte aligned memory
#define MOLEC_INTERNAL_MALLOC_ALIGN_32(ptr, size)                                                  \
    {                                                                                              \
        (ptr) = _aligned_malloc((size), 32lu);                                                     \
    }

#endif

/**
 * Count and report the cell lists interaction miss-rate
 *
 * @see http://tiny.cc/1avaay
 */
#ifndef MOLEC_CELLLIST_COUNT_INTERACTION
#define MOLEC_CELLLIST_COUNT_INTERACTION 0
#endif

/**
 * Dump molecule coordinates to file
 */
#ifndef MOLEC_DUMP_COORDINATES
#define MOLEC_DUMP_COORDINATES 0
#endif

#endif