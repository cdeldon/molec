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
 
#ifndef MOLEC_TIMER_H
#define MOLEC_TIMER_H

#include <molec/Common.h>

#ifdef MOLEC_PLATFORM_WINDOWS
#include <intrin.h>
#else /* MOLEC_PLATFORM_POSIX */
#include <unistd.h>
#include <sys/types.h>
#include <sys/syscall.h>
#endif

/**
 * Time stamp counter (TSC)
 */
typedef union
{
    molec_uint64_t int64;
    struct { molec_uint32_t lo, hi; } int32;
} molec_TSC;

#define MOLEC_TSC_VAL(a) ((a).int64)

/**
 * @brief Count the number of cycles since last reset
 * @param cpu_c     molec_TSC to store the current tsc
 */
#ifdef MOLEC_PLATFORM_WINDOWS
#define MOLEC_RDTSC(cpu_c) (cpu_c).int64 = __rdtsc()
#else /* MOLEC_PLATFORM_POSIX */
#define MOLEC_RDTSC(cpu_c) \
    MOLEC_ASM MOLEC_VOLATILE ("rdtsc" : "=a" ((cpu_c).int32.lo), "=d"((cpu_c).int32.hi))
#endif

/**
 * Query cpu-id (this is used to serialize the pipeline)
 */ 
#ifdef MOLEC_PLATFORM_WINDOWS
#define MOLEC_CPUID() { int __cpuInfo__[4]; __cpuid(__cpuInfo__, 0x0); }
#else /* MOLEC_PLATFORM_POSIX */
#define MOLEC_CPUID() \
    MOLEC_ASM MOLEC_VOLATILE ("cpuid" : : "a" (0) : "bx", "cx", "dx" )  
#endif

/**
 * Serialize the pipeline and query the TSC
 * @return 64-bit unsigned integer representing a tick count
 */ 
molec_uint64_t molec_start_tsc();

/**
 * @brief Stop the RDTSC timer and return diffrence from start (in cycles)
 *
 * @param start     64-bit unsigned integer representing a tick count
 * @return number of elapsed cycles since start
 */ 
molec_uint64_t molec_stop_tsc(molec_uint64_t start);

/**
 * Infrastructure to measure runtime using the Time Stamp Counter (TSC)
 *
 * To time an exection:
 * @code{.c}
 *     molec_measurement_init(100); // Allocate space for 100 measurements
 *     
 *     for(int i = 0; i < 100; ++i)
 *     {
 *         molec_measurement_start();
 *
 *         // Do something intresting ...
 *
 *         molec_measurement_stop();
 *     }
 *
 *     printf("Meadian of elapsed cycles: %llu\n, molec_measurement_finish());
 * @endcode
 */
typedef struct molec_Measurement
{
    /** Measured runtimes (in cycles) */
    molec_uint64_t* values;
    
    /** Number of measurements */
    int num_measurements; 
    
    /** Current iteration */
    int iteration;
    
    /** Current tick count returned by molec_start_tsc() */
    molec_uint64_t start;
    
} molec_Measurement_t;

/**
 * Start the measurement by allocating the molec_Measurement_t struct
 *
 * @param num_measurements  Number of measurements    
 */
void molec_measurement_init(const int num_measurements);

/**
 * Start the TSC
 */
void molec_measurement_start();

/**
 * Stop the TSC and register the value
 */
void molec_measurement_stop();

/**
 * Compute the median of all measurements (in cycles)
 *
 * @return meadian of all measurements 
 */
molec_uint64_t molec_measurement_finish();

#endif