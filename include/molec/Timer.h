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

#ifndef MOLEC_TIMER_H
#define MOLEC_TIMER_H

#include <molec/Common.h>

#ifdef MOLEC_PLATFORM_WINDOWS
#include <intrin.h>
#else /* MOLEC_PLATFORM_POSIX */
#include <sys/syscall.h>
#include <sys/types.h>
#include <unistd.h>
#endif

/**
 * Time stamp counter (TSC)
 */
typedef union {
    molec_uint64_t int64;
    struct
    {
        molec_uint32_t lo, hi;
    } int32;
} molec_TSC;
#define MOLEC_TSC_VAL(a) ((a).int64)

/**
 *
 */
typedef struct molec_Measurement_Node
{
    molec_uint64_t value;
    struct molec_Measurement_Node* next;
} molec_Measurement_Node_t;

/**
 * @brief Count the number of cycles since last reset
 * @param cpu_c     molec_TSC to store the current tsc
 */
#ifdef MOLEC_PLATFORM_WINDOWS
#define MOLEC_RDTSC(cpu_c) (cpu_c).int64 = __rdtsc()
#else /* MOLEC_PLATFORM_POSIX */
#define MOLEC_RDTSC(cpu_c)                                                                         \
    MOLEC_ASM MOLEC_VOLATILE("rdtsc" : "=a"((cpu_c).int32.lo), "=d"((cpu_c).int32.hi))
#endif

/**
 * Query cpu-id (this is used to serialize the pipeline)
 */
#ifdef MOLEC_PLATFORM_WINDOWS
#define MOLEC_CPUID()                                                                              \
    {                                                                                              \
        int __cpuInfo__[4];                                                                        \
        __cpuid(__cpuInfo__, 0x0);                                                                 \
    }
#else /* MOLEC_PLATFORM_POSIX */
#define MOLEC_CPUID() MOLEC_ASM MOLEC_VOLATILE("cpuid" : : "a"(0) : "bx", "cx", "dx")
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
 *     molec_measurement_init(2); // Number of independent timers
 *
 *     molec_measurement_start(0); // Start timer 0
 *
 *     for(int i = 0; i < 100; ++i)
 *     {
 *         molec_measurement_start(1); // Start timer 1
 *
 *         // Do something intresting ...
 *
 *         molec_measurement_stop(1); // Stop timer 1
 *     }
 *
 *     molec_measurement_start(0); // Stop timer 0
 *
 *
 *     printf("Meadian of elapsed cycles of 0: %llu\n", molec_measurement_get_median(0));
 *
 *     molec_measurement_finish();
 *
 * @endcode
 */
typedef struct molec_Measurement
{
    /** Measured runtimes (in cycles) */
    molec_Measurement_Node_t** value_list_heads;
    molec_Measurement_Node_t** value_list_tails;

    /** Number of concurrent timers */
    int num_timers;

    /** Number of measurements */
    int* num_measurements;

    /** Current tick count returned by molec_start_tsc() */
    molec_uint64_t* start;

} molec_Measurement_t;

/**
 * Start the measurement by allocating the molec_Measurement_t struct
 *
 * @param num_timers        Number of timers
 */
void molec_measurement_init(const int num_timers);

/**
 * Start the TSC
 *
 * @param timer_index   Index of the timer to start
 */
void molec_measurement_start(int timer_index);

/**
 * Stop the TSC and register the value
 *
 * @param timer_index   Index of the timer to stop
 */
void molec_measurement_stop(int timer_index);

/**
 * @brief Compute the median of all measurements (in cycles) for timer Index
 *
 * @param timer_index   Index of the timer
 * @return meadian of all measurement of timer timer_index
 */
molec_uint64_t molec_measurement_get_median(int timer_index);

/**
 * @brief Prints the measured timers in readeable format to the command line
 */
void molec_measurement_print();

/**
 * @brief cleans the timing infrastructure
 */
void molec_measurement_finish();


/**********************************************************************************/
/*                                   R E A D M E                                  */
/**********************************************************************************/
/*
 *  When adding a new timer, check the following:
 *
 *   = MOLEC_MAX_NUM_TIMERS has to be at least as large as the total number of timers
 *
 *   = add the following lines to this file (assuming name of timer to be "POTATO":
 *
 *      "  #ifndef MOLEC_TIME_POTATO
 *          MOLEC_INTERNAL_IGNORE_TIMER(POTATO, <next-number-of-the-sequence>)
 *         #else
 *          MOLEC_INTERNAL_MAKE_TIMER(POTATO, <same-number-as-above>)
 *         #endif
 *      "
 *
 *   = add a line of code in the function "MOLEC_MEASUREMENT_GET_TIMER" corresponding
 *     to the timer you added
 *
 *   = modify the CMakeLists file according to the other examples
 */

#define MOLEC_MAX_NUM_TIMERS 10

#ifdef MOLEC_TIME
#define MOLEC_MEASUREMENT_INIT molec_measurement_init(MOLEC_MAX_NUM_TIMERS)
#define MOLEC_MEASUREMENT_FINISH molec_measurement_finish()
#define MOLEC_MEASUREMENT_PRINT molec_measurement_print()
#define MOLEC_INTERNAL_START_MEASUREMENT(id) molec_measurement_start((id))
#define MOLEC_INTERNAL_STOP_MEASUREMENT(id) molec_measurement_stop((id))
#define MOLEC_INTERNAL_GET_MEDIAN(id) molec_measurement_get_median((id))
#else // MOLEC_TIME
#define MOLEC_MEASUREMENT_INIT (void) 0
#define MOLEC_MEASUREMENT_FINISH (void) 0
#define MOLEC_MEASUREMENT_PRINT (void) 0
#define MOLEC_INTERNAL_START_MEASUREMENT(id) (void) 0
#define MOLEC_INTERNAL_STOP_MEASUREMENT(id) (void) 0
#endif // MOLEC_TIME

#define xstr(s) str(s)
     #define str(s) #s

#define MOLEC_INTERNAL_MAKE_TIMER(name, id)                                                        \
    MOLEC_INLINE void MOLEC_MEASUREMENT_##name##_START() { MOLEC_INTERNAL_START_MEASUREMENT(id); } \
    MOLEC_INLINE void MOLEC_MEASUREMENT_##name##_STOP() { MOLEC_INTERNAL_STOP_MEASUREMENT(id); }   \
    MOLEC_INLINE molec_uint64_t MOLEC_MEASUREMENT_##name##_GET_MEDIAN()                            \
        { return MOLEC_INTERNAL_GET_MEDIAN(id); }                                                  \
    MOLEC_INLINE char* MOLEC_MEASUREMENT_GET_TIMER_##id##_() { return xstr(name); }

#define MOLEC_INTERNAL_IGNORE_TIMER(name, id)                                                      \
    MOLEC_INLINE void MOLEC_MEASUREMENT_##name##_START() { (void) 0; }                             \
    MOLEC_INLINE void MOLEC_MEASUREMENT_##name##_STOP() { (void) 0; }                              \
    MOLEC_INLINE molec_uint64_t MOLEC_MEASUREMENT_##name##_GET_MEDIAN()                            \
        { return 1ul; }                                                                            \
    MOLEC_INLINE char* MOLEC_MEASUREMENT_GET_TIMER_##id##_() { return ""; }


#ifndef MOLEC_TIME_FORCE
MOLEC_INTERNAL_IGNORE_TIMER(FORCE, 0)
#else
MOLEC_INTERNAL_MAKE_TIMER(FORCE, 0)
#endif

#ifndef MOLEC_TIME_INTEGRATOR
MOLEC_INTERNAL_IGNORE_TIMER(INTEGRATOR, 1)
#else
MOLEC_INTERNAL_MAKE_TIMER(INTEGRATOR, 1)
#endif

#ifndef MOLEC_TIME_PERIODIC
MOLEC_INTERNAL_IGNORE_TIMER(PERIODIC, 2)
#else
MOLEC_INTERNAL_MAKE_TIMER(PERIODIC, 2)
#endif

#ifndef MOLEC_TIME_SIMULATION
MOLEC_INTERNAL_IGNORE_TIMER(SIMULATION, 3)
#else
MOLEC_INTERNAL_MAKE_TIMER(SIMULATION, 3)
#endif

MOLEC_INTERNAL_IGNORE_TIMER(CELL_CONSTRUCTION, 4)


MOLEC_INLINE char* MOLEC_MEASUREMENT_GET_TIMER(int id)
{
    switch(id){
        case 0: return MOLEC_MEASUREMENT_GET_TIMER_0_();
        case 1: return MOLEC_MEASUREMENT_GET_TIMER_1_();
        case 2: return MOLEC_MEASUREMENT_GET_TIMER_2_();
        case 3: return MOLEC_MEASUREMENT_GET_TIMER_3_();
        case 4: return MOLEC_MEASUREMENT_GET_TIMER_4_();
        default: molec_error("Index %d does not correspond to any timer\n", id);
    }
    molec_error("Index %d does not correspond to any timer\n", id);
    return "ERROR";
}


#endif
