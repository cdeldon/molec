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

#include <molec/Timer.h>

static molec_Measurement_t* measurement = NULL;

molec_uint64_t molec_start_tsc() 
{
    molec_TSC start;
    MOLEC_CPUID();
    MOLEC_RDTSC(start);
    return MOLEC_TSC_VAL(start);
}

molec_uint64_t molec_stop_tsc(molec_uint64_t start) 
{
    molec_TSC end;
    MOLEC_RDTSC(end);
    MOLEC_CPUID();
    return MOLEC_TSC_VAL(end) - start;
}

void molec_measurement_init(const int num_timers = 1)
{
    if(measurement)
        molec_error("Multiple measurements ongoing!");

    // Allocate
    measurement = (molec_Measurement_t*) malloc(sizeof(molec_Measurement_t));

    measurement->value_list_heads = malloc(sizeof(molec_Measurement_Node_t*) * num_timers);
    measurement->num_timers = num_timers;
    measurement->num_measurements = malloc(sizeof(int) * num_timers);
    measurement->start = malloc(sizeof(molec_uint64_t) * num_timers);

    // Initialize
    for(int i = 0; i < num_timers; ++i)
    {
        measurement->value_list_heads[i] = NULL;
        measurement->num_measurements[i] = 0;
        measurement->start[i] = 0;
    }
}

void molec_measurement_start()
{
    molec_TSC start;
    MOLEC_CPUID();
    MOLEC_RDTSC(start);
    measurement->start = MOLEC_TSC_VAL(start);   
}

void molec_measurement_stop()
{
    molec_TSC end;
    MOLEC_RDTSC(end);
    MOLEC_CPUID();
    measurement->values[measurement->iteration++] = MOLEC_TSC_VAL(end) - measurement->start;
}

int compare_uint64(const void * a, const void * b) 
{
    molec_uint64_t* ia = (molec_uint64_t*) a;
    molec_uint64_t* ib = (molec_uint64_t*) b;

    if ( *ia == *ib ) 
    {
        return 0;
    } 
    else if (*ia < *ib) 
    {
        return -1;
    } 
    else 
    {
        return 1;
    }
}

molec_uint64_t molec_measurement_finish()
{
    qsort(measurement->values, measurement->iteration, sizeof(molec_uint64_t), &compare_uint64);
    molec_uint64_t ret = measurement->values[measurement->iteration / 2];
    
    free(measurement);
    measurement = NULL;
    
    return ret;
}
