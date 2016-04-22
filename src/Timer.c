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

void molec_measurement_init(const int num_timers)
{
    if(measurement)
        molec_error("Multiple measurements ongoing!");

    // Allocate
    measurement = (molec_Measurement_t*) malloc(sizeof(molec_Measurement_t));

    measurement->value_list_heads = malloc(sizeof(molec_Measurement_Node_t*) * num_timers);
    measurement->value_list_tails = malloc(sizeof(molec_Measurement_Node_t*) * num_timers);

    measurement->num_timers = num_timers;
    measurement->num_measurements = malloc(sizeof(int) * num_timers);
    measurement->start = malloc(sizeof(molec_uint64_t) * num_timers);

    // Initialize
    for(int i = 0; i < num_timers; ++i)
    {
        measurement->value_list_heads[i] = NULL;
        measurement->value_list_tails[i] = NULL;
        measurement->num_measurements[i] = 0;
        measurement->start[i] = 0;
    }
}

void molec_measurement_start(int timer_index)
{
    molec_TSC start;
    MOLEC_CPUID();
    MOLEC_RDTSC(start);
    measurement->start[timer_index] = MOLEC_TSC_VAL(start);
}

void molec_measurement_stop(int timer_index)
{
    molec_TSC end;
    MOLEC_RDTSC(end);
    MOLEC_CPUID();

    // Construct node
    molec_Measurement_Node_t* node = malloc(sizeof(molec_Measurement_Node_t));
    node->value = MOLEC_TSC_VAL(end) - measurement->start[timer_index];
    node->next = NULL;

    // Set node
    if(measurement->value_list_heads[timer_index] == NULL)
    {
        measurement->value_list_heads[timer_index] = node;
        measurement->value_list_tails[timer_index] = node;
    }
    else
    {
        measurement->value_list_tails[timer_index]->next = node;
        measurement->value_list_tails[timer_index] = node;
    }

    measurement->num_measurements[timer_index]++;
}

int compare_uint64(const void* a, const void* b)
{
    molec_uint64_t* ia = (molec_uint64_t*) a;
    molec_uint64_t* ib = (molec_uint64_t*) b;

    if(*ia == *ib)
    {
        return 0;
    }
    else if(*ia < *ib)
    {
        return -1;
    }
    else
    {
        return 1;
    }
}

molec_uint64_t molec_measurement_get_median(int timer_index)
{
    int len = measurement->num_measurements[timer_index];
    molec_uint64_t* values = malloc(sizeof(molec_uint64_t) * len);

    // Copy
    molec_Measurement_Node_t* node = measurement->value_list_heads[timer_index];
    for(int i = 0; i < len; ++i)
    {
        values[i] = node->value;
        node = node->next;
    }

    qsort(values, len, sizeof(molec_uint64_t), &compare_uint64);

    molec_uint64_t ret = values[len / 2];

    return ret;
}

void molec_measurement_print()
{
    printf("\n      ==================== MOLEC - Timers ====================\n\n");
    const char* prefix[5]
        = {"cycles", "x 10^3 cycles", "x 10^6 cycles", "x 10^9 cycles", "x 10^12 cycles"};
    for(int timer_index = 0; timer_index < measurement->num_timers; ++timer_index)
    {
        // check wheter this timer has been used
        if(measurement->value_list_heads[timer_index] != NULL)
        {
            double cycles = (double) molec_measurement_get_median(timer_index);
            int prefix_idx = 0;
            while(cycles > 1024 && prefix_idx < 5)
            {
                cycles /= 1024;
                ++prefix_idx;
            }

            printf("      Timer %-15s %9.6f %s\n",
                   MOLEC_MEASUREMENT_GET_TIMER(timer_index), cycles, prefix[prefix_idx]);
        }
    }
}

void molec_measurement_finish()
{
    // iterate over the timers, delete all measurements nodes
    for(int timer_index = 0; timer_index < measurement->num_timers; ++timer_index)
    {
        molec_Measurement_Node_t* current = measurement->value_list_heads[timer_index];
        molec_Measurement_Node_t* next = current->next;

        free(current);
        current = next;
        next = current->next;
    }

    free(measurement->value_list_heads);
    free(measurement->value_list_tails);
    free(measurement->num_measurements);
    free(measurement->start);

    free(measurement);
}
