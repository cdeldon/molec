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

#include "Unittest.h"
#include <molec/Timer.h>

enum timer_t
{
   TIMER_FORCE = 0,
   TIMER_INTEGRATOR
};


/**
 * Test the timer implementation of molec
 */
TEST_CASE(molec_UnittestTimer)
{
    molec_measurement_init(2);

    for(int i = 0; i < 100; ++i)
    {
        molec_measurement_start(TIMER_FORCE);
        molec_measurement_start(TIMER_INTEGRATOR);
        molec_measurement_stop(TIMER_FORCE);
        molec_measurement_stop(TIMER_INTEGRATOR);
    }

    //printf("Meadian of elapsed cycles: %llu\n", molec_measurement_get_median(0));
    //printf("Meadian of elapsed cycles: %llu\n", molec_measurement_get_median(1));

    molec_measurement_finish();
}
