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

/**
 * Test the timer implementation of molec
 */
TEST_CASE(molec_UnittestTimer)
{
    molec_measurement_init(100); // Allocate space for 100 measurements

    for(int i = 0; i < 100; ++i)
    {
        molec_measurement_start();

        molec_measurement_stop();
    }

    printf("Meadian of elapsed cycles: %llu\n", molec_measurement_finish());
}
