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

/**
 * Test if data is properly aligned by the MOLEC_ALIGNAS macro
 */
TEST_CASE(molec_UnittestAlign)
{

    for (int i = 1; i < 100; i += 1)
    {
        // 16 byte alignment for SSE instructions
        MOLEC_ALIGNAS(16) float x[10 * i];
        CHECK(((unsigned long)x % 16) == 0);

        // 32 byte alignment for AVX instructions
        MOLEC_ALIGNAS(32) float y[10 * i];
        CHECK(((unsigned long)y % 32) == 0);
    }
}
