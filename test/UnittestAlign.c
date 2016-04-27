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
 * Test if data is properly aligned by the MOLEC macros
 */
TEST_CASE(molec_UnittestAlign)
{
    // 16 byte alignment for SSE instructions
    MOLEC_ALIGNAS(16) float array16_00[4];
    MOLEC_ALIGNAS(16) float array16_01[16];
    MOLEC_ALIGNAS(16) float array16_02[32];
    MOLEC_ALIGNAS(16) float array16_03[64];
    CHECK(((unsigned long)array16_00 % 16lu) == 0);
    CHECK(((unsigned long)array16_01 % 16lu) == 0);
    CHECK(((unsigned long)array16_02 % 16lu) == 0);
    CHECK(((unsigned long)array16_03 % 16lu) == 0);

    // 32 byte alignment for AVX instructions
    MOLEC_ALIGNAS(32) float array32_00[4];
    MOLEC_ALIGNAS(32) float array32_01[16];
    MOLEC_ALIGNAS(32) float array32_02[32];
    MOLEC_ALIGNAS(32) float array32_03[64];
    CHECK(((unsigned long)array32_00 % 32lu) == 0);
    CHECK(((unsigned long)array32_01 % 32lu) == 0);
    CHECK(((unsigned long)array32_02 % 32lu) == 0);
    CHECK(((unsigned long)array32_03 % 32lu) == 0);

#if defined(MOLEC_ALIGN_16)
    float* array;

    MOLEC_MALLOC(array, sizeof(float) * 4);
    CHECK(((unsigned long)array % 16lu) == 0);
    MOLEC_FREE(array);

    MOLEC_MALLOC(array, sizeof(float) * 16);
    CHECK(((unsigned long)array % 16lu) == 0);
    MOLEC_FREE(array);

    MOLEC_MALLOC(array, sizeof(float) * 32);
    CHECK(((unsigned long)array % 16lu) == 0);
    MOLEC_FREE(array);

    MOLEC_MALLOC(array, sizeof(float) * 64);
    CHECK(((unsigned long)array % 16lu) == 0);
    MOLEC_FREE(array);

#endif

#if defined(MOLEC_ALIGN_32)
    float* array;

    MOLEC_MALLOC(array, sizeof(float) * 4);
    CHECK(((unsigned long)array % 32lu) == 0);
    MOLEC_FREE(array);

    MOLEC_MALLOC(array, sizeof(float) * 16);
    CHECK(((unsigned long)array % 32lu) == 0);
    MOLEC_FREE(array);

    MOLEC_MALLOC(array, sizeof(float) * 32);
    CHECK(((unsigned long)array % 32lu) == 0);
    MOLEC_FREE(array);

    MOLEC_MALLOC(array, sizeof(float) * 64);
    CHECK(((unsigned long)array % 32lu) == 0);
    MOLEC_FREE(array);
#endif
}
