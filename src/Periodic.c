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

#include <molec/Parameter.h>
#include <molec/Periodic.h>
#include <math.h>

void molec_periodic_refrence(float* x, const int N, const float L)
{
    for(int i = 0; i < N; ++i)
    {
        x[i] = fmod(x[i] + 1e-7, L);
        if(x[i] < 0)
            x[i] += L;
    }
}

void molec_periodic_close(float* x, const int N, const float L)
{
    for (int i = 0; i < N; ++i)
    {
        if (x[i] < 0)
            x[i] += L;
        else if (x[i] > L)
            x[i] -= L;
    }
}

/** Assume that under a not too large timestep,
 * particles can only lie in [-L, 2L]
 */
void molec_periodic_close4(float* x, const int N, const float L)
{
    float x1, x2, x3, x4;

    int i;
    for(i = 0; i < N; i += 4)
    {
        x1 = x[i];
        x2 = x[i + 1];
        x3 = x[i + 2];
        x4 = x[i + 3];

        float is_low1 = x1 < 0;
        float is_high1 = x1 > L;

        float is_low2 = x2 < 0;
        float is_high2 = x2 > L;

        float is_low3 = x3 < 0;
        float is_high3 = x3 > L;

        float is_low4 = x4 < 0;
        float is_high4 = x4 > L;

        x1 = x1 + L * (is_low1 - is_high1);
        x2 = x2 + L * (is_low2 - is_high2);
        x3 = x3 + L * (is_low3 - is_high3);
        x4 = x4 + L * (is_low4 - is_high4);

        x[i] = x1;
        x[i + 1] = x2;
        x[i + 2] = x3;
        x[i + 3] = x4;
    }

    for(int j = i-4 ; j < N; ++j)
    {
        x1 = x[j];

        float is_low = x1 < 0;
        float is_high = x1 > L;

        x1 = x1 + L * (is_low - is_high);

        x[j] = x1;
    }
}

