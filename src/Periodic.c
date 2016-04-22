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

void molec_periodic_refrence(float* x, const int N)
{
    const float L = molec_parameter->L;

    for(int i = 0; i < N; ++i)
    {
        x[i] = fmod(x[i], L);
        if(x[i] < 0)
            x[i] += L;
    }
}

void molec_periodic_v3(float* x, const int N)
{
    const float L = molec_parameter->L;
    int rep;
    for(rep = 0; rep < 10; ++rep)
    {
        for(int i = 0; i < N; ++i)
        {
            if(x[i] < 0)
                x[i] += L;
            if(x[i] > L)
                x[i] -= L;
        }
    }
}

void molec_periodic_v4(float* x, const int N)
{
    const float L = molec_parameter->L;

    for(int i = 0; i < N; ++i)
    {
        if(x[i] < 0)
            while(x[i] < 0)
                x[i] += L;
        if(x[i] > L)
            while(x[i] > L)
                x[i] -= L;
    }
}

void molec_periodic_v5(float* x, const int N)
{
    const float L = molec_parameter->L;

    float x1, x2;
    for(int i = 0; i < N; i += 2)
    {
        x1 = x[i];
        x2 = x[i + 1];

        if(x1 < 0)
            while(x1 < 0)
                x1 += L;
        if(x2 < 0)
            while(x2 < 0)
                x2 += L;

        if(x1 > L)
            while(x1 > L)
                x1 -= L;
        if(x2 > L)
            while(x2 > L)
                x2 -= L;
        x[i] = x1;
        x[i + 1] = x2;
    }
}
