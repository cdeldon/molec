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

#include <molec/Parameter.h>
#include <molec/Periodic.h>
#include <math.h>

void molec_periodic_refrence(Real* x, const int N)
{
    const Real L = molec_parameter->L;

    for(int i = 0; i < N; ++i)
    {
        x[i] = fmod(x[i], L);
        if(x[i] < 0)
            x[i] += L;
    }
}

