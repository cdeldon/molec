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


#ifndef MOLEC_QUADRANT_H
#define MOLEC_QUADRANT_H

typedef struct molec_Quadrant
{
    float* x;
    float* y;
    float* z;

    float* f_x;
    float* f_y;
    float* f_z;

    int N;

} molec_Quadrant_t;

#endif // MOLEC_QUADRANT_H
