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
#include <molec/CellVector.h>

const static float molec_ref_CellLookupTable[27][3] =
     {
       {-1,-1,-1},
       { 0,-1,-1},
       { 1,-1,-1},

       {-1, 0,-1},
       { 0, 0,-1},
       { 1, 0,-1},

       {-1, 1,-1},
       { 0, 1,-1},
       { 1, 1,-1},


       {-1,-1, 0},
       { 0,-1, 0},
       { 1,-1, 0},

       {-1, 0, 0},
       { 0, 0, 0},
       { 1, 0, 0},

       {-1, 1, 0},
       { 0, 1, 0},
       { 1, 1, 0},


       {-1,-1, 1},
       { 0,-1, 1},
       { 1,-1, 1},

       {-1, 0, 1},
       { 0, 0, 1},
       { 1, 0, 1},

       {-1, 1, 1},
       { 0, 1, 1},
       { 1, 1, 1}
      };


MOLEC_INLINE void cross(const float a[3], const float b[3])
{
    float cross[3];
    float zero[3] = {0., 0., 0.};
    cross[0] = a[1] * b[2] - a[2] * b[1];
    cross[1] = a[2] * b[0] - a[0] * b[2];
    cross[2] = a[0] * b[1] - a[1] * b[0];

    ALLCLOSE_FLOAT(cross, zero, 3, MOLEC_ATOL, MOLEC_RTOL)
}

/**
 * @brief Test CellVector directions
 *
 * Tests each direction of @c molec_CellLookupTable
 * checking that the cross product with a reference direction gives a zero vector
 */
TEST_CASE(molec_UnittestCellVectorDirections)
{
    // Check each vector direction
    for(int i = 0; i < 27; ++i)
        cross(molec_ref_CellLookupTable[i], molec_CellLookupTable[i]);
}

