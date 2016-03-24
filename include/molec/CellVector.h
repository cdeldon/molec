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

#ifndef MOLEC_CELL_NORM_H
#define MOLEC_CELL_NORM_H

#include <molec/CellListParam.h>

//  The cell indices are arranged as shown in the illustration below
//  Note that the cell list is traversed with slowest running index 'z', and
//  fastest running index 'x'.
//
//                     _________________________
//                    / _____________________  /|
//                   / / ___________________/ / |
//                  / / /| |               / /  |
//                 / / / | |              / / . |
//                / / /| | |             / / /| |
//               / / / | | |            / / / | |  ^
//              / / /  | | |           / / /| | |  | Z direction
//             / /_/__________________/ / / | | |  |
//            /________________________/ /  | | |
//            | ______________________ | |  | | |
//            | | |    | | |_________| | |__| | |
//            | | |    | |___________| | |____| |
//            | | |   / / ___________| | |_  / /
//            | | |  / / /           | | |/ / /
//            | | | / / /            | | | / /
//            | | |/ / /             | | |/ /
//            | | | / /              | | ' /   ^
//            | | |/_/_______________| |  /   / Y direction
//            | |____________________| | /   /
//            |________________________|/
//                  
//                  --> X direction

#define MOLEC_DIST_PARALL 1.0
#define MOLEC_DIST_DIAG_1 0.70710678118 // 1/sqrt(2)
#define MOLEC_DIST_DIAG_2 0.57735026919 // 1/sqrt(3)
#define MOLEC_DIST_NULLL 0.0

/**
 * @brief Contains normal vectors between two cells
 *
 * The normal vector between two cells is only dependend on the
 * relative position of the two cells, i.e. from the triplet (dx,dy,dz)
 * which belongs to {-1,0,1}x{-1,0,1}x{-1,0,1}.
 *
 * The normal vectors in @c molec_CellLookupTable are ordered in such a way
 * that traversing the neighbour cells as done in @c molec_force_cellList
 * corresponds to iterating the normal vectors in a descending order:
 *
 * @code
 * {-1,-1,-1} --> molec_CellLookupTable[0];
 * {0 ,-1,-1} --> molec_CellLookupTable[1];
 * ...
 * {1 , 1, 1} --> molec_CellLookupTable[26];
 * @endcode
 *
 */
const static Real molec_CellLookupTable[27][3] =
     {
       {-MOLEC_DIST_DIAG_2, -MOLEC_DIST_DIAG_2, -MOLEC_DIST_DIAG_2},
       { MOLEC_DIST_NULLL,  -MOLEC_DIST_DIAG_1, -MOLEC_DIST_DIAG_1},
       { MOLEC_DIST_DIAG_2, -MOLEC_DIST_DIAG_2, -MOLEC_DIST_DIAG_2},

       {-MOLEC_DIST_DIAG_1,  MOLEC_DIST_NULLL,  -MOLEC_DIST_DIAG_1},
       { MOLEC_DIST_NULLL,   MOLEC_DIST_NULLL,  -MOLEC_DIST_PARALL},
       { MOLEC_DIST_DIAG_1,  MOLEC_DIST_NULLL,  -MOLEC_DIST_DIAG_1},

       {-MOLEC_DIST_DIAG_2,  MOLEC_DIST_DIAG_2, -MOLEC_DIST_DIAG_2},
       { MOLEC_DIST_NULLL,   MOLEC_DIST_DIAG_1, -MOLEC_DIST_DIAG_1},
       { MOLEC_DIST_DIAG_2,  MOLEC_DIST_DIAG_2, -MOLEC_DIST_DIAG_2},



       {-MOLEC_DIST_DIAG_1, -MOLEC_DIST_DIAG_1,  MOLEC_DIST_NULLL},
       { MOLEC_DIST_NULLL,  -MOLEC_DIST_PARALL,  MOLEC_DIST_NULLL},
       { MOLEC_DIST_DIAG_1, -MOLEC_DIST_DIAG_1,  MOLEC_DIST_NULLL},

       {-MOLEC_DIST_PARALL, -MOLEC_DIST_NULLL,   MOLEC_DIST_NULLL},
       { MOLEC_DIST_NULLL,  -MOLEC_DIST_NULLL,   MOLEC_DIST_NULLL},
       { MOLEC_DIST_PARALL, -MOLEC_DIST_NULLL,   MOLEC_DIST_NULLL},

       {-MOLEC_DIST_DIAG_1,  MOLEC_DIST_DIAG_1,  MOLEC_DIST_NULLL},
       { MOLEC_DIST_NULLL,   MOLEC_DIST_PARALL,  MOLEC_DIST_NULLL},
       { MOLEC_DIST_DIAG_1,  MOLEC_DIST_DIAG_1,  MOLEC_DIST_NULLL},


       {-MOLEC_DIST_DIAG_2, -MOLEC_DIST_DIAG_2,  MOLEC_DIST_DIAG_2},
       { MOLEC_DIST_NULLL,  -MOLEC_DIST_DIAG_1,  MOLEC_DIST_DIAG_1},
       { MOLEC_DIST_DIAG_2, -MOLEC_DIST_DIAG_2,  MOLEC_DIST_DIAG_2},

       {-MOLEC_DIST_DIAG_1,  MOLEC_DIST_NULLL,   MOLEC_DIST_DIAG_1},
       { MOLEC_DIST_NULLL,   MOLEC_DIST_NULLL,   MOLEC_DIST_PARALL},
       { MOLEC_DIST_DIAG_1,  MOLEC_DIST_NULLL,   MOLEC_DIST_DIAG_1},

       {-MOLEC_DIST_DIAG_2,  MOLEC_DIST_DIAG_2,  MOLEC_DIST_DIAG_2},
       { MOLEC_DIST_NULLL,   MOLEC_DIST_DIAG_1,  MOLEC_DIST_DIAG_1},
       { MOLEC_DIST_DIAG_2,  MOLEC_DIST_DIAG_2,  MOLEC_DIST_DIAG_2}
      };

/**
 * @brief Return the normal vector between 2 cells
 *
 * Computes the normalizing vector of cells and number with all indices
 */
MOLEC_INLINE void
molec_cell_vector(int idx_x, int idx_y, int idx_z, int n_idx_x, int n_idx_y, int n_idx_z)
{
}

#endif

