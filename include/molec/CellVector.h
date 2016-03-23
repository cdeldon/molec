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

#define MOLEC_DIST_PARALL 1.0
#define MOLEC_DIST_DIAG_1 0.70710678118 // 1/sqrt(2)
#define MOLEC_DIST_DIAG_2 0.57735026919 // 1/sqrt(3)
#define MOLEC_DIST_NULLL 0.0

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

