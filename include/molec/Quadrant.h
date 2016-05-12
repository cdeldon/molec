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

#include <molec/Common.h>
#include <molec/Parameter.h>
#include <molec/Simulation.h>

#include <string.h>

#ifndef MOLEC_QUADRANT_H
#define MOLEC_QUADRANT_H

typedef struct molec_Quadrant
{
    float* x;
    float* y;
    float* z;

    float* v_x;
    float* v_y;
    float* v_z;

    float* f_x;
    float* f_y;
    float* f_z;

    int N;

} molec_Quadrant_t;


/**
 * Calculate distance between x and y taking periodic boundaries into account
 *
 * @param      x  coordinate of first particle
 * @param      y  coordinate of second particle
 * @param      L  bounding box size
 * @param      L  bounding box size * 0.5
 *
 * @return     distance between two particles taking into
 *             account periodicity of bounding box
 */
MOLEC_INLINE float distL2(float x, float y, float L, float L2)
{
    float r = x - y;
    if(r < -L2)
        r += L;
    else if(r > L2)
        r -= L;
    return r;
}

/**
 * Calculate the positive modulo between two integers, used for periodic BC
 */
MOLEC_INLINE int mod(int b, int m)
{
    return (b + m) % m;
}

/**
 * Build cellList neighbor structure, where neighbor_cells[i] is a 27-integer array
 * holding the indices of the neighbors of cell i
 */
void molec_build_cell_neighbors(int** neighbor_cells,
                                molec_CellList_Parameter_t cellList_parameter);

/**
  * Initializes the quadrant array by allocating memory and copying the data from the simulation
  * SOA. Always uses in pair with @c molec_quadrants_finalize()
  */
molec_Quadrant_t* molec_quadrant_init(const int N,
                                      molec_CellList_Parameter_t cellList_parameter,
                                      molec_Simulation_SOA_t* sim);

/**
 * Writes back data from the quadrant data structure to the simulation SOA. Memory is freed by this
 * routine. Always to use in pair with @c molec_quadrant_init()
 */
void molec_quadrants_finalize(molec_Quadrant_t* quadrants,
                              molec_CellList_Parameter_t cellList_parameter,
                              molec_Simulation_SOA_t* sim);

#endif // MOLEC_QUADRANT_H
