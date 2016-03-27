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

#include <molec/Sort.h>
#include <molec/Parameter.h>
#include <string.h>

typedef struct molec_Sort_Pair {
    int key, value;
} molec_Sort_Pair_t;

int molec_compare(const void *x, const void *y)
{
    Real a = ((const molec_Sort_Pair_t*)x)->key;
    Real b = ((const molec_Sort_Pair_t*)y)->key;
    return a < b ? -1 : a > b;
}

void molec_sort_qsort(molec_Simulation_SOA_t* sim)
{
    // Local aliases
    Real *x = sim->x;
    Real *y = sim->y;
    Real *z = sim->z;

    Real *v_x = sim->v_x;
    Real *v_y = sim->v_y;
    Real *v_z = sim->v_z;

    const int N = molec_parameter->N;

    molec_Sort_Pair_t *key;
    MOLEC_MALLOC(key, N * sizeof(molec_Sort_Pair_t));

    for(int i=0;i<N;++i){
        key[i].key  = x[i];
        key[i].value=i;
    }
    qsort(key, N, sizeof(molec_Sort_Pair_t), molec_compare);
    Real *x_temp, *y_temp, *z_temp;
    Real *v_x_temp, *v_y_temp, *v_z_temp;

    MOLEC_MALLOC(x_temp, N * sizeof(Real));
    MOLEC_MALLOC(y_temp, N * sizeof(Real));
    MOLEC_MALLOC(z_temp, N * sizeof(Real));
    MOLEC_MALLOC(v_x_temp, N * sizeof(Real));
    MOLEC_MALLOC(v_y_temp, N * sizeof(Real));
    MOLEC_MALLOC(v_z_temp, N * sizeof(Real));

    memcpy(x_temp, x, N * sizeof(Real));
    memcpy(y_temp, y, N * sizeof(Real));
    memcpy(z_temp, z, N * sizeof(Real));
    memcpy(v_x_temp, v_x, N * sizeof(Real));
    memcpy(v_y_temp, v_y, N * sizeof(Real));
    memcpy(v_z_temp, v_z, N * sizeof(Real));

    for(int i=0;i<N;++i){
        x[i] = x_temp[key[i].value];
        y[i] = y_temp[key[i].value];
        z[i] = z_temp[key[i].value];
        v_x[i] = v_x_temp[key[i].value];
        v_y[i] = v_y_temp[key[i].value];
        v_z[i] = v_z_temp[key[i].value];
    }

    MOLEC_FREE(x_temp);
    MOLEC_FREE(y_temp);
    MOLEC_FREE(z_temp);
    MOLEC_FREE(v_x_temp);
    MOLEC_FREE(v_y_temp);
    MOLEC_FREE(v_z_temp);

    MOLEC_FREE(key);

    sim->x = x;
    sim->y = y;
    sim->z = z;
    sim->v_x = v_x;
    sim->v_y = v_y;
    sim->v_z = v_z;
}
