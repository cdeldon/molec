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

int molec_compare(const void* pair1, const void* pair2)
{
    Real pair1_key = (Real) ((molec_Sort_Pair_t*) pair1)->key;
    Real pair2_key = (Real) ((molec_Sort_Pair_t*) pair2)->key;
    
    // Equalilty will never happen
    if(pair1_key < pair2_key)
        return -1;
    else if(pair1_key > pair2_key)
        return 1;
    else
        return 0;
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

    molec_Sort_Pair_t* key;
    MOLEC_MALLOC(key, N * sizeof(molec_Sort_Pair_t));

    for(int i = 0; i < N; ++i)
    {
        key[i].key  = x[i];
        key[i].value = i;
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

    for(int i=0; i < N; ++i)
    {
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


void molec_sort_qsort_forces(molec_Simulation_SOA_t* sim)
{
    // Local aliases
    Real *x = sim->x;
    Real *y = sim->y;
    Real *z = sim->z;

    Real *f_x = sim->f_x;
    Real *f_y = sim->f_y;
    Real *f_z = sim->f_z;


    const int N = molec_parameter->N;

    molec_Sort_Pair_t* key;
    MOLEC_MALLOC(key, N * sizeof(molec_Sort_Pair_t));

    for(int i = 0; i < N; ++i)
    {
        key[i].key  = x[i];
        key[i].value = i;
    }

    qsort(key, N, sizeof(molec_Sort_Pair_t), molec_compare);

    Real *x_temp, *y_temp, *z_temp;
    Real *f_x_temp, *f_y_temp, *f_z_temp;

    MOLEC_MALLOC(x_temp, N * sizeof(Real));
    MOLEC_MALLOC(y_temp, N * sizeof(Real));
    MOLEC_MALLOC(z_temp, N * sizeof(Real));
    MOLEC_MALLOC(f_x_temp, N * sizeof(Real));
    MOLEC_MALLOC(f_y_temp, N * sizeof(Real));
    MOLEC_MALLOC(f_z_temp, N * sizeof(Real));

    memcpy(x_temp, x, N * sizeof(Real));
    memcpy(y_temp, y, N * sizeof(Real));
    memcpy(z_temp, z, N * sizeof(Real));
    memcpy(f_x_temp, f_x, N * sizeof(Real));
    memcpy(f_y_temp, f_y, N * sizeof(Real));
    memcpy(f_z_temp, f_z, N * sizeof(Real));

    for(int i=0; i < N; ++i)
    {
        x[i] = x_temp[key[i].value];
        y[i] = y_temp[key[i].value];
        z[i] = z_temp[key[i].value];
        f_x[i] = f_x_temp[key[i].value];
        f_y[i] = f_y_temp[key[i].value];
        f_z[i] = f_z_temp[key[i].value];
    }

    MOLEC_FREE(x_temp);
    MOLEC_FREE(y_temp);
    MOLEC_FREE(z_temp);
    MOLEC_FREE(f_x_temp);
    MOLEC_FREE(f_y_temp);
    MOLEC_FREE(f_z_temp);

    MOLEC_FREE(key);

    sim->x = x;
    sim->y = y;
    sim->z = z;
    sim->f_x = f_x;
    sim->f_y = f_y;
    sim->f_z = f_z;
}
