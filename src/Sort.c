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
    float pair1_key1 = (float) ((molec_Sort_Pair_t*) pair1)->key1;
    float pair2_key1 = (float) ((molec_Sort_Pair_t*) pair2)->key1;

    float pair1_key2 = (float) ((molec_Sort_Pair_t*) pair1)->key2;
    float pair2_key2 = (float) ((molec_Sort_Pair_t*) pair2)->key2;

    int pair1_value = (int)((molec_Sort_Pair_t*) pair1)->value;
    int pair2_value = (int)((molec_Sort_Pair_t*) pair2)->value;

    if(pair1_key1 < pair2_key1)
        return -1;
    else if(pair1_key1 > pair2_key1)
        return 1;
    else // pair1_key1 == pair2_key1
    {
        if(pair1_key2 < pair2_key2)
            return -1;
        else
            return 1;
    }

}

void molec_sort_qsort(molec_Simulation_SOA_t* sim)
{
    // Local aliases
    float *x = sim->x;
    float *y = sim->y;
    float *z = sim->z;

    float *v_x = sim->v_x;
    float *v_y = sim->v_y;
    float *v_z = sim->v_z;

    const int N = molec_parameter->N;

    molec_Sort_Pair_t* key;
    MOLEC_MALLOC(key, N * sizeof(molec_Sort_Pair_t));

    for(int i = 0; i < N; ++i)
    {
        key[i].key1  = x[i];
        key[i].key2  = y[i];
        key[i].value = i;
    }
    
    qsort(key, N, sizeof(molec_Sort_Pair_t), molec_compare);
    
    float *x_temp, *y_temp, *z_temp;
    float *v_x_temp, *v_y_temp, *v_z_temp;

    MOLEC_MALLOC(x_temp, N * sizeof(float));
    MOLEC_MALLOC(y_temp, N * sizeof(float));
    MOLEC_MALLOC(z_temp, N * sizeof(float));
    MOLEC_MALLOC(v_x_temp, N * sizeof(float));
    MOLEC_MALLOC(v_y_temp, N * sizeof(float));
    MOLEC_MALLOC(v_z_temp, N * sizeof(float));

    memcpy(x_temp, x, N * sizeof(float));
    memcpy(y_temp, y, N * sizeof(float));
    memcpy(z_temp, z, N * sizeof(float));
    memcpy(v_x_temp, v_x, N * sizeof(float));
    memcpy(v_y_temp, v_y, N * sizeof(float));
    memcpy(v_z_temp, v_z, N * sizeof(float));

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
    float *x = sim->x;
    float *y = sim->y;
    float *z = sim->z;

    float *f_x = sim->f_x;
    float *f_y = sim->f_y;
    float *f_z = sim->f_z;


    const int N = molec_parameter->N;

    molec_Sort_Pair_t* key;
    MOLEC_MALLOC(key, N * sizeof(molec_Sort_Pair_t));

    for(int i = 0; i < N; ++i)
    {
        key[i].key1  = x[i];
        key[i].key2  = y[i];
        key[i].value = i;
    }

    qsort(key, N, sizeof(molec_Sort_Pair_t), molec_compare);

    float *x_temp, *y_temp, *z_temp;
    float *f_x_temp, *f_y_temp, *f_z_temp;

    MOLEC_MALLOC(x_temp, N * sizeof(float));
    MOLEC_MALLOC(y_temp, N * sizeof(float));
    MOLEC_MALLOC(z_temp, N * sizeof(float));
    MOLEC_MALLOC(f_x_temp, N * sizeof(float));
    MOLEC_MALLOC(f_y_temp, N * sizeof(float));
    MOLEC_MALLOC(f_z_temp, N * sizeof(float));

    memcpy(x_temp, x, N * sizeof(float));
    memcpy(y_temp, y, N * sizeof(float));
    memcpy(z_temp, z, N * sizeof(float));
    memcpy(f_x_temp, f_x, N * sizeof(float));
    memcpy(f_y_temp, f_y, N * sizeof(float));
    memcpy(f_z_temp, f_z, N * sizeof(float));

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
