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

#include <molec/Force.h>
#include <molec/Quadrant.h>
#include <molec/Timer.h>
#include <immintrin.h>

/**
 * @brief neighbor_cells[i] constains the 27 indices of its neighbor cells
 */
static int** neighbor_cells;

/**************************************************************************************************/

void molec_quadrant_neighbor_interaction(molec_Quadrant_t q, molec_Quadrant_t q_n, float* Epot_)
{
    const float sigLJ = molec_parameter->sigLJ;
    const float epsLJ = molec_parameter->epsLJ;
    const float L_x = molec_parameter->L_x;
    const float L_y = molec_parameter->L_y;
    const float L_z = molec_parameter->L_z;
    const float L_x2 = L_x * 0.5f;
    const float L_y2 = L_y * 0.5f;
    const float L_z2 = L_z * 0.5f;
    const float Rcut2 = molec_parameter->Rcut2;

    const int N = q.N;
    const int N_n = q_n.N;

    for(int i = 0; i < N; ++i)
    {
        float xi = q.x[i];
        float yi = q.y[i];
        float zi = q.z[i];

        float f_xi = q.f_x[i];
        float f_yi = q.f_y[i];
        float f_zi = q.f_z[i];

        for(int j = 0; j < N_n; ++j)
        {
            const float xij = distL2(xi, q_n.x[j], L_x, L_x2);
            const float yij = distL2(yi, q_n.y[j], L_y, L_y2);
            const float zij = distL2(zi, q_n.z[j], L_z, L_z2);

            const float r2 = xij * xij + yij * yij + zij * zij;

            if(r2 < Rcut2)
            {
                // V(s) = 4 * eps * (s^12 - s^6) with  s = sig/r
                const float s2 = (sigLJ * sigLJ) / r2;
                const float s6 = s2 * s2 * s2;

                *Epot_ += 4 * epsLJ*(s6 * s6 - s6);

                const float fr = 24 * epsLJ / r2 * (2 * s6 * s6 - s6);

                f_xi += fr * xij;
                f_yi += fr * yij;
                f_zi += fr * zij;

                q_n.f_x[j] -= fr * xij;
                q_n.f_y[j] -= fr * yij;
                q_n.f_z[j] -= fr * zij;
            }
        }
        q.f_x[i] = f_xi;
        q.f_y[i] = f_yi;
        q.f_z[i] = f_zi;
    }
}

void molec_quadrant_self_interaction(molec_Quadrant_t q, float* Epot_)
{
    const float sigLJ = molec_parameter->sigLJ;
    const float epsLJ = molec_parameter->epsLJ;
    const float L_x = molec_parameter->L_x;
    const float L_y = molec_parameter->L_y;
    const float L_z = molec_parameter->L_z;
    const float L_x2 = L_x * 0.5f;
    const float L_y2 = L_y * 0.5f;
    const float L_z2 = L_z * 0.5f;
    const float Rcut2 = molec_parameter->Rcut2;

    for(int i = 0; i < q.N; ++i)
    {
        float xi = q.x[i];
        float yi = q.y[i];
        float zi = q.z[i];

        float f_xi = q.f_x[i];
        float f_yi = q.f_y[i];
        float f_zi = q.f_z[i];

        for(int j = i + 1; j < q.N; ++j)
        {

            const float xij = distL2(xi, q.x[j], L_x, L_x2);
            const float yij = distL2(yi, q.y[j], L_y, L_y2);
            const float zij = distL2(zi, q.z[j], L_z, L_z2);

            const float r2 = xij * xij + yij * yij + zij * zij;

            if(r2 < Rcut2)
            {
                // V(s) = 4 * eps * (s^12 - s^6) with  s = sig/r
                const float s2 = (sigLJ * sigLJ) / r2;
                const float s6 = s2 * s2 * s2;

                *Epot_ += 4 * epsLJ*(s6 * s6 - s6);

                const float fr = 24 * epsLJ / r2 * (2 * s6 * s6 - s6);

                f_xi += fr * xij;
                f_yi += fr * yij;
                f_zi += fr * zij;

                q.f_x[j] -= fr * xij;
                q.f_y[j] -= fr * yij;
                q.f_z[j] -= fr * zij;
            }
        }
        q.f_x[i] = f_xi;
        q.f_y[i] = f_yi;
        q.f_z[i] = f_zi;
    }
}

void molec_force_quadrant(molec_Simulation_SOA_t* sim, float* Epot, const int N)
{
    assert(molec_parameter);

    molec_CellList_Parameter_t cellList_parameter = molec_parameter->cellList;

    float Epot_ = 0;

    // Build the neighbor_cell array only if not initialized before
    if(neighbor_cells == NULL)
    {
        MOLEC_MALLOC(neighbor_cells, cellList_parameter.N * sizeof(int*));
        for(int i = 0; i < cellList_parameter.N; ++i)
            MOLEC_MALLOC(neighbor_cells[i], 27 * sizeof(int));

        molec_build_cell_neighbors(neighbor_cells, cellList_parameter);
    }

    molec_Quadrant_t* quadrants = molec_quadrant_init(N, cellList_parameter, sim);


    // loop over the quadrants
    for(int idx = 0; idx < cellList_parameter.N; ++idx)
    {

        // loop over all the neighbors
        for(int neighbor_cell = 0; neighbor_cell < 27; ++neighbor_cell)
        {
            int n_idx = neighbor_cells[idx][neighbor_cell];

            if(idx > n_idx)
                molec_quadrant_neighbor_interaction(quadrants[idx], quadrants[n_idx], &Epot_);

            else if(idx == n_idx)
                molec_quadrant_self_interaction(quadrants[idx], &Epot_);
        }
    }

    molec_quadrants_finalize(quadrants, cellList_parameter, sim);

    *Epot = Epot_;
}

/**************************************************************************************************/

/**
 * @brief neighbor_cells_ghost[i] constains the 27 indices of its neighbor cells
 */
static int** neighbor_cells_ghost;

MOLEC_INLINE float dist_ghost(float x, float y)
{
    float r = x - y;
    return r;
}

void molec_quadrant_neighbor_interaction_ghost(molec_Quadrant_t q,
                                               molec_Quadrant_t q_n,
                                               float* Epot_)
{
    const float sigLJ = molec_parameter->sigLJ;
    const float epsLJ = molec_parameter->epsLJ;

    const float Rcut2 = molec_parameter->Rcut2;

    const int N = q.N;
    const int N_n = q_n.N;

    for(int i = 0; i < N; ++i)
    {
        float xi = q.x[i];
        float yi = q.y[i];
        float zi = q.z[i];

        float f_xi = q.f_x[i];
        float f_yi = q.f_y[i];
        float f_zi = q.f_z[i];

        for(int j = 0; j < N_n; ++j)
        {
            const float xij = dist_ghost(xi, q_n.x[j]);
            const float yij = dist_ghost(yi, q_n.y[j]);
            const float zij = dist_ghost(zi, q_n.z[j]);

            const float r2 = xij * xij + yij * yij + zij * zij;

            if(r2 < Rcut2)
            {
                // V(s) = 4 * eps * (s^12 - s^6) with  s = sig/r
                const float s2 = (sigLJ * sigLJ) / r2;
                const float s6 = s2 * s2 * s2;

                *Epot_ += 4 * epsLJ*(s6 * s6 - s6);

                const float fr = 24 * epsLJ / r2 * (2 * s6 * s6 - s6);

                f_xi += fr * xij;
                f_yi += fr * yij;
                f_zi += fr * zij;

                q_n.f_x[j] -= fr * xij;
                q_n.f_y[j] -= fr * yij;
                q_n.f_z[j] -= fr * zij;
            }
        }
        q.f_x[i] = f_xi;
        q.f_y[i] = f_yi;
        q.f_z[i] = f_zi;
    }
}

void molec_quadrant_self_interaction_ghost(molec_Quadrant_t q, float* Epot_)
{
    const float sigLJ = molec_parameter->sigLJ;
    const float epsLJ = molec_parameter->epsLJ;

    const float Rcut2 = molec_parameter->Rcut2;

    for(int i = 0; i < q.N; ++i)
    {
        float xi = q.x[i];
        float yi = q.y[i];
        float zi = q.z[i];

        float f_xi = q.f_x[i];
        float f_yi = q.f_y[i];
        float f_zi = q.f_z[i];

        for(int j = i + 1; j < q.N; ++j)
        {
            const float xij = dist_ghost(xi, q.x[j]);
            const float yij = dist_ghost(yi, q.y[j]);
            const float zij = dist_ghost(zi, q.z[j]);

            const float r2 = xij * xij + yij * yij + zij * zij;

            if(r2 < Rcut2)
            {
                // V(s) = 4 * eps * (s^12 - s^6) with  s = sig/r
                const float s2 = (sigLJ * sigLJ) / r2;
                const float s6 = s2 * s2 * s2;

                *Epot_ += 4 * epsLJ*(s6 * s6 - s6);

                const float fr = 24 * epsLJ / r2 * (2 * s6 * s6 - s6);

                f_xi += fr * xij;
                f_yi += fr * yij;
                f_zi += fr * zij;

                q.f_x[j] -= fr * xij;
                q.f_y[j] -= fr * yij;
                q.f_z[j] -= fr * zij;
            }
        }
        q.f_x[i] = f_xi;
        q.f_y[i] = f_yi;
        q.f_z[i] = f_zi;
    }
}

void molec_force_quadrant_ghost(molec_Simulation_SOA_t* sim, float* Epot, const int N)
{
    assert(molec_parameter);

    molec_CellList_Parameter_t cellList_parameter = molec_parameter->cellList;

    const int N_x_ghost = cellList_parameter.N_x + 2;
    const int N_y_ghost = cellList_parameter.N_y + 2;
    const int N_z_ghost = cellList_parameter.N_z + 2;

    float Epot_ = 0;

    MOLEC_MEASUREMENT_CELL_CONSTRUCTION_START();
    // Build the neighbor_cell array only if not initialized before
    if(neighbor_cells_ghost == NULL)
    {
        const int N_ghost = N_x_ghost * N_y_ghost * N_z_ghost;
        MOLEC_MALLOC(neighbor_cells_ghost, N_ghost * sizeof(int*));
        // allocate the neighors for the inner quadrants
        for(int idx_z = 1; idx_z <= cellList_parameter.N_z; ++idx_z)
            for(int idx_y = 1; idx_y <= cellList_parameter.N_y; ++idx_y)
                for(int idx_x = 1; idx_x <= cellList_parameter.N_x; ++idx_x)
                {
                    // linear index of cell
                    const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
                    MOLEC_MALLOC(neighbor_cells_ghost[idx], 27 * sizeof(int));
                }

        molec_build_cell_neighbors_ghost(neighbor_cells_ghost, cellList_parameter);
    }

    molec_Quadrant_t* quadrants = molec_quadrant_init_ghost(N, cellList_parameter, sim);

    MOLEC_MEASUREMENT_CELL_CONSTRUCTION_STOP();

    // loop over the internal quadrants
    for(int idx_z = 1; idx_z <= cellList_parameter.N_z; ++idx_z)
        for(int idx_y = 1; idx_y <= cellList_parameter.N_y; ++idx_y)
            for(int idx_x = 1; idx_x <= cellList_parameter.N_x; ++idx_x)
            {
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
                // loop over all the neighbors
                for(int neighbor_cell = 0; neighbor_cell < 27; ++neighbor_cell)
                {
                    int n_idx = neighbor_cells_ghost[idx][neighbor_cell];

                    if(quadrants[idx].idx < quadrants[n_idx].idx)
                        molec_quadrant_neighbor_interaction_ghost(quadrants[idx], quadrants[n_idx],
                                                                  &Epot_);

                    else if(idx == n_idx)
                        molec_quadrant_self_interaction_ghost(quadrants[idx], &Epot_);
                }
            }

    molec_quadrants_finalize_ghost(quadrants, cellList_parameter, sim);

    *Epot = Epot_;
}


/**************************************************************************************************/

void molec_quadrant_neighbor_interaction_avx(molec_Quadrant_t q, molec_Quadrant_t q_n, float* Epot_)
{
    const __m256 sigLJ = _mm256_set1_ps(molec_parameter->sigLJ);
    const __m256 epsLJ = _mm256_set1_ps(molec_parameter->epsLJ);

    const __m256 Rcut2 = _mm256_set1_ps(molec_parameter->Rcut2);

    const int N = q.N;
    const int N_n = q_n.N_pad;

    __m256 Epot8 = _mm256_setzero_ps();
    __m256 _1 = _mm256_set1_ps(1.f);
    __m256 _2 = _mm256_set1_ps(2.f);
    __m256 _24epsLJ = _mm256_mul_ps(_mm256_set1_ps(24.f), epsLJ);

    for(int i = 0; i < N; ++i)
    {
        const __m256 xi = _mm256_set1_ps(q.x[i]);
        const __m256 yi = _mm256_set1_ps(q.y[i]);
        const __m256 zi = _mm256_set1_ps(q.z[i]);

        __m256 f_xi = _mm256_setzero_ps();
        __m256 f_yi = _mm256_setzero_ps();
        __m256 f_zi = _mm256_setzero_ps();

        for(int j = 0; j < N_n; j += 8)
        {
            // load coordinates and fores into AVX vectors
            const __m256 xj = _mm256_load_ps(&q_n.x[j]);
            const __m256 yj = _mm256_load_ps(&q_n.y[j]);
            const __m256 zj = _mm256_load_ps(&q_n.z[j]);

            __m256 f_xj = _mm256_load_ps(&q_n.f_x[j]);
            __m256 f_yj = _mm256_load_ps(&q_n.f_y[j]);
            __m256 f_zj = _mm256_load_ps(&q_n.f_z[j]);


            // distance computation
            const __m256 xij = _mm256_sub_ps(xi, xj);
            const __m256 xij2 = _mm256_mul_ps(xij, xij);

            const __m256 yij = _mm256_sub_ps(yi, yj);
            const __m256 yij2 = _mm256_mul_ps(yij, yij);

            const __m256 zij = _mm256_sub_ps(zi, zj);
            const __m256 zij2 = _mm256_mul_ps(zij, zij);

            const __m256 r2 = _mm256_add_ps(_mm256_add_ps(xij2, yij2), zij2);


            // r2 < Rcut2
            const __m256 mask = _mm256_cmp_ps(r2, Rcut2, _CMP_LT_OQ);

            // if( any(r2 < R2) )
            if(_mm256_movemask_ps(mask))
            {
                const __m256 r2inv = _mm256_div_ps(_1, r2);

                const __m256 s2 = _mm256_mul_ps(_mm256_mul_ps(sigLJ, sigLJ), r2inv);
                const __m256 s6 = _mm256_mul_ps(_mm256_mul_ps(s2, s2), s2);
                const __m256 s12 = _mm256_mul_ps(s6, s6);

                const __m256 s12_minus_s6 = _mm256_sub_ps(s12, s6);
                const __m256 two_s12_minus_s6 = _mm256_sub_ps(_mm256_mul_ps(_2, s12), s6);

                Epot8 = _mm256_add_ps(Epot8, _mm256_and_ps(s12_minus_s6, mask));

                const __m256 fr = _mm256_mul_ps(_mm256_mul_ps(_24epsLJ, r2inv), two_s12_minus_s6);
                const __m256 fr_mask = _mm256_and_ps(fr, mask);

                const __m256 fr_x = _mm256_mul_ps(fr_mask, xij);
                const __m256 fr_y = _mm256_mul_ps(fr_mask, yij);
                const __m256 fr_z = _mm256_mul_ps(fr_mask, zij);

                // update forces
                f_xi = _mm256_add_ps(f_xi, fr_x);
                f_yi = _mm256_add_ps(f_yi, fr_y);
                f_zi = _mm256_add_ps(f_zi, fr_z);

                f_xj = _mm256_sub_ps(f_xj, fr_x);
                f_yj = _mm256_sub_ps(f_yj, fr_y);
                f_zj = _mm256_sub_ps(f_zj, fr_z);

                // store back j-forces
                _mm256_store_ps(&q_n.f_x[j], f_xj);
                _mm256_store_ps(&q_n.f_y[j], f_yj);
                _mm256_store_ps(&q_n.f_z[j], f_zj);
            }
        }

        // update i-forces
        float MOLEC_ALIGNAS(32) f_array[8];
        _mm256_store_ps(f_array, f_xi);
        q.f_x[i] += f_array[0] + f_array[1] + f_array[2] + f_array[3] + f_array[4] + f_array[5]
                    + f_array[6] + f_array[7];
        _mm256_store_ps(f_array, f_yi);
        q.f_y[i] += f_array[0] + f_array[1] + f_array[2] + f_array[3] + f_array[4] + f_array[5]
                    + f_array[6] + f_array[7];
        _mm256_store_ps(f_array, f_zi);
        q.f_z[i] += f_array[0] + f_array[1] + f_array[2] + f_array[3] + f_array[4] + f_array[5]
                    + f_array[6] + f_array[7];
    }

    float MOLEC_ALIGNAS(32) E_pot_array[8];
    _mm256_store_ps(E_pot_array, Epot8);

    // perform reduction of potential energy
    *Epot_ += 4
              * molec_parameter->epsLJ*(E_pot_array[0] + E_pot_array[1] + E_pot_array[2]
                                        + E_pot_array[3] + E_pot_array[4] + E_pot_array[5]
                                        + E_pot_array[6] + E_pot_array[7]);
}

void molec_quadrant_self_interaction_avx(molec_Quadrant_t q, float* Epot_)
{
    const float sigLJ = molec_parameter->sigLJ;
        const float epsLJ = molec_parameter->epsLJ;

        const float Rcut2 = molec_parameter->Rcut2;

        for(int i = 0; i < q.N; ++i)
        {
            float xi = q.x[i];
            float yi = q.y[i];
            float zi = q.z[i];

            float f_xi = q.f_x[i];
            float f_yi = q.f_y[i];
            float f_zi = q.f_z[i];

            for(int j = i + 1; j < q.N; ++j)
            {
                const float xij = xi - q.x[j];
                const float yij = yi - q.y[j];
                const float zij = zi - q.z[j];

                const float r2 = xij * xij + yij * yij + zij * zij;

                if(r2 < Rcut2)
                {
                    // V(s) = 4 * eps * (s^12 - s^6) with  s = sig/r
                    const float s2 = (sigLJ * sigLJ) / r2;
                    const float s6 = s2 * s2 * s2;

                    *Epot_ += 4 * epsLJ*(s6 * s6 - s6);

                    const float fr = 24 * epsLJ / r2 * (2 * s6 * s6 - s6);

                    f_xi += fr * xij;
                    f_yi += fr * yij;
                    f_zi += fr * zij;

                    q.f_x[j] -= fr * xij;
                    q.f_y[j] -= fr * yij;
                    q.f_z[j] -= fr * zij;
                }
            }
            q.f_x[i] = f_xi;
            q.f_y[i] = f_yi;
            q.f_z[i] = f_zi;
        }
    /*
    const float sigLJ = molec_parameter->sigLJ;
    const float epsLJ = molec_parameter->epsLJ;

    const float Rcut2 = molec_parameter->Rcut2;

    float Epot = 0.f;

    for(int i = 0; i < q.N; ++i)
    {
        float xi = q.x[i];
        float yi = q.y[i];
        float zi = q.z[i];

        float f_xi = q.f_x[i];
        float f_yi = q.f_y[i];
        float f_zi = q.f_z[i];

        for(int j = i + 1; j < q.N; j+=2)
        {
            const float xij_1 = xi - q.x[j+0];
            const float yij_1 = yi - q.y[j+0];
            const float zij_1 = zi - q.z[j+0];

            const float xij_2 = xi - q.x[j+1];
            const float yij_2 = yi - q.y[j+1];
            const float zij_2 = zi - q.z[j+1];

            const float r2_1 = xij_1 * xij_1 + yij_1 * yij_1 + zij_1 * zij_1;
            const float r2_2 = xij_2 * xij_2 + yij_2 * yij_2 + zij_2 * zij_2;

            int r2_1_smaller_Rcut2 = r2_1 < Rcut2;
            int r2_2_smaller_Rcut2 = r2_2 < Rcut2;

            if(r2_1_smaller_Rcut2 && r2_2_smaller_Rcut2 )
            {
                const float r2_1_inv = 1.f/r2_1;
                const float r2_2_inv = 1.f/r2_2;

                const float s2_1 = (sigLJ * sigLJ) * r2_1_inv;
                const float s2_2 = (sigLJ * sigLJ) * r2_2_inv;
                const float s6_1 = s2_1 * s2_1 * s2_1;
                const float s6_2 = s2_2 * s2_2 * s2_2;
                const float s12_1 = s6_1 * s6_1;
                const float s12_2 = s6_2 * s6_2;

                Epot += (s12_1 - s6_1) + (s12_2 - s6_2);

                const float fr_1 = 24 * epsLJ * r2_1_inv * (2 * s12_1 - s6_1);
                const float fr_2 = 24 * epsLJ * r2_2_inv * (2 * s12_2 - s6_2);

                f_xi += fr_1 * xij_1 + fr_2 * xij_2;
                f_yi += fr_1 * yij_1 + fr_2 * yij_2;
                f_zi += fr_1 * zij_1 + fr_2 * zij_2;

                q.f_x[j + 0] -= fr_1 * xij_1;
                q.f_y[j + 0] -= fr_1 * yij_1;
                q.f_z[j + 0] -= fr_1 * zij_1;

                q.f_x[j + 1] -= fr_2 * xij_2;
                q.f_y[j + 1] -= fr_2 * yij_2;
                q.f_z[j + 1] -= fr_2 * zij_2;
            }
            else if(r2_1_smaller_Rcut2)
            {
                // V(s) = 4 * eps * (s^12 - s^6) with  s = sig/r
                const float s2 = (sigLJ * sigLJ) / r2_1;
                const float s6 = s2 * s2 * s2;

                *Epot_ += 4 * epsLJ*(s6 * s6 - s6);

                const float fr = 24 * epsLJ / r2_1 * (2 * s6 * s6 - s6);

                f_xi += fr * xij_1;
                f_yi += fr * yij_1;
                f_zi += fr * zij_1;

                q.f_x[j] -= fr * xij_1;
                q.f_y[j] -= fr * yij_1;
                q.f_z[j] -= fr * zij_1;
            }
            else if(r2_2_smaller_Rcut2)
            {
                // V(s) = 4 * eps * (s^12 - s^6) with  s = sig/r
                const float s2 = (sigLJ * sigLJ) / r2_2;
                const float s6 = s2 * s2 * s2;

                *Epot_ += 4 * epsLJ*(s6 * s6 - s6);

                const float fr = 24 * epsLJ / r2_2 * (2 * s6 * s6 - s6);

                f_xi += fr * xij_2;
                f_yi += fr * yij_2;
                f_zi += fr * zij_2;

                q.f_x[j+1] -= fr * xij_2;
                q.f_y[j+1] -= fr * yij_2;
                q.f_z[j+1] -= fr * zij_2;
            }
        }
        q.f_x[i] = f_xi;
        q.f_y[i] = f_yi;
        q.f_z[i] = f_zi;
    }

    Epot *= 4 * epsLJ;

    *Epot_ += Epot;
    */
}

void molec_force_quadrant_avx(molec_Simulation_SOA_t* sim, float* Epot, const int N)
{
    assert(molec_parameter);

    molec_CellList_Parameter_t cellList_parameter = molec_parameter->cellList;

    const int N_x_ghost = cellList_parameter.N_x + 2;
    const int N_y_ghost = cellList_parameter.N_y + 2;
    const int N_z_ghost = cellList_parameter.N_z + 2;

    float Epot_ = 0;

    MOLEC_MEASUREMENT_CELL_CONSTRUCTION_START();
    // Build the neighbor_cell array only if not initialized before
    if(neighbor_cells_ghost == NULL)
    {
        const int N_ghost = N_x_ghost * N_y_ghost * N_z_ghost;
        MOLEC_MALLOC(neighbor_cells_ghost, N_ghost * sizeof(int*));
        // allocate the neighors for the inner quadrants
        for(int idx_z = 1; idx_z <= cellList_parameter.N_z; ++idx_z)
            for(int idx_y = 1; idx_y <= cellList_parameter.N_y; ++idx_y)
                for(int idx_x = 1; idx_x <= cellList_parameter.N_x; ++idx_x)
                {
                    // linear index of cell
                    const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
                    MOLEC_MALLOC(neighbor_cells_ghost[idx], 27 * sizeof(int));
                }

        molec_build_cell_neighbors_ghost(neighbor_cells_ghost, cellList_parameter);
    }

    molec_Quadrant_t* quadrants = molec_quadrant_init_ghost(N, cellList_parameter, sim);

    MOLEC_MEASUREMENT_CELL_CONSTRUCTION_STOP();

    // loop over the internal quadrants
    for(int idx_z = 1; idx_z <= cellList_parameter.N_z; ++idx_z)
        for(int idx_y = 1; idx_y <= cellList_parameter.N_y; ++idx_y)
            for(int idx_x = 1; idx_x <= cellList_parameter.N_x; ++idx_x)
            {
                const int idx = idx_x + N_x_ghost * (idx_y + N_y_ghost * idx_z);
                // loop over all the neighbors
                for(int neighbor_cell = 0; neighbor_cell < 27; ++neighbor_cell)
                {
                    int n_idx = neighbor_cells_ghost[idx][neighbor_cell];

                    if(quadrants[idx].idx < quadrants[n_idx].idx)
                        molec_quadrant_neighbor_interaction_avx(quadrants[idx], quadrants[n_idx],
                                                                &Epot_);

                    else if(idx == n_idx)
                        molec_quadrant_self_interaction_avx(quadrants[idx], &Epot_);
                }
            }

    molec_quadrants_finalize_ghost(quadrants, cellList_parameter, sim);

    *Epot = Epot_;
}
