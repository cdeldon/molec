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

#include <molec/Parameter.h>
#include <molec/Integrator.h>

void molec_integrator_leapfrog_refrence(
    float* x, float* v, const float* f, float* Ekin, const int N)
{
    assert(molec_parameter);
    const float dt = molec_parameter->dt;
    const float m = molec_parameter->mass;
    const float m0125 = 0.125 * m;

    float v_old = 0;
    float Ekin_ = 0;

    // Integrate velocity: v_{i+1/2} = v_{i-1/2} + dt * f_i / m
    for(int i = 0; i < N; ++i)
    {
        v_old = v[i];
        v[i] = v[i] + dt * f[i] / m;

        // Lineraly interpolate v_i with v_{i-1/2} and v_{i+1/2}
        Ekin_ += m0125 * (v[i] + v_old) * (v[i] + v_old);
    }

    // Integrate position: x_i = x_{i-1} + dt * v_{i-1/2}
    for(int i = 0; i < N; ++i)
        x[i] = x[i] + dt * v[i];

    *Ekin = Ekin_;
}

void molec_integrator_leapfrog_unroll_2(
    float* x, float* v, const float* f, float* Ekin, const int N)
{
    assert(molec_parameter);
    const float dt = molec_parameter->dt;
    const float m = molec_parameter->mass;
    const float m0125 = 0.125 * m;
    const float minv = 1.0 / molec_parameter->mass;

    // Loop logic
    int i = 0;
    const int N2 = N / 2;
    const int N2_upper = N2 * 2;

    // Temporaries
    float x_00, x_01;
    float v_00, v_01;
    float v_2_00, v_2_01;
    float v_old_00, v_old_01;
    float f_minv_00, f_minv_01;

    float Ekin_00 = 0.0, Ekin_01 = 0.0;

    for(i = 0; i < N2_upper; i += 2)
    {
        // Load
        v_old_00 = v[i + 0];
        v_old_01 = v[i + 1];

        v_00 = v[i + 0];
        v_01 = v[i + 1];

        // Compute
        f_minv_00 = f[i + 0] * minv;
        f_minv_01 = f[i + 1] * minv;

        v_00 = v_00 + dt * f_minv_00;
        v_01 = v_01 + dt * f_minv_01;

        v_2_00 = (v_00 + v_old_00) * (v_00 + v_old_00);
        v_2_01 = (v_01 + v_old_01) * (v_01 + v_old_01);

        Ekin_00 = Ekin_00 + m0125 * v_2_00;
        Ekin_01 = Ekin_01 + m0125 * v_2_01;

        // Store
        v[i + 0] = v_00;
        v[i + 1] = v_01;
    }
    for(i = N2_upper; i < N; ++i)
    {
        v_old_00 = v[i];
        v[i] = v[i] + dt * f[i] * minv;
        Ekin_00 = Ekin_00 + m0125 * (v[i] + v_old_00) * (v[i] + v_old_00);
    }

    *Ekin = Ekin_00 + Ekin_01;

    for(i = 0; i < N2_upper; i += 2)
    {
        // Load
        x_00 = x[i + 0];
        x_01 = x[i + 1];

        // Compute
        x_00 = x_00 + dt * v[i + 0];
        x_01 = x_01 + dt * v[i + 1];

        // Store
        x[i + 0] = x_00;
        x[i + 1] = x_01;
    }
    for(i = N2_upper; i < N; ++i)
        x[i] = x[i] + dt * v[i];
}

void molec_integrator_leapfrog_unroll_4(
    float* x, float* v, const float* f, float* Ekin, const int N)
{
    assert(molec_parameter);
    const float dt = molec_parameter->dt;
    const float m = molec_parameter->mass;
    const float m0125 = 0.125 * m;
    const float minv = 1.0 / molec_parameter->mass;

    // Loop logic
    int i = 0;
    const int N4 = N / 4;
    const int N4_upper = N4 * 4;

    // Temporaries
    float x_00, x_01, x_02, x_03;
    float v_00, v_01, v_02, v_03;
    float vi_00, vi_01, vi_02, vi_03;
    float v_2_00, v_2_01, v_2_02, v_2_03;
    float v_old_00, v_old_01, v_old_02, v_old_03;
    float a_00, a_01, a_02, a_03;

    float Ekin_00 = 0.0, Ekin_01 = 0.0, Ekin_02 = 0.0, Ekin_03 = 0.0;

    for(i = 0; i < N4_upper; i += 4)
    {
        // Load
        v_old_00 = v[i + 0];
        v_old_01 = v[i + 1];
        v_old_02 = v[i + 2];
        v_old_03 = v[i + 3];

        v_00 = v[i + 0];
        v_01 = v[i + 1];
        v_02 = v[i + 2];
        v_03 = v[i + 3];

        // Compute
        a_00 = f[i + 0] * minv;
        a_01 = f[i + 1] * minv;
        a_02 = f[i + 2] * minv;
        a_03 = f[i + 3] * minv;

        v_00 = v_00 + dt * a_00;
        v_01 = v_01 + dt * a_01;
        v_02 = v_02 + dt * a_02;
        v_03 = v_03 + dt * a_03;

        vi_00 = v_00 + v_old_00;
        vi_01 = v_01 + v_old_01;
        vi_02 = v_02 + v_old_02;
        vi_03 = v_03 + v_old_03;

        v_2_00 = vi_00 * vi_00;
        v_2_01 = vi_01 * vi_01;
        v_2_02 = vi_02 * vi_02;
        v_2_03 = vi_03 * vi_03;

        Ekin_00 = Ekin_00 + m0125 * v_2_00;
        Ekin_01 = Ekin_01 + m0125 * v_2_01;
        Ekin_02 = Ekin_02 + m0125 * v_2_02;
        Ekin_03 = Ekin_03 + m0125 * v_2_03;

        // Store
        v[i + 0] = v_00;
        v[i + 1] = v_01;
        v[i + 2] = v_02;
        v[i + 3] = v_03;
    }
    for(i = N4_upper; i < N; ++i)
    {
        v_old_00 = v[i];
        v[i] = v[i] + dt * f[i] * minv;
        Ekin_00 = Ekin_00 + m0125 * (v[i] + v_old_00) * (v[i] + v_old_00);
    }

    *Ekin = Ekin_00 + Ekin_01 + Ekin_02 + Ekin_03;

    for(i = 0; i < N4_upper; i += 4)
    {
        // Load
        x_00 = x[i + 0];
        x_01 = x[i + 1];
        x_02 = x[i + 2];
        x_03 = x[i + 3];

        // Compute
        x_00 = x_00 + dt * v[i + 0];
        x_01 = x_01 + dt * v[i + 1];
        x_02 = x_02 + dt * v[i + 2];
        x_03 = x_03 + dt * v[i + 3];

        // Store
        x[i + 0] = x_00;
        x[i + 1] = x_01;
        x[i + 2] = x_02;
        x[i + 3] = x_03;
    }
    for(i = N4_upper; i < N; ++i)
        x[i] = x[i] + dt * v[i];
}
