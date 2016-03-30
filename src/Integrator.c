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

void molec_integrator_leapfrog_refrence(Real* x, Real* v, const Real* f, Real* Ekin, const int N)
{
    assert(molec_parameter);
    const Real dt = molec_parameter->dt;
    const Real m = molec_parameter->mass;
    
    Real v_old = 0;
    Real Ekin_ = 0;

    // Integrate velocity: v_{i+1/2} = v_{i-1/2} + dt * f_i / m
    for(int i = 0; i < N; ++i)
    {
        v_old = v[i];
        v[i] = v[i] + dt * f[i] / m;

        // Lineraly interpolate v_i with v_{i-1/2} and v_{i+1/2}
        Ekin_ += 0.125 * m * (v[i] + v_old) * (v[i] + v_old);
    }
    
    // Integrate position: x_i = x_{i-1} + dt * v_{i-1/2}
    for(int i = 0; i < N; ++i)
        x[i] = x[i] + dt * v[i];
    
    *Ekin = Ekin_;
}

