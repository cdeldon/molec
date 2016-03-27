#ifndef MOLEC_INLINE_H
#define MOLEC_INLINE_H

#include <molec/Common.h>

/**
 * Calculate distance between x and y taking periodic boundaries into account
 */
MOLEC_INLINE Real dist(Real x, Real y, Real L)
{
    Real r = x - y;
    if(r < -L / 2)
        r += L;
    else if(r > L / 2)
        r -= L;
    return r;
}


/**
 * Calculate the positive modulo between two integers, used for periodic BC
 */
MOLEC_INLINE int mod(int b, int m)
{
    return (b % m + m) % m;
}
/**
  *Update routine for the innermost force calculation between two particles
  */
MOLEC_INLINE Real update(Real xij,Real yij,Real zij,Real* f_xi,Real* f_yi,Real* f_zi,Real* Epot_)
{
    assert(molec_parameter);
    const Real sigLJ = molec_parameter->sigLJ;
    const Real epsLJ = molec_parameter->epsLJ;
    const Real Rcut2 = molec_parameter->Rcut2;

    const Real r2 = xij * xij + yij * yij + zij * zij;

    Real fr=0.;

    if(r2 < Rcut2)
    {
        // V(s) = 4 * eps * (s^12 - s^6) with  s = sig/r
        const Real s2 = (sigLJ * sigLJ) / r2;
        const Real s6 = s2 * s2 * s2;

        *Epot_ += 4 * epsLJ * (s6 * s6 - s6);

        fr = 24 * epsLJ / r2 * (2 * s6 * s6 - s6);

        *f_xi += fr * xij;
        *f_yi += fr * yij;
        *f_zi += fr * zij;
    }
    return fr;
}



#endif
