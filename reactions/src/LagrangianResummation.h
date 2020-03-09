#ifndef LAGRANGIAN_RESUMMATION_H
#define LAGRANGIAN_RESUMMATION_H

#include "Common.h"
#include "LinearPS.h"

/*******************************************************************************
 * LagrangianResummation
 *
 * Matsubara's Lagrangian resummation scheme, where the real-space power
 * spectrum is given to 1-loop order by
 *   P(k) = e^{-Ak^2} [P_L(k) + P22(k) + P13(k)].
 * For details see T. Matsubara, arXiv:0711.2521.
 ******************************************************************************/

class LagrangianResummation {
public:
    LagrangianResummation(const Cosmology& C, real z, real epsrel = 1e-4);
    ~LagrangianResummation();

    /* Exponential damping coefficient */
    real A;

/***** Real space *****/

    /* 1-loop power spectrum in Lagrangian resummation theory */
    real P(real k) const;

    /* Separate mode-coupling terms */
    real P22(real k) const;
    real P13(real k) const;

    /* Nonlinear propagator computed from Lagrangian resummation formalism
     * (normalized to 1 on large scales) */
    real G(real k) const;


/***** Redshift space *****/

    /* 1-loop redshift-space power spectrum */
    real Ps(real k, real mu) const;

    /* Angular moment of the 1-loop redshift-space power spectrum */
    real Ps_ell(real k, int ell) const;

protected:
    const Cosmology& C;
    real z;
    LinearPS P_L;
    real epsrel;
    real f;
};

#endif // LAGRANGIAN_RESUMMATION_H
