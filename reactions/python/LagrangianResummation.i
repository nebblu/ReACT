%{
#include "LagrangianResummation.h"
%}

class LagrangianResummation {
public:
    LagrangianResummation(const Cosmology& C, real z, real epsrel = 1e-4);

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
};
