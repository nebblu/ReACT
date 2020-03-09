%{
#include "LargeN.h"
%}

class LargeN {
public:
    LargeN(const Cosmology& C, const PowerSpectrum& P_L0, real z0, real zi, int Neta, real epsrel = 1e-4);

    real Sigma_11(real k) const;
    real Sigma_12(real k) const;
    real Sigma_21(real k) const;
    real Sigma_22(real k) const;

    real Pi_11(real k) const;
    real Pi_12(real k) const;
    real Pi_22(real k) const;

    /* Compute PS */
    real P(real k, real z) const;
    void P_ab(real k, real z, real& OUTPUT, real& OUTPUT, real& OUTPUT) const;

    /* Compute response function (i.e. the nonlinear propagator) */
    void ComputeR(real k, real eta1, int Neta, array& R_11, array& R_12, array& R_21, array& R_22) const;
};
