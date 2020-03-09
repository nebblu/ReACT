#ifndef LARGEN_H
#define LARGEN_H

#include "array.h"
#include "Common.h"
#include "GrowthFunction.h"


class LargeN {
public:
    /* The time coordinate \eta is defined so that \eta=0 corresponds to redshift z0.
     * So z0 should be the redshift at which you want to calculate the power spectrum,
     * and P_L0 should be the linear power spectrum at that redshift. */
    LargeN(const Cosmology& C, const PowerSpectrum& P_L0, real z0, real zi, int Neta, real epsrel = 1e-4);

    /* Eq. (60): $\Sigma^+_{0;i_1i_2}(k)$ [use Eq. (61) to get $\Sigma^-$] */
    real Sigma_11(real k) const;
    real Sigma_12(real k) const;
    real Sigma_21(real k) const;
    real Sigma_22(real k) const;

    /* Eq. (66): $\Pi_{0;i_1i_2}(k)$ [note $\Pi_{0;21} = \Pi_{0;12}$] */
    real Pi_11(real k) const;
    real Pi_12(real k) const;
    real Pi_22(real k) const;

    real P(real k, real z) const;
    void P_ab(real k, real z, real& p11, real& p12, real& p22) const;

    /* Computes the response function R_ab(eta1,eta2) (a.k.a. the nonlinear
     * propagator) for Neta discrete values eta2 between etai and eta1 */
    void ComputeR(real k, real eta1, int Neta, array& R_11, array& R_12, array& R_21, array& R_22) const;

protected:
    const Cosmology& C;
    GrowthFunction D;
    const PowerSpectrum& P_L0;
    real z0;
    real zi;
    real etai;
    int Neta;

    real epsrel;
};

#endif // LARGEN_H
