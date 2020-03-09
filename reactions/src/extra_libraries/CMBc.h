#ifndef CMBc_H
#define CMBc_H


#include "Common.h"

class CMBc {
public:
    CMBc(const Cosmology& C, const PowerSpectrum& P_L, real epsrel = 1e-3);

    // analytic forms

    // Initialise z(xi) and xi_star (a==2 we initialise non-linear scale for HALOFIT)
    void comdist_init(int a, double omega0, double mg1);

    // BLSS prediction
    double BLSS(int a, int b, double params[], double l1, double l2, double x) const;
    // BLSS lensing kernel
    double BLSS_kernel(int a, int b, double params[], double l1, double l2, double x, double xi) const;

    // Post Born prediction
    double BPB(int a, double omega0, double l1, double l2, double x) const;
    // CMB lensing prediction (sum of BLSSa and BPB with wigner prefactor )
    double CMBla(int a, int b,  double params[],double l1, double l2, double x) const;

private:
    const Cosmology& C;
    const PowerSpectrum& P_L;
    real epsrel;

};

#endif // CMBc_H
