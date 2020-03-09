#ifndef BSPT_H
#define BSPT_H

#include "Common.h"

// splines needed for n_effective and linear power spectrum
extern Spline neffehu;
extern Spline neffnw;
extern Spline mylinearps;

class BSPT {
public:
    BSPT(const Cosmology& C, const PowerSpectrum& P_L, real epsrel = 1e-3);

    // tree level analytic (LCDM + DGP)
    real Btree(int a, double k1, double k2, real x) const;
    // 1-loop analytic (LCDM and DGP)
    real Bloop(int a, double k1, double k2, double x) const;
    // 1-loop terms except for B123
    real Bloopterms(int a, real k1, real k2, real k3, real x) const;


    // tree level numerical
    real Btreen(double vars[], double k1, double k2, double x) const;
    // 1-loop numerical
    real Bloopn(double vars[], double k1, double k2, double x) const;



    // Gil-marin/Namikawa fitting formula
    real Bfit(double k1, double k2, double x) const;
    // Initialize the NW power spectrum, A^{nw,0l}, sigma8 and knl - used in fitting formula
    void mypspline(double scalef, double omega0, int j) const;

private:
    const Cosmology& C;
    const PowerSpectrum& P_L;
    real epsrel;

};

#endif // BSPT_H
