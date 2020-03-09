#ifndef CDE_H
#define CDE_H

#include "Common.h"

class CDE {
public:
    CDE(const Cosmology& C, const PowerSpectrum& P_l, real epsrel = 1e-3);

    void F1_init(int a, double scalef, double omega0, double par1, double par2, double par3) const;

    double Bcde(int a, double vars[], real k, real p, real x) const;

    void initcde(int a, double scalef, double omega0, double par1, double par2, double par3) const;

    double P22nres(double kmin, double kmax,  int a, real k) const;

    double P13nres(double kmin, double kmax,  int a, real k) const;

    double  PNWloopn(double kmin, double kmax,  double k, int a ) const;

    double Presumn(double kmin, double kmax, double k, int a) const;

    private:
        const Cosmology& C;
        const PowerSpectrum& P_l;
        real epsrel;

    };

    #endif // CDE_H
