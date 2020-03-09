%{
#include "SPT.h"
%}

class SPT {
public:
    SPT(const Cosmology& C, const PowerSpectrum& P_L, real epsrel = 1e-5);

    /* Full 1-loop power spectrum */
    real P(real k, int a = 1, int b = 1) const;

    /* Separate 1-loop contributions to power spectrum */
    real P13(real k, int a = 1, int b = 1) const;
    real P22(real k, int a = 1, int b = 1) const;

    real P22_dd(real k) const;
    real P22_dt(real k) const;
    real P22_tt(real k) const;
    real P13_dd(real k) const;
    real P13_dt(real k) const;
    real P13_tt(real k) const;

    /* 1-loop propagator (normalized to 1 on large scales) */
    real G(real k) const;
};
