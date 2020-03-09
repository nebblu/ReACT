%{
#include "RPT.h"
%}

class RPT {
public:
    RPT(const Cosmology& C, const PowerSpectrum& P_i, real z_i, real z, int Neta, int Nk, real kcut);

    /* Compute nonlinear propagator */
    void ComputeG_ab(real k, real eta1, real eta2, real& OUTPUT, real& OUTPUT, real& OUTPUT, real& OUTPUT) const;

    /* Accessor functions for precomputed propagator */
    real G_11(real k, int t = -1, int tp = 0) const;
    real G_12(real k, int t = -1, int tp = 0) const;
    real G_21(real k, int t = -1, int tp = 0) const;
    real G_22(real k, int t = -1, int tp = 0) const;
    real G_1(real k, int t = -1) const;
    real G_2(real k, int t = -1) const;

    /* Tree-level PS */
    real P1(real k, int a = 1, int b = 1) const;
    real P1_11(real k) const;
    real P1_12(real k) const;
    real P1_22(real k) const;

    /* 1-loop contribution to PS */
    real P2(real k, int a = 1, int b = 1) const;
    real P2_11(real k) const;
    real P2_12(real k) const;
    real P2_22(real k) const;

    /* 2-loop contribution to PS */
    real P3(real k, int a = 1, int b = 1) const;
    real P3_11(real k) const;
    real P3_12(real k) const;
    real P3_22(real k) const;

    /* Linear PS */
    real P_L(real k) const;
};
