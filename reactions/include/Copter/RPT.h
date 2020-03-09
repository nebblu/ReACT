#ifndef RPT_H
#define RPT_H

#include "Common.h"
#include "GrowthFunction.h"
#include "MonteCarlo.h"
#include "PowerSpectrum.h"
#include "array.h"


/******************************************************************************
 * RPT
 *
 * Crocce and Scoccimarro's renormalized perturbation theory
 ******************************************************************************/
class RPT {
public:
    RPT(const Cosmology& C, const PowerSpectrum& P_i, real z_i, real z, int Neta, int Nk = 1000, real kcut = 10.);
    ~RPT();

    /* Nonlinear propagator */
//    real G(int a, int b, real k, real z1, real z2) const;
    void ComputeG_ab(real k, real eta1, real eta2, real& g11, real& g12, real& g21, real& g22) const;

    /* Accessor functions for interpolated propagator */
    real G_11(real k, int t = -1, int tp = 0) const;
    real G_12(real k, int t = -1, int tp = 0) const;
    real G_21(real k, int t = -1, int tp = 0) const;
    real G_22(real k, int t = -1, int tp = 0) const;
    real G_1(real k, int t = -1) const;
    real G_2(real k, int t = -1) const;

    /* Tree-level power spectrum */
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

protected:
    const Cosmology& C;
    const PowerSpectrum& P_i;   // linear PS at z = z_i
    real z_i;                   // initial redshift
    real z;                     // redshift at which to evaluate PS
    GrowthFunction D;

    int Neta;                   // number of time divisions for interpolated propagators
    int Nk;                     // number of k-values for interpolated propagators
    vector<Spline> g11, g12, g21, g22; // nonlinear propagator
    real deta;
    real kcut;

    void PrecomputePropagator();

    /* Renormalized propagator factors */
    real G_alpha(real eta, real etap = 0) const;
    real G_beta_g(real eta, real etap = 0) const;
    real G_beta_d(real eta, real etap = 0) const;
    real G_gamma_g(real eta, real etap = 0) const;
    real G_gamma_d(real eta, real etap = 0) const;
    real G_delta(real eta, real etap = 0) const;
    real G_f(real k) const;
    real G_g(real k) const;
    real G_h(real k) const;
    real G_i(real k) const;

    /* Vertices */
    static real alpha(real k, real q, real r);
    static real beta(real k, real q, real r);

    /* Mode-coupling integrands */
//    real f2_11(real k, real q, real r) const;
//    real f2_12(real k, real q, real x) const;
//    real f2_22(real k, real q, real x) const;
    real F2_11(real k, real x, real y) const;
    real F2_12(real k, real x, real y) const;
    real F2_22(real k, real x, real y) const;

    struct ThirdOrderIntegral : public MonteCarloIntegral {
        ThirdOrderIntegral(const RPT& rpt);
        virtual void Integrand(const double x[], double* f, double* param) const;

        const RPT& rpt;
    } p3int;

    /* Mode-coupling integral factors */
    real I_1(real k, real q, real kq, int t = -1) const;
    real I_2(real k, real q, real kq, int t = -1) const;
    real J_1(real k, real p, real kp, real q, real pq, int t = -1) const;
    real J_2(real k, real p, real kp, real q, real pq, int t = -1) const;
};

#endif // RPT_H
