#if HAVE_CONFIG_H
# include <config.h>
#endif

#include "ODE.h"
#include "PowerSpectrum.h"
#include "Quadrature.h"
#include "LargeN.h"
#include <functional>

using std::cref;
using std::bind;

const real QMIN = 1e-4;
const real QMAX = 1e1;
const real XMAX = 0.99999;


LargeN::LargeN(const Cosmology& C_, const PowerSpectrum& P, real z0_, real zi_, int Neta_, real epsrel_)
    : C(C_), D(C_), P_L0(P)
{
    z0 = z0_;
    zi = zi_;
    Neta = Neta_;

    etai = log(D(zi)/D(z0));

    epsrel = epsrel_;
}

/* Integrands for Sigma_ab */
static real f_11(const PowerSpectrum& P_L0, real k, real q) {
    real r = q/k;
    real s;
    if(r < 1e-3)
        s = -0.4 + 0.8*pow2(r) - 0.16*pow4(r) - 4/175.*pow6(r);
    else if(fabs(r-1) < 1e-10)
        s = 0.2 + 0.6*(r-1);
    else if(r > 50)
        s = 0.4 - 0.16/pow2(r) - 4/175./pow4(r) - 4/525./pow6(r);
    else
        s = -0.1 + 0.3*r*r - 0.15/r * pow2(r*r-1) * log((1+r)/fabs(1-r));
    return P_L0(q) * s;
}

static real f_12(const PowerSpectrum& P_L0, real k, real q) {
    real r = q/k;
    real s;
    if(r < 1e-3)
        s = -4/15. + 0.8*pow2(r) - 0.8*pow4(r) + 0.16*pow6(r);
    else if(fabs(r-1) < 1e-10)
        s = -1/15. - 0.2*(r-1);
    else if(r > 50)
        s = -8/75. + 4/175./pow2(r) + 4/525./pow4(r) + 4/1155./pow6(r);
    else
        s = -4./15. + 0.5*r*r - 0.3*pow4(r) + 0.15*r*pow2(r*r-1) * log((1+r)/fabs(1-r));
    return P_L0(q) * s;
}

static real f_21(const PowerSpectrum& P_L0, real k, real q) {
    real r = q/k;
    real s;
    if(r < 1e-3)
        s = -0.4 - 0.16*pow2(r) - 4/175.*pow4(r) - 4/525*pow6(r);
    else if(fabs(r-1) < 1e-10)
        s = -0.6 - 0.6*(r-1);
    else if(r > 50)
        s = -1.2 + 0.8/pow2(r) - 0.16/pow4(r) - 4/175./pow6(r);
    else
        s = 0.3/pow2(r) - 0.9 - 0.15/pow3(r) * pow2(r*r-1) * log((1+r)/fabs(1-r));
    return P_L0(q) * s;
}

static real f_22(const PowerSpectrum& P_L0, real k, real q) {
    real r = q/k;
    real s;
    if(r < 1e-3)
        s = -4/15. - 0.8*pow2(r) + 0.16*pow4(r) + 4/175*pow6(r);
    else if(fabs(r-1) < 1e-10)
        s = -13/15. - 0.6*(r-1);
    else if(r > 50)
        s = -16/15. + 0.16/pow2(r) + 4/175./pow4(r) + 4/525./pow6(r);
    else
        s = -17./30. - 0.3*r*r + 0.15/r * pow2(r*r-1) * log((1+r)/fabs(1-r));
    return P_L0(q) * s;
}

real LargeN::Sigma_11(real k) const {
    return k*k/(4*M_PI*M_PI) * Integrate<ExpSub>(bind(f_11, cref(P_L0), k, std::placeholders::_1), QMIN, QMAX, epsrel);
}

real LargeN::Sigma_12(real k) const {
    return k*k/(4*M_PI*M_PI) * Integrate<ExpSub>(bind(f_12, cref(P_L0), k, std::placeholders::_1), QMIN, QMAX, epsrel);
}

real LargeN::Sigma_21(real k) const {
    return k*k/(4*M_PI*M_PI) * Integrate<ExpSub>(bind(f_21, cref(P_L0), k, std::placeholders::_1), QMIN, QMAX, epsrel);
}

real LargeN::Sigma_22(real k) const {
    return k*k/(4*M_PI*M_PI) * Integrate<ExpSub>(bind(f_22, cref(P_L0), k, std::placeholders::_1), QMIN, QMAX, epsrel);
}

/* Inner and outer integrands for Pi_ab mode-coupling terms */
static real g_11(const PowerSpectrum& P_L0, real k, real q, real x) {
    real r = q/k;
    real d = 1 + r*r - 2*r*x;
    return P_L0(q) * P_L0(k*sqrt(d)) * pow2(r + x - 2*r*x*x) / pow2(d);
}

static real g_12(const PowerSpectrum& P_L0, real k, real q, real x) {
    real r = q/k;
    real d = 1 + r*r - 2*r*x;
    return P_L0(q) * P_L0(k*sqrt(d)) * (r + x - 2*r*x*x)*(x - r) / pow2(d);
}

static real g_22(const PowerSpectrum& P_L0, real k, real q, real x) {
    real r = q/k;
    real d = 1 + r*r - 2*r*x;
    return P_L0(q) * P_L0(k*sqrt(d)) * pow2(x - r) / pow2(d);
}

real LargeN::Pi_11(real k) const {
    real a[2] = { QMIN, -1 };
    real b[2] = { QMAX, XMAX };
    return k*k/(8*M_PI*M_PI) * Integrate<2>(bind(g_11, cref(P_L0), k, std::placeholders::_1, std::placeholders::_2), a, b, epsrel);
}

real LargeN::Pi_12(real k) const {
    real a[2] = { QMIN, -1 };
    real b[2] = { QMAX, XMAX };
    return k*k/(8*M_PI*M_PI) * Integrate<2>(bind(g_12, cref(P_L0), k, std::placeholders::_1, std::placeholders::_2), a, b, epsrel);
}

real LargeN::Pi_22(real k) const {
    real a[2] = { QMIN, -1 };
    real b[2] = { QMAX, XMAX };
    return k*k/(8*M_PI*M_PI) * Integrate<2>(bind(g_22, cref(P_L0), k, std::placeholders::_1, std::placeholders::_2), a, b, epsrel);
}

real LargeN::P(real k, real z) const {
    real eta = log(D(z)/D(z0));

    real pi11 = Pi_11(k);
    real pi12 = Pi_12(k);
    real pi22 = Pi_22(k);

    array R_11, R_12, R_21, R_22;
    ComputeR(k, eta, Neta, R_11, R_12, R_21, R_22);

    /* Integrate R_ab(eta,eta1) with respect to eta1 */
    real deta = (eta - etai)/(Neta-1);
    array s11(Neta), s12(Neta), s21(Neta), s22(Neta);
    for(int n = 0; n < Neta; n++) {
        real eta1 = etai + n*deta;
        s11[n] = exp(2*eta1) * R_11[n];
        s12[n] = exp(2*eta1) * R_12[n];
        s21[n] = exp(2*eta1) * R_21[n];
        s22[n] = exp(2*eta1) * R_22[n];
    }
    /* S_{ab} = \int_{\eta_i}^\eta e^{\eta_1} R_{ab}(k;\eta_1,\eta_i) ~ d\eta_1 */
    real S11 = DiscreteIntegrate(Neta, s11, deta);
    real S12 = DiscreteIntegrate(Neta, s12, deta);
//    real S21 = DiscreteIntegrate(Neta, s21, deta);
//    real S22 = DiscreteIntegrate(Neta, s22, deta);

    return exp(2*etai)*pow2(R_11[0] + R_12[0])*P_L0(k) + S11*pi11*S11 + 2*S11*pi12*S12 + S12*pi22*S12;
}

void LargeN::P_ab(real k, real z, real& p11, real& p12, real& p22) const {
    real eta = log(D(z)/D(z0));

    real pi11 = Pi_11(k);
    real pi12 = Pi_12(k);
    real pi22 = Pi_22(k);

    array R_11, R_12, R_21, R_22;
    ComputeR(k, eta, Neta, R_11, R_12, R_21, R_22);

    /* Integrate R_ab(eta,eta1) with respect to eta1 */
    real deta = (eta - etai)/(Neta-1);
    array s11(Neta), s12(Neta), s21(Neta), s22(Neta);
    for(int n = 0; n < Neta; n++) {
        real eta1 = etai + n*deta;
        s11[n] = exp(2*eta1) * R_11[n];
        s12[n] = exp(2*eta1) * R_12[n];
        s21[n] = exp(2*eta1) * R_21[n];
        s22[n] = exp(2*eta1) * R_22[n];
    }
    /* S_{ab} = \int_{\eta_i}^\eta e^{\eta_1} R_{ab}(k;\eta_1,\eta_i) ~ d\eta_1 */
    real S11 = DiscreteIntegrate(Neta, s11, deta);
    real S12 = DiscreteIntegrate(Neta, s12, deta);
    real S21 = DiscreteIntegrate(Neta, s21, deta);
    real S22 = DiscreteIntegrate(Neta, s22, deta);

    p11 = exp(2*etai)*pow2(R_11[0] + R_12[0])*P_L0(k) + S11*pi11*S11 + 2*S11*pi12*S12 + S12*pi22*S12;
    p12 = exp(2*etai)*(R_11[0] + R_12[0])*(R_21[0] + R_22[0])*P_L0(k) + S11*pi11*S21 + S11*pi12*S22 + S12*pi12*S21 + S12*pi22*S22;
    p22 = exp(2*etai)*pow2(R_21[0] + R_22[0])*P_L0(k) + S21*pi11*S21 + 2*S21*pi12*S22 + S22*pi22*S22;
}

/* Derivative function for integrating R_ab according to Eqs. (70), (71) */
static vector<real> dXdeta1(real eta1, const vector<real>& X, real s11, real s12, real s21, real s22) {
    vector<real> dX(6);
    dX[0] = X[2];
    dX[1] = X[3];
    dX[2] = X[4];
    dX[3] = X[5];
    dX[4] = 1.5*X[4] + X[5] + X[2] - 1.5*X[3] - X[1]
           + exp(2*eta1)*((s11 + s22)*X[2] + 2.5*s11*X[0] + 2.5*s12*X[1]);
    dX[5] = 1.5*X[4] + X[5] - 2.25*X[2] + 1.75*X[3] - 1.5*X[0] + 0.5*X[1]
           + exp(2*eta1)*((s11 + s22)*X[3] + 2.5*s21*X[0] + 2.5*s22*X[1]);
    return dX;
}

void LargeN::ComputeR(real k, real eta1, int Neta, array& R_11, array& R_12, array& R_21, array& R_22) const {
    R_11.resize(Neta);
    R_12.resize(Neta);
    R_21.resize(Neta);
    R_22.resize(Neta);

    real s11 = Sigma_11(k);
    real s12 = Sigma_12(k);
    real s21 = Sigma_21(k);
    real s22 = Sigma_22(k);

    real deta = (eta1 - etai)/(Neta-1);
    real eta2;
    vector<real> x0(6);
    vector<vector<real> > X;
    for(int n = 0; n < Neta; n++) {
        eta2 = etai + n*deta;
        x0[0] = 1;
        x0[1] = 0;
        x0[2] = 0;
        x0[3] = 1.5;
        x0[4] = 1.5 + exp(2*eta2)*(s11 + s22);
        x0[5] = -0.75;
        X = RungeKutta4(eta2, eta1, x0, bind(dXdeta1, std::placeholders::_1, std::placeholders::_2, s11, s12, s21, s22), Neta);
        R_11[n] = X[0].back();
        R_21[n] = X[1].back();

        x0[0] = 0;
        x0[1] = 1;
        x0[2] = 1;
        x0[3] = -0.5;
        x0[4] = -0.5;
        x0[5] = 1.75 + exp(2*eta2)*(s11 + s22);
        X = RungeKutta4(eta2, eta1, x0, bind(dXdeta1, std::placeholders::_1, std::placeholders::_2, s11, s12, s21, s22), Neta);
        R_12[n] = X[0].back();
        R_22[n] = X[1].back();
    }
}
