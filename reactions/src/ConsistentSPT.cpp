#if HAVE_CONFIG_H
# include <config.h>
#endif


#include "ConsistentSPT.h"
#include "ODE.h"
#include "Quadrature.h"
#include <functional>

using std::cref;
using std::bind;


/* Limits of integration for second-order power spectrum */
const real QMIN = 1e-5;
const real QMAX = 1e2;
const real XMAX = 0.99999;


ConsistentSPT::ConsistentSPT(const Cosmology& C_, const PowerSpectrum& P0_, real z0, real z, real epsrel_)
    : C(C_), P0(P0_)
{
    a0 = 1/(1 + z0);
    zdefault = z;
    epsrel = epsrel_;

    r = C.Omega_Lambda/C.Omega_m;
    c0 = (sqrt(1 + r*pow3(a0)) - 1)/(sqrt(1 + r*pow3(a0)) + 1);
    A = P0.VelocityDispersion();

    /* Pre-compute time dependent functions */
    real s1 = sigma(1.0001);
    vector<real> S(Nsteps+1);
    for(int i = 0; i <= Nsteps; i++)
        S[i] = i*s1/Nsteps;

    vector<real> f0(18);
    f0[0] = 1;
    for(int i = 1; i < 9; i++)
        f0[i] = 0;
    f0[9] = -1;
    for(int i = 10; i < 18; i++)
        f0[i] = 0;

    vector<vector<real> > F = RungeKutta4(0, s1, f0, bind(&ConsistentSPT::dfds, this, std::placeholders::_1, std::placeholders::_2), Nsteps);
    U    = CubicSpline(S, F[0]);
    VA   = CubicSpline(S, F[1]);
    VB   = CubicSpline(S, F[2]);
    WA0A = CubicSpline(S, F[3]);
    WAA0 = CubicSpline(S, F[4]);
    WA0B = CubicSpline(S, F[5]);
    WAB0 = CubicSpline(S, F[6]);
    WB0A = CubicSpline(S, F[7]);
    WB0B = CubicSpline(S, F[8]);
    X    = CubicSpline(S, F[9]);
    YA   = CubicSpline(S, F[10]);
    YB   = CubicSpline(S, F[11]);
    ZA0A = CubicSpline(S, F[12]);
    ZAA0 = CubicSpline(S, F[13]);
    ZA0B = CubicSpline(S, F[14]);
    ZAB0 = CubicSpline(S, F[15]);
    ZB0A = CubicSpline(S, F[16]);
    ZB0B = CubicSpline(S, F[17]);
}

real ConsistentSPT::P(real k, int a, int b) const {
    switch(a*b) {
        case 1:
            return P_dd(k);
        case 2:
            return P_dt(k);
        case 4:
            return P_tt(k);
        default:
            warning("ConsistentSPT: invalid indices, a = %d, b = %d\n", a, b);
            return 0;
    }
}

real ConsistentSPT::P_dd(real k, real z) const {
    if(z == -1)
        z = zdefault;
    real a = 1/(1+z);
    real s = sigma(a);
    real u = U(s);
    real vA = VA(s);
    real vB = VB(s);
    real wA0A = WA0A(s);
    real wAA0 = WAA0(s);
    real wA0B = WA0B(s);
    real wAB0 = WAB0(s);
    real wB0A = WB0A(s);
    real wB0B = WB0B(s);
    real pAA = P_AA(k);
    real pAB = P_AB(k);
    real pBB = P_BB(k);
    real pA0A = P_A0A(k);
    real pAA0 = P_AA0(k);
    real pA0B = P_A0B(k);
    real pAB0 = P_AB0(k);
    real pB0A = P_B0A(k);
    real pB0B = P_B0B(k);
    return exp(2*s) *
        ( u*u*P0(k) + vA*vA*pAA + 2*vA*vB*pAB + vB*vB*pBB
          + 2*u*(wA0A*pA0A + wAA0*pAA0 + wA0B*pA0B + wAB0*pAB0 + 2*wB0A*pB0A + 2*wB0B*pB0B) );
}

real ConsistentSPT::P_dt(real k, real z) const {
    if(z == -1)
        z = zdefault;
    real a = 1/(1+z);
    real s = sigma(a);

    real u = U(s);
    real vA = VA(s);
    real vB = VB(s);
    real wA0A = WA0A(s);
    real wAA0 = WAA0(s);
    real wA0B = WA0B(s);
    real wAB0 = WAB0(s);
    real wB0A = WB0A(s);
    real wB0B = WB0B(s);
    real x = X(s);
    real yA = YA(s);
    real yB = YB(s);
    real zA0A = ZA0A(s);
    real zAA0 = ZAA0(s);
    real zA0B = ZA0B(s);
    real zAB0 = ZAB0(s);
    real zB0A = ZB0A(s);
    real zB0B = ZB0B(s);
    real p0 = P0(k);
    real pAA = P_AA(k);
    real pAB = P_AB(k);
    real pBB = P_BB(k);
    real pA0A = P_A0A(k);
    real pAA0 = P_AA0(k);
    real pA0B = P_A0B(k);
    real pAB0 = P_AB0(k);
    real pB0A = P_B0A(k);
    real pB0B = P_B0B(k);

    return u/x * exp(2*s) *
        ( u*x*p0 + vA*yA*pAA + (vA*yB+vB*yA)*pAB + vB*yB*pBB
        + (u*zA0A+x*wA0A)*pA0A + (u*zAA0+x*wAA0)*pAA0
        + (u*zA0B+x*wA0B)*pA0B + (u*zAB0+x*wAB0)*pAB0
        + 2*(u*zB0A+x*wB0A)*pB0A + 2*(u*zB0B+x*wB0B)*pB0B );
}

real ConsistentSPT::P_tt(real k, real z) const {
    if(z == -1)
        z = zdefault;
    real a = 1/(1+z);
    real s = sigma(a);

    real u = U(s);
    real x = X(s);
    real yA = YA(s);
    real yB = YB(s);
    real zA0A = ZA0A(s);
    real zAA0 = ZAA0(s);
    real zA0B = ZA0B(s);
    real zAB0 = ZAB0(s);
    real zB0A = ZB0A(s);
    real zB0B = ZB0B(s);
    real p0 = P0(k);
    real pAA = P_AA(k);
    real pAB = P_AB(k);
    real pBB = P_BB(k);
    real pA0A = P_A0A(k);
    real pAA0 = P_AA0(k);
    real pA0B = P_A0B(k);
    real pAB0 = P_AB0(k);
    real pB0A = P_B0A(k);
    real pB0B = P_B0B(k);

    return pow2(u/x) * exp(2*s) *
        ( x*x*p0 + yA*yA*pAA + 2*yA*yB*pAB + yB*yB*pBB
        + 2*x*(zA0A*pA0A + zAA0*pAA0 + zA0B*pA0B + zAB0*pAB0 + 2*zB0A*pB0A + 2*zB0B*pB0B) );
}

real ConsistentSPT::P22(real k, int a, int b) const {
    switch(a*b) {
        case 1:
            return P22_dd(k);
        case 2:
            return P22_dt(k);
        case 4:
            return P22_tt(k);
        default:
            warning("ConsistentSPT: invalid indices, a = %d, b = %d\n", a, b);
            return 0;
    }
}

real ConsistentSPT::P13(real k, int a, int b) const {
    switch(a*b) {
        case 1:
            return P13_dd(k);
        case 2:
            return P13_dt(k);
        case 4:
            return P13_tt(k);
        default:
            warning("ConsistentSPT: invalid indices, a = %d, b = %d\n", a, b);
            return 0;
    }
}

real ConsistentSPT::P22_dd(real k, real z) const {
    if(z == -1)
        z = zdefault;
    real a = 1/(1+z);
    real s = sigma(a);

    real vA = VA(s);
    real vB = VB(s);
    real pAA = P_AA(k);
    real pAB = P_AB(k);
    real pBB = P_BB(k);

    return exp(2*s) * (vA*vA*pAA + 2*vA*vB*pAB + vB*vB*pBB);
}

real ConsistentSPT::P22_dt(real k, real z) const {
    if(z == -1)
        z = zdefault;
    real a = 1/(1+z);
    real s = sigma(a);

    real u = U(s);
    real vA = VA(s);
    real vB = VB(s);
    real yA = YA(s);
    real yB = YB(s);
    real pAA = P_AA(k);
    real pAB = P_AB(k);
    real pBB = P_BB(k);

    return pow2(u) * exp(2*s) * (vA*yA*pAA + (vA*yB+vB*yA)*pAB + vB*yB*pBB);
}

real ConsistentSPT::P22_tt(real k, real z) const {
    if(z == -1)
        z = zdefault;
    real a = 1/(1+z);
    real s = sigma(a);

    real u = U(s);
    real yA = YA(s);
    real yB = YB(s);
    real pAA = P_AA(k);
    real pAB = P_AB(k);
    real pBB = P_BB(k);

    return pow2(u) * exp(2*s) * (yA*yA*pAA + 2*yA*yB*pAB + yB*yB*pBB);
}

real ConsistentSPT::P13_dd(real k, real z) const {
    if(z == -1)
        z = zdefault;
    real a = 1/(1+z);
    real s = sigma(a);
    real u = U(s);
    real wA0A = WA0A(s);
    real wAA0 = WAA0(s);
    real wA0B = WA0B(s);
    real wAB0 = WAB0(s);
    real wB0A = WB0A(s);
    real wB0B = WB0B(s);
    real pA0A = P_A0A(k);
    real pAA0 = P_AA0(k);
    real pA0B = P_A0B(k);
    real pAB0 = P_AB0(k);
    real pB0A = P_B0A(k);
    real pB0B = P_B0B(k);
    return exp(2*s) * 2*u*(wA0A*pA0A + wAA0*pAA0 + wA0B*pA0B + wAB0*pAB0 + 2*wB0A*pB0A + 2*wB0B*pB0B);
}

real ConsistentSPT::P13_dt(real k, real z) const {
    if(z == -1)
        z = zdefault;
    real a = 1/(1+z);
    real s = sigma(a);

    real u = U(s);
    real x = X(s);
    real wA0A = WA0A(s);
    real wAA0 = WAA0(s);
    real wA0B = WA0B(s);
    real wAB0 = WAB0(s);
    real wB0A = WB0A(s);
    real wB0B = WB0B(s);
    real zA0A = ZA0A(s);
    real zAA0 = ZAA0(s);
    real zA0B = ZA0B(s);
    real zAB0 = ZAB0(s);
    real zB0A = ZB0A(s);
    real zB0B = ZB0B(s);
    real pA0A = P_A0A(k);
    real pAA0 = P_AA0(k);
    real pA0B = P_A0B(k);
    real pAB0 = P_AB0(k);
    real pB0A = P_B0A(k);
    real pB0B = P_B0B(k);

    return u/x * exp(2*s) *
        ( (u*zA0A+x*wA0A)*pA0A + (u*zAA0+x*wAA0)*pAA0
        + (u*zA0B+x*wA0B)*pA0B + (u*zAB0+x*wAB0)*pAB0
        + 2*(u*zB0A+x*wB0A)*pB0A + 2*(u*zB0B+x*wB0B)*pB0B );
}

real ConsistentSPT::P13_tt(real k, real z) const {
    if(z == -1)
        z = zdefault;
    real a = 1/(1+z);
    real s = sigma(a);

    real u = U(s);
    real x = X(s);
    real zA0A = ZA0A(s);
    real zAA0 = ZAA0(s);
    real zA0B = ZA0B(s);
    real zAB0 = ZAB0(s);
    real zB0A = ZB0A(s);
    real zB0B = ZB0B(s);
//    real p0 = P0(k);
//    real pAA = P_AA(k);
//    real pAB = P_AB(k);
//    real pBB = P_BB(k);
    real pA0A = P_A0A(k);
    real pAA0 = P_AA0(k);
    real pA0B = P_A0B(k);
    real pAB0 = P_AB0(k);
    real pB0A = P_B0A(k);
    real pB0B = P_B0B(k);

    return pow2(u/x) * exp(2*s) * 2*x*(zA0A*pA0A + zAA0*pAA0 + zA0B*pA0B + zAB0*pAB0 + 2*zB0A*pB0A + 2*zB0B*pB0B);
}

void ConsistentSPT::ComputeAll(real k, real& p11, real& p22dd, real& p13dd, real& p22tt, real& p13tt, real& p22dt, real& p13dt, real z) const {
    if(z == -1)
        z = zdefault;
    real a = 1/(1+z);
    real s = sigma(a);

    real u = U(s);
    real vA = VA(s);
    real vB = VB(s);
    real wA0A = WA0A(s);
    real wAA0 = WAA0(s);
    real wA0B = WA0B(s);
    real wAB0 = WAB0(s);
    real wB0A = WB0A(s);
    real wB0B = WB0B(s);
    real x = X(s);
    real yA = YA(s);
    real yB = YB(s);
    real zA0A = ZA0A(s);
    real zAA0 = ZAA0(s);
    real zA0B = ZA0B(s);
    real zAB0 = ZAB0(s);
    real zB0A = ZB0A(s);
    real zB0B = ZB0B(s);
    real p0 = P0(k);
    real pAA = P_AA(k);
    real pAB = P_AB(k);
    real pBB = P_BB(k);
    real pA0A = P_A0A(k);
    real pAA0 = P_AA0(k);
    real pA0B = P_A0B(k);
    real pAB0 = P_AB0(k);
    real pB0A = P_B0A(k);
    real pB0B = P_B0B(k);

    p11 = exp(2*s) * u*u*p0;
    p22dd = exp(2*s) * (vA*vA*pAA + 2*vA*vB*pAB + vB*vB*pBB);
    p13dd = exp(2*s) * 2*u*(wA0A*pA0A + wAA0*pAA0 + wA0B*pA0B + wAB0*pAB0 + 2*wB0A*pB0A + 2*wB0B*pB0B);
    p22tt = pow2(u/x) * exp(2*s) * (yA*yA*pAA + 2*yA*yB*pAB + yB*yB*pBB);
    p13tt = pow2(u/x) * exp(2*s) * 2*x*(zA0A*pA0A + zAA0*pAA0 + zA0B*pA0B + zAB0*pAB0 + 2*zB0A*pB0A + 2*zB0B*pB0B);
    p22dt = u/x * exp(2*s) * (vA*yA*pAA + (vA*yB + vB*yA)*pAB + vB*yB*pBB);
    p13dt = u/x * exp(2*s) * ((u*zA0A+x*wA0A)*pA0A + (u*zAA0+x*wAA0)*pAA0 + (u*zA0B+x*wA0B)*pA0B + (u*zAB0+x*wAB0)*pAB0 + 2*(u*zB0A+x*wB0A)*pB0A + 2*(u*zB0B+x*wB0B)*pB0B);
}

/* Integrands for P_22 mode-coupling terms */
static real f_AA(const PowerSpectrum& P0, real k, real q, real x) {
    real r = q/k;
    real d = 1 + r*r - 2*r*x;
    return P0(q) * P0(k*sqrt(d)) * pow2(r + x - 2*r*x*x) / pow2(d);
}

static real f_AB(const PowerSpectrum& P0, real k, real q, real x) {
    real r = q/k;
    real d = 1 + r*r - 2*r*x;
    return P0(q) * P0(k*sqrt(d)) * (r + x - 2*r*x*x)*(x - r) / pow2(d);
}

static real f_BB(const PowerSpectrum& P0, real k, real q, real x) {
    real r = q/k;
    real d = 1 + r*r - 2*r*x;
    return P0(q) * P0(k*sqrt(d)) * pow2(x - r) / pow2(d);
}

/* Integrands for P_13 terms */
static real c_AA0(const PowerSpectrum& P0, real k, real q) {
    real r = q/k;
    real s;
    if(r < 1e-3)
        s = (8./3.)*pow2(r) - (8./5.)*pow4(r) + (8./35.)*pow6(r);
    else if(fabs(r-1) < 1e-10)
        s = 4./3.;
    else if(r > 100)
        s = 8./5. - (8./35.)/pow2(r) - (8./315.)/pow4(r) - (8./1155.)/pow6(r);
    else
        s = 0.5 + (4./3.)*pow2(r) - 0.5*pow4(r) + 1/(4*r) * pow3(r*r - 1) * log((1+r)/fabs(1-r));

    return P0(q) * s;
}

static real c_B0A(const PowerSpectrum& P0, real k, real q) {
    real r = q/k;
    real s;
    if(r < 1e-3)
        s = -(1./3.) - (4./5.)*pow2(r) + (4./35.)*pow4(r) + (4./315.)*pow6(r);
    else if(fabs(r-1) < 1e-10)
        s = -2;
    else if(r > 100)
        s = -(5./3.) + (4./5.)/pow2(r) - (4./35.)/pow4(r) - (4./315.)/pow6(r);
    else
        s = 0.25/pow2(r) - 1 - 0.25*pow2(r) + 1/(8*pow3(r)) * pow3(r*r - 1) * log((1+r)/fabs(1-r));

    return P0(q) * s;
}

real ConsistentSPT::P_AA(real k) const {
    real a[2] = { QMIN, -1 };
    real b[2] = { QMAX, XMAX };
    return k*k/(8*M_PI*M_PI) * Integrate<2>(bind(f_AA, cref(P0), k, std::placeholders::_1, std::placeholders::_2), a, b, epsrel);
}
real ConsistentSPT::P_AB(real k) const {
    real a[2] = { QMIN, -1 };
    real b[2] = { QMAX, XMAX };
    return k*k/(8*M_PI*M_PI) * Integrate<2>(bind(f_AB, cref(P0), k, std::placeholders::_1, std::placeholders::_2), a, b, epsrel);
}
real ConsistentSPT::P_BB(real k) const {
    real a[2] = { QMIN, -1 };
    real b[2] = { QMAX, XMAX };
    return k*k/(8*M_PI*M_PI) * Integrate<2>(bind(f_BB, cref(P0), k, std::placeholders::_1, std::placeholders::_2), a, b, epsrel);
}
real ConsistentSPT::P_A0A(real k) const {
    return -A*k*k*P0(k);
}
real ConsistentSPT::P_AA0(real k) const {
    return k*k/(4*M_PI*M_PI) * P0(k) * Integrate<ExpSub>(bind(c_AA0, cref(P0), k, std::placeholders::_1), QMIN, QMAX, epsrel);
}
real ConsistentSPT::P_AB0(real k) const {
    return 0;
}
real ConsistentSPT::P_A0B(real k) const {
    return -A*k*k*P0(k);
}
real ConsistentSPT::P_B0A(real k) const {
    return k*k/(4*M_PI*M_PI) * P0(k) * Integrate<ExpSub>(bind(c_B0A, cref(P0), k, std::placeholders::_1), QMIN, QMAX, epsrel);
}
real ConsistentSPT::P_B0B(real k) const {
    return -0.5*A*k*k*P0(k);
}

real ConsistentSPT::sigma(real a) const {
    if(r == 0)
        return log(a/a0);
    else {
        real b = sqrt(1 + r*pow3(a));
        return 1./3. * (log((b - 1)/(b + 1)) - log(c0));
    }
}

real ConsistentSPT::epsilon(real s) const {
    return c0/(exp(-3*s) - c0);
}

vector<real> ConsistentSPT::dfds(real s, const vector<real>& f) {
    real eps = epsilon(s);
    vector<real> D(18);
    D[0] = -f[0] - f[9];                                // U' = -U - X
    D[9] = -1.5*f[0] - (1.5 + eps)*f[9];                // X' = -3/2 U - (3/2 + eps) X
    D[1] = -f[1] - f[10] + exp(s)*f[9]*f[0];            // VA' = -VA - YA + e^s X U
    D[10] = -1.5*f[1] - (1.5 + eps)*f[10];              // YA' = -3/2 VA - (3/2 + eps) YA
    D[2] = -f[2] - f[11];                               // VB' = -VB - YB
    D[11] = -1.5*f[2] - (1.5 + eps)*f[11] + exp(s)*f[9]*f[9];   // YB' = -3/2 VB - (3/2 + eps) YB + e^s X^2
    D[3] = -f[3] - f[12] + exp(s)*f[9]*f[1];
    D[12] = -1.5*f[3] - (1.5 + eps)*f[12];
    D[4] = -f[4] - f[13] + exp(s)*f[10]*f[0];
    D[13] = -1.5*f[4] - (1.5 + eps)*f[13];
    D[5] = -f[5] - f[14] + exp(s)*f[9]*f[2];
    D[14] = -1.5*f[5] - (1.5 + eps)*f[14];
    D[6] = -f[6] - f[15] + exp(s)*f[11]*f[0];
    D[15] = -1.5*f[6] - (1.5 + eps)*f[15];
    D[7] = -f[7] - f[16];
    D[16] = -1.5*f[7] - (1.5 + eps)*f[16] + exp(s)*f[9]*f[10];
    D[8] = -f[8] - f[17];
    D[17] = -1.5*f[8] - (1.5 + eps)*f[17] + exp(s)*f[9]*f[11];
    return D;
}
