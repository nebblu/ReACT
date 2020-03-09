#if HAVE_CONFIG_H
# include <config.h>
#endif

#include "MonteCarlo.h"
#include "Quadrature.h"
#include "RPT.h"
#include <functional>
using std::cref;
using std::bind;

/* Limits and tolerance for scale-dependent factors (f,g,h,i) in nonlinear propagator */
const real QMIN = 1e-4;
const real QMAX = 1e1;
const real EPSREL = 1e-5;


RPT::RPT(const Cosmology& C_, const PowerSpectrum& P_i_, real z_i_, real z_, int Neta_, int Nk_, real kcut_)
    : C(C_), P_i(P_i_), p3int(*this)
{
    z_i = z_i_;
    z = z_;
    Neta = Neta_;
    Nk = Nk_;
    kcut = kcut_;

    D = GrowthFunction(C);
    deta = log(D(z)/D(z_i))/(Neta-1);

    PrecomputePropagator();

    MonteCarloIntegral::Range domain[5] = {
        { log(QMIN), log(kcut) },
        { log(QMIN), log(kcut) },
        { -1, 1 },
        { -1, 1 },
        { 0, 2*M_PI }
    };
    p3int.SetDomain(domain);
}

RPT::~RPT() {
}


/***************************
 * Renormalized propagator *
 ***************************/

/* Time-dependent factors */
real RPT::G_alpha(real eta, real etap) const {
    return exp(2*etap) * (exp(2*(eta-etap)) - 1.4*exp(eta-etap) + 0.4*exp(-1.5*(eta-etap)));
}
real RPT::G_beta_g(real eta, real etap) const {
    return exp(2*etap) * (0.6*exp(eta-etap) - 1 + 0.4*exp(-1.5*(eta-etap)));
}
real RPT::G_beta_d(real eta, real etap) const {
    return exp(2*etap) * (0.6*exp(3.5*(eta-etap)) - exp(2.5*(eta-etap)) + 0.4*exp(eta-etap));
}
real RPT::G_gamma_g(real eta, real etap) const {
    return exp(2*etap) * (0.4*exp(eta-etap) - exp(-0.5*(eta-etap)) + 0.6*exp(-1.5*(eta-etap)));
}
real RPT::G_gamma_d(real eta, real etap) const {
    return exp(2*etap) * (0.4*exp(3.5*(eta-etap)) - exp(2*(eta-etap)) + 0.6*exp(eta-etap));
}
real RPT::G_delta(real eta, real etap) const {
    return exp(2*etap) * (0.4*exp(3.5*(eta-etap)) - 1.4*exp(eta-etap) + 1);
}

/* k-dependent integrands */
static real B_f(const PowerSpectrum& P, real k, real q) {
    real r = q/k;
    real r2 = pow2(r), r4 = pow4(r);
    if(r < 1e-4)
        return P(q) * (-84 + 464./5.*r2 - 2256./35.*r4);
    else if(fabs(1 - r) < 1e-8)
        return P(q) * (-44 - 4*(1-r));
    else if(r > 10)
        return P(q) * (-244./5. + 48./5/r2 - 80./21./r4);
    else
        return P(q) * (6/r2 - 79 + 50*r2 - 21*r4 + 0.75/pow3(r)*pow3(1-r2)*(2+7*r2)*log(pow2((1-r)/(1+r))));
}
static real B_g(const PowerSpectrum& P, real k, real q) {
    real r = q/k;
    real r2 = pow2(r), r4 = pow4(r);
    if(r < 1e-4)
        return P(q) * (-28 - 16./5.*r2 - 48./7.*r4);
    else if(fabs(1 - r) < 1e-8)
        return P(q) * (-36 + 20*(1-r));
    else if(r > 10)
        return P(q) * (-252./5. + 624./35/r2 - 304./105./r4);
    else
        return P(q) * (6/r2 - 41 + 2*r2 - 3*r4 + 0.75/pow3(r)*pow3(1-r2)*(2+r2)*log(pow2((1-r)/(1+r))));
}
static real B_h(const PowerSpectrum& P, real k, real q) {
    real r = q/k;
    real r2 = pow2(r), r4 = pow4(r);
    if(r < 1e-4)
        return P(q) * (-4 + 64./5.*r2 + 544./35.*r4);
    else if(fabs(1 - r) < 1e-8)
        return P(q) * (16 - 24*(1-r));
    else if(r > 10)
        return P(q) * (76./5. + 256./35/r2 - 32./7./r4);
    else
        return P(q) * (6/r2 + 1 + 9*r4 + 0.75/pow3(r)*pow2(1-r2)*(2+5*r2+3*r4)*log(pow2((1-r)/(1+r))));
}
static real B_i(const PowerSpectrum& P, real k, real q) {
    real r = q/k;
    real r2 = pow2(r), r4 = pow4(r);
    if(r < 1e-4)
        return P(q) * (12 - 16./5.*r2 + 400./7.*r4);
    else if(fabs(1 - r) < 1e-8)
        return P(q) * (44 - 60*(1-r));
    else if(r > 10)
        return P(q) * (268./5. - 16./35/r2 - 208./35./r4);
    else
        return P(q) * (6/r2 + 29 - 18*r2 + 27*r4 + 0.75/pow3(r)*pow2(1-r2)*(2+9*r2+9*r4)*log(pow2((1-r)/(1+r))));
}

/* k-dependent factors */
real RPT::G_f(real k) const {
    return k*k/(4*M_PI*M_PI) * 1/252. * Integrate<ExpSub>(bind(B_f, cref(P_i), k, std::placeholders::_1), QMIN, QMAX, 1e-5);
}
real RPT::G_g(real k) const {
    return k*k/(4*M_PI*M_PI) * 1/84. * Integrate<ExpSub>(bind(B_g, cref(P_i), k, std::placeholders::_1), QMIN, QMAX, 1e-5);
}
real RPT::G_h(real k) const {
    return k*k/(4*M_PI*M_PI) * 1/12. * Integrate<ExpSub>(bind(B_h, cref(P_i), k, std::placeholders::_1), QMIN, QMAX, 1e-5);
}
real RPT::G_i(real k) const {
    return -k*k/(4*M_PI*M_PI) * 1/36. * Integrate<ExpSub>(bind(B_i, cref(P_i), k, std::placeholders::_1), QMIN, QMAX, 1e-5);
}

void RPT::ComputeG_ab(real k, real eta1, real eta2, real& g11, real& g12, real& g21, real& g22) const {
    real eta = eta1 - eta2;
    real a = G_alpha(eta1, eta2);
    real bg = G_beta_g(eta1, eta2);
    real bd = G_beta_d(eta1, eta2);
    real gg = G_gamma_g(eta1, eta2);
    real gd = G_gamma_d(eta1, eta2);
    real d = G_delta(eta1, eta2);
    real f = G_f(k);
    real g = G_g(k);
    real h = G_h(k);
    real i = G_i(k);

    g11 = 0.6*exp(eta)*exp(a*f - bg*i) + 0.4*exp(-1.5*eta)*exp(d*g - gd*h);
    g12 = 0.4*exp(eta)*exp(a*f - bg*h) - 0.4*exp(-1.5*eta)*exp(d*f - gd*h);
    g21 = 0.6*exp(eta)*exp(a*g + gg*h) - 0.6*exp(-1.5*eta)*exp(d*g + bd*i);
    g22 = 0.4*exp(eta)*exp(a*g - 3/2.*gg*i) + 0.6*exp(-1.5*eta)*exp(d*f - 2/3.*bd*h);
}

void RPT::PrecomputePropagator() {
    g11.resize(Neta*Neta);
    g12.resize(Neta*Neta);
    g21.resize(Neta*Neta);
    g22.resize(Neta*Neta);

    /* Compute scale-dependent functions */
    array k(Nk), f(Nk), g(Nk), h(Nk), i(Nk);
    int n;
    #pragma omp parallel default(shared) private(n)
    #pragma omp for schedule(dynamic)
    for(n = 0; n < Nk; n++) {
        k[n] = n*kcut/(Nk-1);
        f[n] = G_f(k[n]);
        g[n] = G_g(k[n]);
        h[n] = G_h(k[n]);
        i[n] = G_i(k[n]);
    }

    /* Spline the propagator at each point in time */
    array g11tmp(Nk), g12tmp(Nk), g21tmp(Nk), g22tmp(Nk);
    for(int t = 0; t < Neta; t++) {
        real eta = t*deta;
        for(int tp = 0; tp < Neta; tp++) {
            real etap = tp*deta;
            real a = G_alpha(eta, etap);
            real bg = G_beta_g(eta, etap);
            real bd = G_beta_d(eta, etap);
            real gg = G_gamma_g(eta, etap);
            real gd = G_gamma_d(eta, etap);
            real d = G_delta(eta, etap);

            for(int n = 0; n < Nk; n++) {
                g11tmp[n] = 0.6*exp(eta-etap)*exp(a*f[n] - bg*i[n]) + 0.4*exp(-1.5*(eta-etap))*exp(d*g[n] - gd*h[n]);
                g12tmp[n] = 0.4*exp(eta-etap)*exp(a*f[n] - bg*h[n]) - 0.4*exp(-1.5*(eta-etap))*exp(d*f[n] - gd*h[n]);
                g21tmp[n] = 0.6*exp(eta-etap)*exp(a*g[n] + gg*i[n]) - 0.6*exp(-1.5*(eta-etap))*exp(d*g[n] + bd*i[n]);
                g22tmp[n] = 0.4*exp(eta-etap)*exp(a*g[n] - 3/2.*gg*i[n]) + 0.6*exp(-1.5*(eta-etap))*exp(d*f[n] - 2/3.*bd*h[n]);
            }
            g11[Neta*t + tp] = CubicSpline(Nk, k, g11tmp);
            g12[Neta*t + tp] = CubicSpline(Nk, k, g12tmp);
            g21[Neta*t + tp] = CubicSpline(Nk, k, g21tmp);
            g22[Neta*t + tp] = CubicSpline(Nk, k, g22tmp);
        }
    }
}

real RPT::G_11(real k, int t, int tp) const {
    if(t == -1)
        t += Neta;
    return g11[Neta*t + tp](k);
}

real RPT::G_12(real k, int t, int tp) const {
    if(t == -1)
        t += Neta;
    return g12[Neta*t + tp](k);
}

real RPT::G_21(real k, int t, int tp) const {
    if(t == -1)
        t += Neta;
    return g21[Neta*t + tp](k);
}

real RPT::G_22(real k, int t, int tp) const {
    if(t == -1)
        t += Neta;
    return g22[Neta*t + tp](k);
}

real RPT::G_1(real k, int t) const {
    return G_11(k, t) + G_12(k, t);
}

real RPT::G_2(real k, int t) const {
    return G_21(k, t) + G_22(k, t);
}


/************
 * Vertices *
 ************/

/* $\alpha(\vec{q}, \vec{k}-\vec{q}) = \frac{\vec{k}\cdot\vec{q}}{q^2}$ */
real RPT::alpha(real k, real q, real r) {
    real k2 = k*k, q2 = q*q, r2 = r*r;
    return (k2 + q2 - r2)/(2*q2);
}
/* $\beta(\vec{q}, \vec{k}-\vec{q}) = \frac{k^2 \vec{q}\cdot(\vec{k}-\vec{q})}{2 q^2 |\vec{k}-\vec{q}|^2}$ */
real RPT::beta(real k, real q, real r) {
    real k2 = k*k, q2 = q*q, r2 = r*r;
    return k2*(k2 - q2 - r2)/(4*q2*r2);
}


/***************************
 * Mode-coupling integrals *
 ***************************/

/* $I_a(\vec{k},\vec{q};\eta) = \int_0^\eta d\eta_1 ~ G_{ab}(k;\eta,\eta_1) \gamma_{bcd}(\vec{q},\vec{k}-\vec{q}) \tilde G_c(q;\eta_1) \tilde G_d(|\vec{k}-\vec{q}|;\eta_1)$ */
real RPT::I_1(real k, real q, real kq, int t) const {
    if(t == -1)
        t += Neta;
    real gamma112 = alpha(k, kq, q)/2;
    real gamma121 = alpha(k, q, kq)/2;
    real gamma222 = beta(k, q, kq);
    real i1[t+1];
    for(int t1 = 0; t1 <= t; t1++)
        i1[t1] = G_11(k,t,t1)*(gamma112*G_1(q,t1)*G_2(kq,t1) + gamma121*G_2(q,t1)*G_1(kq,t1))
               + G_12(k,t,t1)*gamma222*G_2(q,t1)*G_2(kq,t1);
    return DiscreteIntegrate(t+1, i1, deta);
}

real RPT::I_2(real k, real q, real kq, int t) const {
    if(t == -1)
        t += Neta;
    real gamma112 = alpha(k, kq, q)/2;
    real gamma121 = alpha(k, q, kq)/2;
    real gamma222 = beta(k, q, kq);
    real i2[t+1];
    for(int t1 = 0; t1 <= t; t1++)
        i2[t1] = G_21(k,t,t1)*(gamma112*G_1(q,t1)*G_2(kq,t1) + gamma121*G_2(q,t1)*G_1(kq,t1))
               + G_22(k,t,t1)*gamma222*G_2(q,t1)*G_2(kq,t1);
    return DiscreteIntegrate(t+1, i2, deta);
}

/* $J_a(\vec{k},\vec{p},\vec{q};\eta) = \int_0^\eta d\eta_1 ~ G_{ab}(k;\eta,\eta_1) \gamma_{bcd}(\vec{p},\vec{k}-\vec{p}) I_c(\vec{p},\vec{q};\eta_1) \tilde G_d(|\vec{k}-\vec{p}|;\eta_1)$ */
real RPT::J_1(real k, real p, real kp, real q, real pq, int t) const {
    if(t == -1)
        t += Neta;
    real gamma112 = alpha(k, kp, p)/2;
    real gamma121 = alpha(k, p, kp)/2;
    real gamma222 = beta(k, p, kp);
    real j1[t+1];
    for(int t1 = 0; t1 <= t; t1++) {
        real i1 = I_1(p, q, pq, t1);
        real i2 = I_2(p, q, pq, t1);
        j1[t1] = G_11(k,t,t1)*(gamma112*i1*G_2(kp,t1) + gamma121*i2*G_1(kp,t1))
               + G_12(k,t,t1)*gamma222*i2*G_2(kp,t1);
    }
    return DiscreteIntegrate(t+1, j1, deta);
}

real RPT::J_2(real k, real p, real kp, real q, real pq, int t) const {
    if(t == -1)
        t += Neta;
    real gamma112 = alpha(k, kp, p)/2;
    real gamma121 = alpha(k, p, kp)/2;
    real gamma222 = beta(k, p, kp);
    real j2[t+1];
    for(int t1 = 0; t1 <= t; t1++) {
        real i1 = I_1(p, q, pq, t1);
        real i2 = I_2(p, q, pq, t1);
        j2[t1] = G_21(k,t,t1)*(gamma112*i1*G_2(kp,t1) + gamma121*i2*G_1(kp,t1))
               + G_22(k,t,t1)*gamma222*i2*G_2(kp,t1);
    }
    return DiscreteIntegrate(t+1, j2, deta);
}


/*****************************
 * Tree-level power spectrum *
 *****************************/

/* $P^{(1)}_{ab}(k;\eta) = \tilde G_a(k;\eta) \tilde G_b(k;\eta) P_i(k)$ */
real RPT::P1_11(real k) const {
    return (k == 0) ? 0 : pow2(G_1(k)) * P_i(k);
}
real RPT::P1_12(real k) const {
    return (k == 0) ? 0 : G_1(k)*G_2(k) * P_i(k);
}
real RPT::P1_22(real k) const {
    return (k == 0) ? 0 : pow2(G_2(k)) * P_i(k);
}

real RPT::P1(real k, int a, int b) const {
    switch(a*b) {
        case 1:
            return P1_11(k);
        case 2:
            return P1_12(k);
        case 4:
            return P1_22(k);
        default:
            warning("RPT: invalid indices, a = %d, b = %d\n", a, b);
            return 0;
    }
}

/*****************************
 * 1-loop mode-coupling term *
 *****************************/

/* $P^{(2)}_{ab}(k;\eta) = 2 \int \frac{d^3q}{(2\pi)^3} I_a(\vec{k},\vec{q};\eta,\eta_i) I_b(\vec{k},\vec{q};\eta,\eta_i) P_i(q) P_i(|\vec{k}-\vec{q}|)$ */
real RPT::P2_11(real k) const {
#if 0
    real a[2] = { QMIN, -1 };
    real b[2] = { kcut, +1 };
    return (k == 0) ? 0 : 2*Integrate<2>(bind(&RPT::f2_11, cref(*this), k, std::placeholders::_1, std::placeholders::_2), a, b, 1e-3, 1e-3*P_L(k));
#endif
    real a[2] = { k/M_SQRT2, 0 };
    real b[2] = { kcut, k/M_SQRT2 };
    return (k == 0) ? 0 : 2*Integrate<2>(bind(&RPT::F2_11, cref(*this), k, std::placeholders::_1, std::placeholders::_2), a, b, 1e-3, 1e-3*P_L(k)/2);
}
real RPT::P2_12(real k) const {
#if 0
    real a[2] = { QMIN, -1 };
    real b[2] = { kcut, +1 };
    return (k == 0) ? 0 : 2*Integrate<2>(bind(&RPT::f2_12, cref(*this), k, std::placeholders::_1, std::placeholders::_2), a, b, 1e-4, 1e-4*P_L(k));
#endif
    real a[2] = { k/M_SQRT2, 0 };
    real b[2] = { kcut, k/M_SQRT2 };
    return (k == 0) ? 0 : 2*Integrate<2>(bind(&RPT::F2_12, cref(*this), k, std::placeholders::_1, std::placeholders::_2), a, b, 1e-3, 1e-3*P_L(k)/2);
}
real RPT::P2_22(real k) const {
#if 0
    real a[2] = { QMIN, -1 };
    real b[2] = { kcut, +1 };
    return (k == 0) ? 0 : 2*Integrate<2>(bind(&RPT::f2_22, cref(*this), k, std::placeholders::_1, std::placeholders::_2), a, b, 1e-4, 1e-4*P_L(k));
#endif
    real a[2] = { k/M_SQRT2, 0 };
    real b[2] = { kcut, k/M_SQRT2 };
    return (k == 0) ? 0 : 2*Integrate<2>(bind(&RPT::F2_22, cref(*this), k, std::placeholders::_1, std::placeholders::_2), a, b, 1e-3, 1e-3*P_L(k)/2);
}

real RPT::P2(real k, int a, int b) const {
    switch(a*b) {
        case 1:
            return P2_11(k);
        case 2:
            return P2_12(k);
        case 4:
            return P2_22(k);
        default:
            warning("RPT: invalid indices, a = %d, b = %d\n", a, b);
            return 0;
    }
}

#if 0
real RPT::f2_11(real k, real q, real x) const {
    real kq = sqrt(k*k + q*q - 2*k*q*x);
    return q*q/pow2(2*M_PI) * pow2(I_1(k, q, kq)) * P_i(q)*P_i(kq);
}
real RPT::f2_12(real k, real q, real x) const {
    real kq = sqrt(k*k + q*q - 2*k*q*x);
    return q*q/pow2(2*M_PI) * I_1(k, q, kq)*I_2(k, q, kq) * P_i(q)*P_i(kq);
}
real RPT::f2_22(real k, real q, real x) const {
    real kq = sqrt(k*k + q*q - 2*k*q*x);
    return q*q/pow2(2*M_PI) * pow2(I_2(k, q, kq)) * P_i(q)*P_i(kq);
}
#endif
real RPT::F2_11(real k, real x, real y) const {
    real q = (x+y)/M_SQRT2;
    real r = (x-y)/M_SQRT2;
    return 1/pow2(2*M_PI) * q*r/k * 2*pow2(I_1(k,q,r)) * P_i(q)*P_i(r);
}
real RPT::F2_12(real k, real x, real y) const {
    real q = (x+y)/M_SQRT2;
    real r = (x-y)/M_SQRT2;
    return 1/pow2(2*M_PI) * q*r/k * 2*I_1(k,q,r)*I_2(k,q,r) * P_i(q)*P_i(r);
}
real RPT::F2_22(real k, real x, real y) const {
    real q = (x+y)/M_SQRT2;
    real r = (x-y)/M_SQRT2;
    return 1/pow2(2*M_PI) * q*r/k * 2*pow2(I_2(k,q,r)) * P_i(q)*P_i(r);
}


/*****************************
 * 2-loop mode-coupling term *
 *****************************/

real RPT::P3_11(real k) const {
    if(k == 0)
        return 0;

    real P_L = pow2(D(z)/D(z_i)) * P_i(k);
    double result, error, kd = k;
    int neval;
    p3int.Integrate(&result, &kd, &error, &neval, 0.005, 0.002*P_L);
    return result;
}

real RPT::P3_12(real k) const {
    if(k == 0)
        return 0;

    warning("RPT::P3_12 not implemented yet\n");
    return 0;
}

real RPT::P3_22(real k) const {
    if(k == 0)
        return 0;

    warning("RPT::P3_22 not implemented yet\n");
    return 0;
}

real RPT::P3(real k, int a, int b) const {
    switch(a*b) {
        case 1:
            return P3_11(k);
        case 2:
            return P3_12(k);
        case 4:
            return P3_22(k);
        default:
            warning("RPT: invalid indices, a = %d, b = %d\n", a, b);
            return 0;
    }
}


RPT::ThirdOrderIntegral::ThirdOrderIntegral(const RPT& rpt_)
    : MonteCarloIntegral(5), rpt(rpt_)
{
    maxeval = 1000000;
}

void RPT::ThirdOrderIntegral::Integrand(const double x[], double* f, double* param) const {
    real k = *param;
    real logp = x[0];
    real logq = x[1];
    real mu_p = x[2];
    real mu_q = x[3];
    real phi_q = x[4];
    real p = exp(logp);
    real q = exp(logq);
    real kp = sqrt(k*k + p*p - 2*k*p*mu_p);
    real kq = sqrt(k*k + q*q - 2*k*q*mu_q);
    real pq = sqrt(p*p + q*q - 2*p*q*(sqrt(1-mu_p*mu_p)*sqrt(1-mu_q*mu_q)*cos(phi_q) + mu_p*mu_q));
    f[0] = 16*pow3(p)*pow3(q)/pow5(2*M_PI) * rpt.J_1(k,p,kp,q,pq) * rpt.J_1(k,kq,q,kp,pq) * rpt.P_i(q)*rpt.P_i(kp)*rpt.P_i(pq);
}

real RPT::P_L(real k) const {
    return pow2(D(z)/D(z_i)) * P_i(k);
}
