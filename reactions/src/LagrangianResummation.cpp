#if HAVE_CONFIG_H
# include <config.h>
#endif

#include "GrowthFunction.h"
#include "LagrangianResummation.h"
#include "Quadrature.h"
#include "SpecialFunctions.h"
#include <functional>

using std::cref;
using std::bind;


/* Limits of integration for second-order power spectrum */
const real QMIN = 1e-4;
const real QMAX = 1e1;


LagrangianResummation::LagrangianResummation(const Cosmology& C_, real z_, real epsrel_)
    : C(C_), z(z_), P_L(C_, z_), epsrel(epsrel_)
{
    A = P_L.VelocityDispersion();
    GrowthFunction D(C);
    f = D.f(z);
}

LagrangianResummation::~LagrangianResummation() {
}

/*************************************************/
/***** Matsubara's real space power spectrum *****/
/*************************************************/

/* Inner integrand for P_22 */
static real f_22(const PowerSpectrum& P_L, real k, real logq, real x) {
    if(k <= 0) return 0;
    real q = exp(logq), r = q/k, d = 1 + r*r - 2*r*x;
    if(d < 1e-6)
        return 0;
    else
        return q * P_L(q) * P_L(k*sqrt(d)) * pow2(3*r + 7*x - 10*r*x*x) / pow2(d);
}

/* Integrand for P_13 */
static real f_13(const PowerSpectrum& P_L, real k, real logq) {
    if(k <= 0) return 0;
    real q = exp(logq), r = q/k, r2 = r*r, r3 = r2*r, r4 = r2*r2, r6 = r3*r3;
    real s;
    if(r < 1e-3)
        s = (928./5.)*r2 - (4512./35.)*r4 + (416./21.)*r6;
    else if(fabs(r-1) < 1e-12)
        s = 80;
    else if(r > 1e2)
        s = (352./5.) + (96./5.)/r2 - (160./21.)/r4 - (1376./1155.)/r6;
    else
        s = 12/r2 + 10 + 100*r2 - 42*r4 + 3/r3 * pow3(r2 - 1)*(7*r2 + 2)*log((1+r)/fabs(1-r));
    return q * P_L(q) * s;
}

real LagrangianResummation::P(real k) const {
    return exp(-A*k*k) * (P_L(k) + P22(k) + P13(k));
}

real LagrangianResummation::P22(real k) const {
    if(k <= 0) return 0;
    real a[2] = { log(QMIN), -1 };
    real b[2] = { log(QMAX), +1 };
    real V = pow2(k)/(98*4*M_PI*M_PI);
    return V * Integrate<2>(bind(f_22, cref(P_L), k, std::placeholders::_1, std::placeholders::_2), a, b, epsrel, epsrel*P_L(k)/V);
}

real LagrangianResummation::P13(real k) const {
    if(k <= 0) return 0;
    real a = log(QMIN);
    real b = log(QMAX);
    real V = pow2(k)/(252*4*M_PI*M_PI);
    return V * P_L(k) * Integrate(bind(f_13, cref(P_L), k, std::placeholders::_1), a, b, epsrel, epsrel/V);
}

real LagrangianResummation::G(real k) const {
    return exp(-0.5*A*k*k) * (1 + 0.5*P13(k)/P_L(k));
}


/*****************************************************/
/***** Matsubara's redshift space power spectrum *****/
/*****************************************************/

/* Lower incomplete gamma function times x^{-a} */
static double my_gamma(double a, double x) {
    if(x < 0.1)
        return 1/a - x/(a+1) + pow2(x)/(2*(a+2)) - pow3(x)/(6*(a+3)) + pow4(x)/(24*(a+4)) - pow5(x)/(120*(a+5)) + pow6(x)/(720*(a+6));
    else
        return pow(x, -a) * LowerGamma(a, x);
}

static real f1_22(const PowerSpectrum& P_L, real k, real q, real x) {
    real r = q/k;
    real d = 1 + r*r - 2*r*x;
    if(d < 1e-6)
        return 0;
    else
        return P_L(q) * P_L(k*sqrt(d)) * r*r*pow2(1-x*x) / pow2(d);
}

static real f2_22(const PowerSpectrum& P_L, real k, real q, real x) {
    real r = q/k;
    real d = 1 + r*r - 2*r*x;
    if(d < 1e-6)
        return 0;
    else
        return P_L(q) * P_L(k*sqrt(d)) * (1-x*x)*r*x*(1-r*x) / pow2(d);
}

static real f3_22(const PowerSpectrum& P_L, real k, real q, real x) {
    real r = q/k;
    real d = 1 + r*r - 2*r*x;
    if(d < 1e-6)
        return 0;
    else
        return P_L(q) * P_L(k*sqrt(d)) * x*x*pow2(1-r*x) / pow2(d);
}

static real f4_22(const PowerSpectrum& P_L, real k, real q, real x) {
    real r = q/k;
    real d = 1 + r*r - 2*r*x;
    if(d < 1e-6)
        return 0;
    else
        return P_L(q) * P_L(k*sqrt(d)) * (1-x*x) / pow2(d);
}

static real f1_13(const PowerSpectrum& P_L, real k, real q) {
    real s;
    real r = q/k;
    if(r < 1e-3)
        s = (256./5.)*pow2(r) - (768./35.)*pow4(r) + (256./105.)*pow6(r);
    else if(fabs(r-1) < 1e-20)
        s = 32;
    else if(r > 10)
        s = (256./5.) - (768./35.)/pow2(r) + (256./105.)/pow4(r) + (256./1155.)/pow6(r);
    else
        s = -2/pow2(r)*(1 + pow2(r))*(3 - 14*pow2(r) + 3*pow4(r)) + 3/pow3(r)*pow4(pow2(r) - 1)*log(fabs((1+r)/(1-r)));
    return P_L(q) * s;
}

static real f2_13(const PowerSpectrum& P_L, real k, real q) {
    real s;
    real r = q/k;
    if(r < 1e-3)
        s = (64./5.)*pow2(r) - (576./35.)*pow4(r) + (64./21.)*pow6(r);
    else if(fabs(r-1) < 1e-20)
        s = 0;
    else if(r > 10)
        s = (-64./5.) + (576./35.)/pow2(r) - (64./21.)/pow4(r) - (64./165.)/pow6(r);
    else
        s = 2/pow2(r)*(1 - pow2(r))*(3 - 2*pow2(r) + 3*pow4(r)) + 3/pow3(r)*pow3(pow2(r) - 1)*(1 + pow2(r))*log(fabs((1+r)/(1-r)));
    return P_L(q) * s;
}

real LagrangianResummation::Ps(real k, real mu) const {
    /* Calculate k-dependent coefficients */
    real P = P_L(k);
    real V = pow2(k)/(4*M_PI*M_PI);
    real a[2] = { QMIN, -1 };
    real b[2] = { QMAX, XMAX };
    real B1 = V*Integrate<2>(bind(f1_22, cref(P_L), k, std::placeholders::_1, std::placeholders::_2), a, b, epsrel, epsrel*P/V);
    real B2 = V*Integrate<2>(bind(f2_22, cref(P_L), k, std::placeholders::_1, std::placeholders::_2), a, b, epsrel, epsrel*P/V);
    real B3 = V*Integrate<2>(bind(f3_22, cref(P_L), k, std::placeholders::_1, std::placeholders::_2), a, b, epsrel, epsrel*P/V);
    real B4 = V*Integrate<2>(bind(f4_22, cref(P_L), k, std::placeholders::_1, std::placeholders::_2), a, b, epsrel, epsrel*P/V);
    real C1 = V*Integrate<ExpSub>(bind(f1_13, cref(P_L), k, std::placeholders::_1), QMIN, QMAX, epsrel, epsrel/V);
    real C2 = V*Integrate<ExpSub>(bind(f2_13, cref(P_L), k, std::placeholders::_1), QMIN, QMAX, epsrel, epsrel/V);

    /* Calculate E_nm coefficients */
    real E[5][5];
    for(int n = 0; n < 5; n++)
        for(int m = 0; m < 5; m++)
            E[n][m] = 0;
    E[0][0] = (9./98)*B1 + (3./7)*B2 + (1./2)*B3 + (P/48)*( (10./21)*C1 + (6./7)*C2 );
    E[1][1] = 4*E[0][0];
    E[1][2] = (-3./14)*B1 - (3./2)*B2 + (1./4)*B4 - (P/48)*(6./7)*C1;
    E[2][2] = (57./98)*B1 + (51./14)*B2 + 3*B3 - (1./4)*B4 + (P/48)*( (16./7)*C1 + (30./7)*C2 );
    E[2][3] = (-3./7)*B1 - 3*B2 + (1./2)*B4 - (P/48)*(6./7)*C1;
    E[2][4] = (3./16)*B1;
    E[3][3] = (3./7)*B1 + (27./7)*B2 + 2*B3 - (1./2)*B4 + (P/48)*( (6./7)*C1 + (12./7)*C2 );
    E[3][4] = (-3./8)*B1 - (3./2)*B2 + (1./4)*B4;
    E[4][4] = (3./16)*B1 + (3./2)*B2 + (1./2)*B3 - (1./4)*B4;

    real mu2n[5] = { 1, pow2(mu), pow4(mu), pow6(mu), pow8(mu) };
    real fm[5] = { 1, f, pow2(f), pow3(f), pow4(f) };

    real S = 0;
    for(int n = 0; n <= 4; n++)
        for(int m = 0; m <= 4; m++)
            S += mu2n[n] * fm[m] * E[n][m];

    return exp(-A*k*k*(1 + f*(f+2)*mu*mu)) * (pow2(1+f*mu*mu)*P + S);
}

/* Factorial function */
static real fact(int n) {
    return (n == 0) ? 1 : n*fact(n-1);
}

real LagrangianResummation::Ps_ell(real k, int ell) const {
    /* Odd moments vanish */
    if((ell % 2) == 1)
        return 0;

    /* Compute even coefficients of Legendre polynomial */
    real L[ell/2 + 1];
    for(int j = 0; j <= ell/2; j++) {
        L[j] = pow(2.0, -ell) * pow(-1.0, j) * fact(2*ell-2*j)/(fact(j)*fact(ell-j)*fact(ell-2*j));
    }

    /* Calculate k-dependent coefficients */
    real P = P_L(k);
    real V = pow2(k)/(4*M_PI*M_PI);
    real a[2] = { QMIN, -1 };
    real b[2] = { QMAX, XMAX };
    real B1 = V*Integrate<2>(bind(f1_22, cref(P_L), k, std::placeholders::_1, std::placeholders::_2), a, b, epsrel, epsrel*P/V);
    real B2 = V*Integrate<2>(bind(f2_22, cref(P_L), k, std::placeholders::_1, std::placeholders::_2), a, b, epsrel, epsrel*P/V);
    real B3 = V*Integrate<2>(bind(f3_22, cref(P_L), k, std::placeholders::_1, std::placeholders::_2), a, b, epsrel, epsrel*P/V);
    real B4 = V*Integrate<2>(bind(f4_22, cref(P_L), k, std::placeholders::_1, std::placeholders::_2), a, b, epsrel, epsrel*P/V);
    real C1 = V*Integrate<ExpSub>(bind(f1_13, cref(P_L), k, std::placeholders::_1), QMIN, QMAX, epsrel, epsrel/V);
    real C2 = V*Integrate<ExpSub>(bind(f2_13, cref(P_L), k, std::placeholders::_1), QMIN, QMAX, epsrel, epsrel/V);

    /* Calculate E_nm */
    real E[5][5];
    for(int n = 0; n < 5; n++)
        for(int m = 0; m < 5; m++)
            E[n][m] = 0;
    E[0][0] = (9./98)*B1 + (3./7)*B2 + (1./2)*B3 + (P/48)*( (10./21)*C1 + (6./7)*C2 );
    E[1][1] = 4*E[0][0];
    E[1][2] = (-3./14)*B1 - (3./2)*B2 + (1./4)*B4 - (P/48)*(6./7)*C1;
    E[2][2] = (57./98)*B1 + (51./14)*B2 + 3*B3 - (1./4)*B4 + (P/48)*( (16./7)*C1 + (30./7)*C2 );
    E[2][3] = (-3./7)*B1 - 3*B2 + (1./2)*B4 - (P/48)*(6./7)*C1;
    E[2][4] = (3./16)*B1;
    E[3][3] = (3./7)*B1 + (27./7)*B2 + 2*B3 - (1./2)*B4 + (P/48)*( (6./7)*C1 + (12./7)*C2 );
    E[3][4] = (-3./8)*B1 - (3./2)*B2 + (1./4)*B4;
    E[4][4] = (3./16)*B1 + (3./2)*B2 + (1./2)*B3 - (1./4)*B4;

    /* Calculate Legendre moments */
    real x = f*(f+2)*A*k*k;
    real Zn[5];
    for(int n = 0; n < 5; n++) {
        real s = 0;
        for(int j = 0; j <= ell/2; j++)
            s += L[j] * my_gamma(n+(2*j+1.)/2., x);
        Zn[n] = s;
    }

    /* Perform sum */
    real S = 0;
    real fm[5] = { 1, f, pow2(f), pow3(f), pow4(f) };
    for(int n = 0; n <= 4; n++)
        for(int m = 0; m <= 4; m++)
            S += Zn[n] * fm[m] * E[n][m];

    return (2*ell+1)/2. * exp(-A*k*k) * ((Zn[0] + 2*Zn[1]*f + Zn[2]*f*f)*P + S);
}
