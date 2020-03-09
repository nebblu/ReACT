#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <algorithm>
#include <cassert>

#include "FlowingWithTime.h"
#include "GrowthFunction.h"
#include "InterpolatedPS.h"
#include "Quadrature.h"
#include "ODE.h"
#include <functional>

using std::cref;
using std::bind;


/* My translation of code <-> theory:
 *   X(0,i)  = $P_11(k[i])$
 *   X(1,i)  = $P_12(k[i])$
 *   X(2,i)  = $P_22(k[i])$
 *   X(3,i)  = $I_{112,111}(k[i])$
 *   X(4,i)  = $I_{112,112}(k[i])$
 *   X(5,i)  = $I_{112,122}(k[i])$
 *   X(6,i)  = $I_{112,211}(k[i])$
 *   X(7,i)  = $I_{112,212}(k[i])$
 *   X(8,i)  = $I_{112,222}(k[i])$
 *   X(9,i)  = $I_{222,111}(k[i])$
 *   X(10,i) = $I_{222,112}(k[i])$
 *   X(11,i) = $I_{222,122}(k[i])$
 *   X(12,i) = $I_{222,211}(k[i])$
 *   X(13,i) = $I_{222,212}(k[i])$
 *   X(14,i) = $I_{222,222}(k[i])$
 *
 * It's not obvious from the definition that $I_{acd,bef}(k)$ is symmetric
 * in its last two indices $ef$.  This result follows from the fact that
 * it is initially symmetric ($I_{acd,bef} = 0$ at $\eta = 0$) and the
 * equations of motion preserve this symmetry.
 *
 * See Appendix B of Pietroni's paper for an explanation of the code below. */

/* Tolerance parameter for mode-coupling integrals [Eq (B.6)] */
const real EPSREL = 1e-4;
const real XMAX = 25;


FlowingWithTime::FlowingWithTime(const Cosmology& C_, real z_i, const PowerSpectrum& P_i_, const array& k_)
    : C(C_), P_i(P_i_), k(k_)
{
    a_i = 1/(1 + z_i);
    Nk = k.size();
}

FlowingWithTime::~FlowingWithTime() {
}

static real zero = 0;

InterpolatedP_ab::InterpolatedP_ab(const Cosmology& C, int N, const real* k, const real* p11, const real* p12, const real* p22)
    : P_11(C, N, k, p11), P_12(C, N, k, p12), P_22(C, N, k, p22)
{
    ka = k[0];
    kb = k[N-1];
    p11a = p11[0];
    p11b = p11[N-1];
    p12a = p12[0];
    p12b = p12[N-1];
    p22a = p22[0];
    p22b = p22[N-1];
}

real InterpolatedP_ab::operator()(int a, int b, real k) const {
    switch(a*b) {
        case 1:
            if(k <= ka)
                return exp(-pow2((k-ka)/k)) * p11a;
            else if(k >= kb)
                return exp(-4*pow2((k-kb)/kb)) * p11b;
            else
                return P_11(k);
        case 2:
            if(k <= ka)
                return exp(-pow2((k-ka)/k)) * p12a;
            else if(k >= kb)
                return exp(-4*pow2((k-kb)/kb)) * p12b;
            else
                return P_12(k);
        case 4:
            if(k <= ka)
                return exp(-pow2((k-ka)/k)) * p22a;
            else if(k >= kb)
                return exp(-4*pow2((k-kb)/kb)) * p22b;
            else
                return P_22(k);
//        case 1: return P_11(k);
//        case 2: return P_12(k);
//        case 4: return P_22(k);
        default:
            warning("[InterpolatedP_ab] index error: a = %d, b = %d\n", a, b);
            return 0;
    }
}

/* $P_{ab}(k[i])$ */
struct IndexedP_ab {
    vector<real>& X;
    int Nk;

    IndexedP_ab(vector<real>& X_) : X(X_) {
        Nk = X.size()/15;
    }

#if 0
    real operator()(int a, int b, int i) const {
        switch(a*b) {
            case 1: return exp(X[i]);
            case 2: return exp(X[Nk+i]);
            case 4: return exp(X[2*Nk+i]);
            default:
                warning("[IndexedP_ab] index error: a = %d, b = %d\n", a, b);
                return 0;
        }
    }
#endif

    real& operator()(int a, int b, int i) {
        switch(a*b) {
            case 1: return X[i];
            case 2: return X[Nk+i];
            case 4: return X[2*Nk+i];
            default:
                warning("[IndexedP_ab] index error: a = %d, b = %d\n", a, b);
                return zero;
        }
    }
};

/* Mapping of $(acd,bef)$ into X[n*Nk+i] */
static int m_to_n[64] = {
    /* 111,111 */  -1,  /* 0 */
    /* 111,112 */  -1,  /* 0 */
    /* 111,121 */  -1,  /* 0 */
    /* 111,122 */  -1,  /* 0 */
    /* 111,211 */  -1,  /* 0 */
    /* 111,212 */  -1,  /* 0 */
    /* 111,221 */  -1,  /* 0 */
    /* 111,222 */  -1,  /* 0 */
    /* 112,111 */  3,   /* X(3) */
    /* 112,112 */  4,   /* X(4) */
    /* 112,121 */  4,   /* = 112,112 */
    /* 112,122 */  5,   /* X(5) */
    /* 112,211 */  6,   /* X(6) */
    /* 112,212 */  7,   /* X(7) */
    /* 112,221 */  7,   /* = 112,212 */
    /* 112,222 */  8,   /* X(8) */
    /* 121,111 */  3,   /* = 112,111 */
    /* 121,112 */  4,   /* = 112,112 */
    /* 121,121 */  4,   /* = 112,112 */
    /* 121,122 */  5,   /* = 112,122 */
    /* 121,211 */  6,   /* = 112,211 */
    /* 121,212 */  7,   /* = 112,212 */
    /* 121,221 */  7,   /* = 112,212 */
    /* 121,222 */  8,   /* = 112,222 */
    /* 122,111 */  -1,  /* 0 */
    /* 122,112 */  -1,  /* 0 */
    /* 122,121 */  -1,  /* 0 */
    /* 122,122 */  -1,  /* 0 */
    /* 122,211 */  -1,  /* 0 */
    /* 122,212 */  -1,  /* 0 */
    /* 122,221 */  -1,  /* 0 */
    /* 122,222 */  -1,  /* 0 */
    /* 211,111 */  -1,  /* 0 */
    /* 211,112 */  -1,  /* 0 */
    /* 211,121 */  -1,  /* 0 */
    /* 211,122 */  -1,  /* 0 */
    /* 211,211 */  -1,  /* 0 */
    /* 211,212 */  -1,  /* 0 */
    /* 211,221 */  -1,  /* 0 */
    /* 211,222 */  -1,  /* 0 */
    /* 212,111 */  -1,  /* 0 */
    /* 212,112 */  -1,  /* 0 */
    /* 212,121 */  -1,  /* 0 */
    /* 212,122 */  -1,  /* 0 */
    /* 212,211 */  -1,  /* 0 */
    /* 212,212 */  -1,  /* 0 */
    /* 212,221 */  -1,  /* 0 */
    /* 212,222 */  -1,  /* 0 */
    /* 221,111 */  -1,  /* 0 */
    /* 221,112 */  -1,  /* 0 */
    /* 221,121 */  -1,  /* 0 */
    /* 221,122 */  -1,  /* 0 */
    /* 221,211 */  -1,  /* 0 */
    /* 221,212 */  -1,  /* 0 */
    /* 221,221 */  -1,  /* 0 */
    /* 221,222 */  -1,  /* 0 */
    /* 222,111 */  9,   /* X(9) */
    /* 222,112 */  10,  /* X(10) */
    /* 222,121 */  10,  /* = 222,112 */
    /* 222,122 */  11,  /* X(11) */
    /* 222,211 */  12,  /* X(12) */
    /* 222,212 */  13,  /* X(13) */
    /* 222,221 */  13,  /* = 222,212 */
    /* 222,222 */  14   /* X(14) */
};

/* $I_{acd,bef}(k[i])$ */
struct IndexedI_acdbef {
    vector<real>& X;
    int Nk;

    IndexedI_acdbef(vector<real>& X_) : X(X_) {
        Nk = X.size()/15;
    }

    real& operator()(int a, int c, int d, int b, int e, int f, int i) {
        assert(zero == 0);      // make sure we haven't somehow changed zero

        a--; b--; c--; d--; e--; f--;
        assert((a | b | c | d | e | f) < 2);  // make sure indices are 1 or 2
        int m = (a << 5) + (c << 4) + (d << 3) + (b << 2) + (e << 1) + (f << 0);
        int n = m_to_n[m];
        if(n < 0)
            return zero;
        else
            return X[n*Nk+i];
    }
};

vector<InterpolatedP_ab> FlowingWithTime::CalculateP_ab(int Nz, const real zin[]) const {
    /* Copy and sort redshift list */
    array z(Nz);
    for(int i = 0; i < Nz; i++)
        z[i] = zin[i];
    std::sort(&z[0], &z[Nz]);
    std::reverse(&z[0], &z[Nz]);

    /* Make sure redshifts are now sequential */
    for(int j = 0; j < Nz-1; j++)
        assert(z[j] > z[j+1]);

    /* Output list */
    vector<InterpolatedP_ab> output;

    GrowthFunction D(C);
    array a(Nz);
    for(int j = 0; j < Nz; j++)
        a[j] = 1/(1 + z[j]);

    /* Initialize X */
    vector<real> X(15*Nk);
    for(int i = 0; i < Nk; i++)
        X[i] = X[Nk+i] = X[2*Nk+i] = log(P_i(k[i]));

    real h_i = log(a[0]/a_i) / 100;   // initial step size

    /* Start off with fixed time steps, to avoid some weirdness where the
       adaptive stepper can't get the integration started */
    vector<vector<real> > Y = RungeKutta4(0, 10*h_i, X, bind(&FlowingWithTime::dXdeta, cref(*this), std::placeholders::_1, std::placeholders::_2), 10);
    for(int j = 0; j < 15*Nk; j++)
        X[j] = Y[j][9];

    real eta = 10*h_i;
    array p11(Nk), p12(Nk), p22(Nk);
    for(int j = 0; j < Nz; j++) {
        /* Integrate forward to time a[j] */
        real etaf = log(a[j]/a_i);
        int err = RKDP(eta, etaf, X, bind(&FlowingWithTime::dXdeta, cref(*this), std::placeholders::_1, std::placeholders::_2), X, 1e-3, h_i);

        if(err < 0)
            error("Aborting following error in RKDP.\n");

        /* Save power spectrum */
        real f = D.f(z[j]);
        IndexedP_ab logP(X);
        for(int i = 0; i < Nk; i++) {
            p11[i] = exp(logP(1,1, i)) * pow2(a[j]/a_i);
            p12[i] = exp(logP(1,2, i)) * pow2(a[j]/a_i) / f;
            p22[i] = exp(logP(2,2, i)) * pow2(a[j]/a_i) / (f*f);
        }
        output.push_back(InterpolatedP_ab(C, Nk, &k[0], p11, p12, p22));

        /* Prepare for next step */
        eta = etaf;
    }

    return output;
}

template<int m, int n>
struct Matrix {
    real a[m*n];

    real& operator()(int i, int j) { return a[(j-1)*n+(i-1)]; }
};

vector<real> FlowingWithTime::dXdeta(real eta, vector<real>& X) const {
    real a = a_i * exp(eta);
    debug("\neta = %g, a = %g\n", eta, a);
    Matrix<2,2> Omega;
    Omega(1,1) = 1;
    Omega(1,2) = -1;
    Omega(2,1) = -1.5 * C.Omega_m*pow2(C.H0)/(pow3(a)*pow2(C.H(a)));
    Omega(2,2) = 3 + a/C.H(a) * C.dHda(a);
    debug("Omega = %g %g %g %g\n", Omega(1,1), Omega(1,2), Omega(2,1), Omega(2,2));

    vector<real> dXdeta(15*Nk);
    IndexedP_ab logP(X);
    IndexedI_acdbef I(X);
    IndexedP_ab dlogPdeta(dXdeta);
    IndexedI_acdbef dIdeta(dXdeta);

    array p11(Nk), p12(Nk), p22(Nk);
    for(int i = 0; i < Nk; i++) {
        p11[i] = exp(logP(1,1, i));
        p12[i] = exp(logP(1,2, i));
        p22[i] = exp(logP(2,2, i));
    }
    InterpolatedP_ab Pint(C, Nk, &k[0], p11, p12, p22);

    int i, c, d, g;
    #pragma omp parallel default(shared) private(i,c,d,g)
    #pragma omp for schedule(dynamic)
    for(i = 0; i < Nk; i++) {
        for(int n = 0; n < 3; n++)
            debug("%.2e ", exp(X[n*Nk+i]));
//        for(int n = 3; n < 15; n++)
//            debug("%.2e ", X[n*Nk+i]);
        debug("\n");

        /* $\partial_\eta P_{ab}(k) = -\Omega_{ac} P_{cb}(k) - \Omega_{bc} P_{ac}(k) + e^\eta \frac{4\pi}{k} [I_{acd,bcd}(k) + I_{bcd,acd}(k)]$ */
        dlogPdeta(1,1, i) = dlogPdeta(1,2, i) = dlogPdeta(2,2, i) = 0;
        for(c = 1; c <= 2; c++) {
            dlogPdeta(1,1, i) += -Omega(1,c)*exp(logP(c,1, i)) - Omega(1,c)*exp(logP(1,c, i));
            dlogPdeta(1,2, i) += -Omega(1,c)*exp(logP(c,2, i)) - Omega(2,c)*exp(logP(1,c, i));
            dlogPdeta(2,2, i) += -Omega(2,c)*exp(logP(c,2, i)) - Omega(2,c)*exp(logP(2,c, i));
            for(d = 1; d <= 2; d++) {
                dlogPdeta(1,1, i) += exp(eta) * 4*M_PI/k[i] * (I(1,c,d,1,c,d, i) + I(1,c,d,1,c,d, i));
                dlogPdeta(1,2, i) += exp(eta) * 4*M_PI/k[i] * (I(1,c,d,2,c,d, i) + I(2,c,d,1,c,d, i));
                dlogPdeta(2,2, i) += exp(eta) * 4*M_PI/k[i] * (I(2,c,d,2,c,d, i) + I(2,c,d,2,c,d, i));
            }
        }
        dlogPdeta(1,1, i) /= exp(logP(1,1, i));
        dlogPdeta(1,2, i) /= exp(logP(1,2, i));
        dlogPdeta(2,2, i) /= exp(logP(2,2, i));

        /* $\partial_\eta I_{acd,bef}(k) = -\Omega_{bg} I_{acd,gef}(k) - \Omega_{eg} I_{acd,bgf}(k) - \Omega_{fg} I_{acd,beg}(k) + 2 e^\eta A_{acd,bef}(k)$ */
        dIdeta(1,1,2,1,1,1, i) = 2*exp(eta) * A(Pint, 1,1,2,1,1,1, k[i]);
        dIdeta(1,1,2,1,1,2, i) = 2*exp(eta) * A(Pint, 1,1,2,1,1,2, k[i]);
        dIdeta(1,1,2,1,2,2, i) = 2*exp(eta) * A(Pint, 1,1,2,1,2,2, k[i]);
        dIdeta(1,1,2,2,1,1, i) = 2*exp(eta) * A(Pint, 1,1,2,2,1,1, k[i]);
        dIdeta(1,1,2,2,1,2, i) = 2*exp(eta) * A(Pint, 1,1,2,2,1,2, k[i]);
        dIdeta(1,1,2,2,2,2, i) = 2*exp(eta) * A(Pint, 1,1,2,2,2,2, k[i]);
        dIdeta(2,2,2,1,1,1, i) = 2*exp(eta) * A(Pint, 2,2,2,1,1,1, k[i]);
        dIdeta(2,2,2,1,1,2, i) = 2*exp(eta) * A(Pint, 2,2,2,1,1,2, k[i]);
        dIdeta(2,2,2,1,2,2, i) = 2*exp(eta) * A(Pint, 2,2,2,1,2,2, k[i]);
        dIdeta(2,2,2,2,1,1, i) = 2*exp(eta) * A(Pint, 2,2,2,2,1,1, k[i]);
        dIdeta(2,2,2,2,1,2, i) = 2*exp(eta) * A(Pint, 2,2,2,2,1,2, k[i]);
        dIdeta(2,2,2,2,2,2, i) = 2*exp(eta) * A(Pint, 2,2,2,2,2,2, k[i]);
        for(g = 1; g <= 2; g++) {
            dIdeta(1,1,2,1,1,1, i) += -Omega(1,g)*I(1,1,2,g,1,1, i) - Omega(1,g)*I(1,1,2,1,g,1, i) - Omega(1,g)*I(1,1,2,1,1,g, i);
            dIdeta(1,1,2,1,1,2, i) += -Omega(1,g)*I(1,1,2,g,1,2, i) - Omega(1,g)*I(1,1,2,1,g,2, i) - Omega(2,g)*I(1,1,2,1,1,g, i);
            dIdeta(1,1,2,1,2,2, i) += -Omega(1,g)*I(1,1,2,g,2,2, i) - Omega(2,g)*I(1,1,2,1,g,2, i) - Omega(2,g)*I(1,1,2,1,2,g, i);
            dIdeta(1,1,2,2,1,1, i) += -Omega(2,g)*I(1,1,2,g,1,1, i) - Omega(1,g)*I(1,1,2,2,g,1, i) - Omega(1,g)*I(1,1,2,2,1,g, i);
            dIdeta(1,1,2,2,1,2, i) += -Omega(2,g)*I(1,1,2,g,1,2, i) - Omega(1,g)*I(1,1,2,2,g,2, i) - Omega(2,g)*I(1,1,2,2,1,g, i);
            dIdeta(1,1,2,2,2,2, i) += -Omega(2,g)*I(1,1,2,g,2,2, i) - Omega(2,g)*I(1,1,2,2,g,2, i) - Omega(2,g)*I(1,1,2,2,2,g, i);
            dIdeta(2,2,2,1,1,1, i) += -Omega(1,g)*I(2,2,2,g,1,1, i) - Omega(1,g)*I(2,2,2,1,g,1, i) - Omega(1,g)*I(2,2,2,1,1,g, i);
            dIdeta(2,2,2,1,1,2, i) += -Omega(1,g)*I(2,2,2,g,1,2, i) - Omega(1,g)*I(2,2,2,1,g,2, i) - Omega(2,g)*I(2,2,2,1,1,g, i);
            dIdeta(2,2,2,1,2,2, i) += -Omega(1,g)*I(2,2,2,g,2,2, i) - Omega(2,g)*I(2,2,2,1,g,2, i) - Omega(2,g)*I(2,2,2,1,2,g, i);
            dIdeta(2,2,2,2,1,1, i) += -Omega(2,g)*I(2,2,2,g,1,1, i) - Omega(1,g)*I(2,2,2,2,g,1, i) - Omega(1,g)*I(2,2,2,2,1,g, i);
            dIdeta(2,2,2,2,1,2, i) += -Omega(2,g)*I(2,2,2,g,1,2, i) - Omega(1,g)*I(2,2,2,2,g,2, i) - Omega(2,g)*I(2,2,2,2,1,g, i);
            dIdeta(2,2,2,2,2,2, i) += -Omega(2,g)*I(2,2,2,g,2,2, i) - Omega(2,g)*I(2,2,2,2,g,2, i) - Omega(2,g)*I(2,2,2,2,2,g, i);
        }
    }

    return dXdeta;
}

static real gamma(int a, int b, int c, real k2, real q2, real p2) {
    if(a == 1) {
        if(b == 1 && c == 2)
            return (k2 - q2 + p2)/(4*p2);
        else if(b == 2 && c == 1)
            return (k2 + q2 - p2)/(4*q2);
        else
            return 0;
    }
    else if(a == 2 && b == 2 && c == 2)
        return k2*(k2 - q2 - p2)/(4*q2*p2);
    else
        return 0;
}

static bool nonzero(int a, int b, int c) {
    return (a == 1 && b*c == 2) || (a*b*c == 8);
}

static real F_A(const InterpolatedP_ab& P, int a, int c, int d, int b, int e, int f, real k, real q, real p) {
    real k2 = k*k, q2 = q*q, p2 = p*p;
    real S = 0;
    for(int g = 1; g <= 2; g++)
        for(int h = 1; h <= 2; h++)
            S += (nonzero(b,g,h) ? gamma(b,g,h, k2,q2,p2) * P(g,e, q) * P(h,f, p) : 0)
               + (nonzero(e,g,h) ? gamma(e,g,h, q2,p2,k2) * P(g,f, p) * P(h,b, k) : 0)
               + (nonzero(f,g,h) ? gamma(f,g,h, p2,k2,q2) * P(g,b, k) * P(h,e, q) : 0);
    return 0.5*gamma(a,c,d, k2,q2,p2) * S;
}

static real f_A(const InterpolatedP_ab& P, int* acdbef, real k, real x, real y) {
    int a = acdbef[0];
    int c = acdbef[1];
    int d = acdbef[2];
    int b = acdbef[3];
    int e = acdbef[4];
    int f = acdbef[5];
    real q = (x+y)/M_SQRT2;
    real p = (x-y)/M_SQRT2;
    return q*p * ( F_A(P, a,c,d,b,e,f, k,q,p) + F_A(P, a,c,d,b,e,f, k,p,q) );
}

real FlowingWithTime::A(const InterpolatedP_ab& P, int a, int c, int d, int b, int e, int f, real k1) const {
    int acdbef[6] = { a, c, d, b, e, f };
    real min[2] = { k1/M_SQRT2, 0 };
    real max[2] = { XMAX, k1/M_SQRT2 - 1e-5 };
    return 1/pow3(2*M_PI) * Integrate<2>(bind(f_A, cref(P), acdbef, k1, std::placeholders::_1, std::placeholders::_2), min, max, EPSREL);
}
