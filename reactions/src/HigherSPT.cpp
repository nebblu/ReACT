#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <algorithm>

#include "HigherSPT.h"
#include "MonteCarlo.h"
#include "PowerSpectrum.h"


/* The P2 and P3 contributions are computed with a Monte-Carlo integrator.
 * The bounds on the error are
 *   |\Delta P2| < max{ 0.005*|P2|, 0.001*|P_L| }
 *   |\Delta P3| < max{ 0.005*|P3|, 0.001*|P_L| }
 * This allows the integration to exit early if the higher order terms are
 * already being swamped by the linear order.
 */

/* vec3: a simple 3-dimensional vector class */
struct vec3 {
    double x, y, z;

    vec3() { x = y = z = 0; }
    vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    vec3(const vec3& q) : x(q.x), y(q.y), z(q.z) {}

    void operator+=(const vec3& p) {
        x += p.x; y += p.y; z += p.z;
    }
    void operator-=(const vec3& p) {
        x -= p.x; y -= p.y; z -= p.z;
    }
    void operator*=(double a) {
        x *= a; y *= a; z *= a;
    }
    void operator/=(double a) {
        x /= a; y /= a; z /= a;
    }
    double square() const {
        return x*x + y*y + z*z;
    }
    double mag() const {
        return sqrt(square());
    }
};

vec3 operator-(const vec3& p) {
    return vec3(-p.x, -p.y, -p.z);
}

vec3 operator+(const vec3& p, const vec3& q) {
    return vec3(p.x + q.x, p.y + q.y, p.z + q.z);
}

vec3 operator-(const vec3& p, const vec3& q) {
    return vec3(p.x - q.x, p.y - q.y, p.z - q.z);
}

vec3 operator*(real a, const vec3& p) {
    return vec3(a*p.x, a*p.y, a*p.z);
}

vec3 operator*(const vec3& p, real a) {
    return vec3(p.x*a, p.y*a, p.z*a);
}

vec3 operator/(const vec3& p, real a) {
    return vec3(p.x/a, p.y/a, p.z/a);
}

double operator*(const vec3& p, const vec3& q) {
    return p.x*q.x + p.y*q.y + p.z*q.z;
}


/* n! up to n = 10 */
static const int nfactorial[] = { 1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800 };

/* Recursion relations for SPT kernels */
static double F(int n, vec3 q[], int order[], const vec3& k);
static double G(int n, vec3 q[], int order[], const vec3& k);

static double F(int n, vec3 q[], int order[], const vec3& k) {
    if(n == 1)
        return 1;

    double k1k1, k1k2, k2k2, kk = k.square();
    double S = 0;
    vec3 k1, k2;
    for(int m = 1; m <= n-1; m++) {
        k1 += q[order[m-1]];
        k2 = k - k1;
        k1k1 = k1.square();
        k2k2 = k2.square();
        if(k1k1 < 1e-10 || k2k2 < 1e-10)
            continue;
        k1k2 = 0.5*(kk - k1k1 - k2k2);
        S += G(m, q, order, k1) * ( (1+2*n)*(k*k1)/k1k1 * F(n-m, q, order+m, k2)
                                  + kk*k1k2/(k1k1*k2k2) * G(n-m, q, order+m, k2) );
    }
    return S/((2*n+3)*(n-1));
}

static double G(int n, vec3 q[], int order[], const vec3& k) {
    if(n == 1)
        return 1;

    double k1k1, k1k2, k2k2, kk = k.square();
    double S = 0;
    vec3 k1, k2;
    for(int m = 1; m <= n-1; m++) {
        k1 += q[order[m-1]];
        k2 = k - k1;
        k1k1 = k1.square();
        k2k2 = k2.square();
        if(k1k1 < 1e-10 || k2k2 < 1e-10)
            continue;
        k1k2 = 0.5*(kk - k1k1 - k2k2);
        S += G(m, q, order, k1) * ( 3*(k*k1)/k1k1 * F(n-m, q, order+m, k2)
                                  + n*kk*k1k2/(k1k1*k2k2) * G(n-m, q, order+m, k2) );
    }
    return S/((2*n+3)*(n-1));
}

/* Symmetrized kernels */
static double Fs(int n, vec3 q[], const vec3& k) {
    double S = 0;               // accumulates sum of F_n over all permutations
    int order[n];               // changes in place as permutations are iterated over
    for(int i = 0; i < n; i++)
        order[i] = i;
    do {
        S += F(n, q, order, k);
    } while(std::next_permutation(order, order+n));

    return S/nfactorial[n];
}

static double Gs(int n, vec3 q[], const vec3& k) {
    double S = 0;               // accumulates sum of G_n over all permutations
    int order[n];               // changes in place as permutations are iterated over
    for(int i = 0; i < n; i++)
        order[i] = i;
    do {
        S += G(n, q, order, k);
    } while(std::next_permutation(order, order+n));

    return S/nfactorial[n];
}

static double F2s(const vec3& k, const vec3& q1, const vec3& q2) {
    vec3 q[2] = { q1, q2 };
    return Fs(2, q, k);
}

static double F3s(const vec3& k, const vec3& q1, const vec3& q2, const vec3& q3) {
    vec3 q[3] = { q1, q2, q3 };
    return Fs(3, q, k);
}

static double F4s(const vec3& k, const vec3& q1, const vec3& q2, const vec3& q3, const vec3& q4) {
    vec3 q[4] = { q1, q2, q3, q4 };
    return Fs(4, q, k);
}

static double F5s(const vec3& k, const vec3& q1, const vec3& q2, const vec3& q3, const vec3& q4, const vec3& q5) {
    vec3 q[5] = { q1, q2, q3, q4, q5 };
    return Fs(5, q, k);
}

/* Second order mode-coupling contribution */
HigherSPT::P2_Integral::P2_Integral(const PowerSpectrum& P_L_, real qmin, real qmax)
    : MonteCarloIntegral(2, 2), P_L(P_L_)
{
    Range domain[2] = {
        { log(qmin), log(qmax) },
        { -1, 1 }
    };
    SetDomain(domain);
    maxeval = 10000000; // 10 million
}

/* \begin{align}
 *    P^{(2)}(k) &= \int \frac{d^3q}{(2\pi)^3} \left[ 2 P_0(q) P_0(|\vec{k}-\vec{q}|) [F_2^{(s)}(\vec{q},\vec{k}-\vec{q})]^2 + 6 P_0(k) P_0(q) F_3^{(s)}(\vec{q},-\vec{q},\vec{k} \right] \\
 *               &= \int_{-\infty}^\infty d\log q \int_{-1}^1 d\mu \frac{q^3}{4\pi^2} M^{(2)}(\vec{k},\vec{q})
 * \end{align} */
void HigherSPT::P2_Integral::Integrand(const double x[], double* f, double* param) const {
    double qmag = exp(x[0]);
    double mu = x[1];
    double kmag = *param;

    vec3 k(0, 0, kmag);
    vec3 q(qmag*sqrt(1-mu*mu), 0, qmag*mu);
    f[0] = 6*P_L(kmag)*P_L(qmag) * F3s(k, k, q, -q);              // P^{(13)}
    f[1] = 2*P_L(qmag)*P_L((k-q).mag()) * pow2(F2s(k, q, k-q));   // P^{(22)}

    f[0] *= pow3(qmag)/pow2(2*M_PI);
    f[1] *= pow3(qmag)/pow2(2*M_PI);
}


/* Third order mode-coupling contribution */
HigherSPT::P3_Integral::P3_Integral(const PowerSpectrum& P_L_, real qmin, real qmax)
    : MonteCarloIntegral(5, 4), P_L(P_L_)
{
    Range domain[5] = {
        { log(qmin), log(qmax) },
        { log(qmin), log(qmax) },
        { -1, 1 },
        { -1, 1 },
        { 0, 2*M_PI }
    };
    SetDomain(domain);
    maxeval = 100000000; // 100 million
}

/* \begin{align}
 *    P^{(3)}(k) &= \int \frac{d^3q}{(2\pi)^3} \frac{d^3p}{(2\pi)^3} M^{(3)}(\vec{k},\vec{q},\vec{p}) \\
 *               &= \frac{1}{(2\pi)^5} \int_{-\infty}^\infty d\log q \int_{-\infty}^\infty d\log p \int_{-1}^1 d\mu_q \int_{-1}^1 d\mu_p \int_0^{2\pi} d\phi_p  M^{(3)}(\vec{k},\vec{q},\vec{p})
 * \end{align} */
void HigherSPT::P3_Integral::Integrand(const double x[], double* f, double* param) const {
    double qmag = exp(x[0]);
    double pmag = exp(x[1]);
    double mu_q = x[2];
    double mu_p = x[3];
    double phi_p = x[4];
    double kmag = *param;

    vec3 k(0, 0, kmag);
    vec3 q(qmag*sqrt(1-mu_q*mu_q), 0, qmag*mu_q);
    vec3 p(pmag*sqrt(1-mu_p*mu_p)*cos(phi_p), pmag*sqrt(1-mu_p*mu_p)*sin(phi_p), pmag*mu_p);
    vec3 kq = k-q;
    vec3 kqp = k-q-p;
    double pk = P_L(kmag);
    double pq = P_L(qmag);
    double pp = P_L(pmag);
    double pkq = P_L(kq.mag());
    double pkqp = P_L(kqp.mag());
    
    f[0] = 30 * F5s(k, k, q, -q, p, -p) * pk * pq * pp;
    f[1] = 24 * F2s(k, q, kq) * F4s(-k, -q, -kq, p, -p) * pq * pp * pkq;
    f[2] = 9 * F3s(k, k, q, -q) * F3s(-k, -k, p, -p) * pk * pq * pp;
    f[3] = 6 * F3s(k, q, p, kqp) * F3s(-k, -q, -p, -kqp) * pq * pp * pkqp;
    for(int j = 0; j < 4; j++)
        f[j] *= pow3(qmag)*pow3(pmag)/pow5(2*M_PI);
}


HigherSPT::HigherSPT(const PowerSpectrum& P_L_, real qmin, real qmax, int order_)
    : P_L(P_L_), p2int(P_L_, qmin, qmax), p3int(P_L_, qmin, qmax)
{
    order = order_;
}

real HigherSPT::P1(real k) {
    return P_L(k);
}

real HigherSPT::P2(real k, real* P13, real* P22) {
    double I[2], err[2];
    double kmag = k;
    int neval;
    p2int.Integrate(I, &kmag, err, &neval, 0.001, 0.001*P_L(kmag));

    if(P13)
        *P13 = I[0];
    if(P22)
        *P22 = I[1];
    return I[0] + I[1];
}

real HigherSPT::P3(real k, real* P15, real* P24, real* P33a, real* P33b) {
    double I[4], err[4];
    double kmag = k;
    int neval;
    p3int.Integrate(I, &kmag, err, &neval, 0.001, 0.001*P_L(kmag));

    if(P15)
        *P15 = I[0];
    if(P24)
        *P24 = I[1];
    if(P33a)
        *P33a = I[2];
    if(P33b)
        *P33b = I[3];
    return I[0] + I[1] + I[2] + I[3];
}

real HigherSPT::P(real k) {
    real pk = P1(k);
    if(order >= 2)
        pk += P2(k);
    if(order >= 3)
        pk += P3(k);
    if(order >= 4)
        warning("HigherSPT only implemented up through order 3\n");
    return pk;
}
