#include <assert.h>
#include <math.h>
#include <string.h>

#include "rng.h"

/******************************************************************************
 * Mersenne Twister pseudo-random number generator
 ******************************************************************************/

/* Note that we're assuming here that ints are 32-bit.  This should be safe on
 * almost all target systems. */

static rng_state default_rng;

static inline unsigned int Twist(unsigned int a, unsigned int b) {
    unsigned int mixbits = (a & 0x80000000) | (b & 0x7fffffff);
    unsigned int matrixA = (-(b & 1)) & 0x9908b0df;
    return (mixbits >> 1) ^ matrixA;
}

static inline void MersenneReload(rng_state *rng) {
    unsigned int *s = rng->state;
    int j;

    for(j = MERSENNE_N - MERSENNE_M + 1; --j; ++s)
        *s = s[MERSENNE_M] ^ Twist(s[0], s[1]);
    for(j = MERSENNE_M; --j; ++s)
        *s = s[MERSENNE_M - MERSENNE_N] ^ Twist(s[0], s[1]);
    *s = s[MERSENNE_M - MERSENNE_N] ^ Twist(s[0], rng->state[0]);
    rng->next = 0;
}

void rng_init_r(rng_state *rng, unsigned int seed) {
    unsigned int *s = rng->state;
    int j;

    for(j = 1; j <= MERSENNE_N; ++j) {
        *s++ = seed;
        seed = 0x6c078965*(seed ^ (seed >> 30)) + j;
        /* see Knuth TAOCP Vol 2, 3rd Ed, p. 106 for multiplier */
    }
    MersenneReload(rng);
}

void rng_init(unsigned int seed) {
    rng_init_r(&default_rng, seed);
}

unsigned int rng_random_r(rng_state *rng) {
    unsigned int s;
    if(rng->next >= MERSENNE_N)
        MersenneReload(rng);
    s = rng->state[rng->next++];
    s ^= s >> 11;
    s ^= (s << 7) & 0x9d2c5680;
    s ^= (s << 15) & 0xefc60000;
    return s ^ (s >> 18);
}

unsigned int rng_random() {
    return rng_random_r(&default_rng);
}

double rng_uniform_r(rng_state *rng) {
#if MERSENNE_NBITS == 53
    unsigned int a = rng_random_r(rng) >> 5;
    unsigned int b = rng_random_r(rng) >> 6;
    return (67108864.*a + b)/9007199254740992.;
#else
    return rng_random_r(rng)/4294967295.;
#endif
}

double rng_uniform() {
    return rng_uniform_r(&default_rng);
}


/******************************************************************************
 * Non-uniform sampling
 ******************************************************************************/

#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI 0.3989422804014327
#endif

/* Compute the quantile function for the standard normal distribution, using
 * the PPND16 routine of Wichura (1988). */
static double qnorm(double p) {
    static double a[8] = { 3.3871328727963666080E0,
                           1.3314166789178437745E2,
                           1.9715909503065514427E3,
                           1.3731693765509461125E4,
                           4.5921953931549871457E4,
                           6.7265770927008700853E4,
                           3.3430575583588128105E4,
                           2.5090809287301226727E3 };
    static double b[8] = { 1,
                           4.2313330701600911252E1,
                           6.8718700749205790830E2,
                           5.3941960214247511077E3,
                           2.1213794301586595867E4,
                           3.9307895800092710610E4,
                           2.8729085735721942674E4,
                           5.2264952788528545610E3 };
    static double c[8] = { 1.42343711074968357734E0,
                           4.63033784615654529590E0,
                           5.76949722146069140550E0,
                           3.64784832476320460504E0,
                           1.27045825245236838258E0,
                           2.41780725177450611770E-1,
                           2.27238449892691845833E-2,
                           7.74545014278341407640E-4 };
    static double d[8] = { 1,
                           2.05319162663775882187E0,
                           1.67638483018380384940E0,
                           6.89767334985100004550E-1,
                           1.48103976427480074590E-1,
                           1.51986665636164571966E-2,
                           5.47593808499534494600E-4,
                           1.05075007164441684324E-9 };
    static double e[8] = { 6.65790464350110377720E0,
                           5.46378491116411436990E0,
                           1.78482653991729133580E0,
                           2.96560571828504891230E-1,
                           2.65321895265761230930E-2,
                           1.24266094738807843860E-3,
                           2.71155556874348757815E-5,
                           2.01033439929228813265E-7 };
    static double f[8] = { 1,
                           5.99832206555887937690E-1,
                           1.36929880922735805310E-1,
                           1.48753612908506148525E-2,
                           7.86869131145613259100E-4,
                           1.84631831751005468180E-5,
                           1.42151175831644588870E-7,
                           2.04426310338993978564E-15 };
    double q, r, s;

    q = p - 0.5;
    s = (q > 0) - (q < 0);      /* s = sgn(q) */
    if(q*s <= 0.425) {          /* 0.075 <= p <= 0.925 */
        r = 0.180625 - q*q;
        return q * (((((((a[7]*r + a[6])*r + a[5])*r + a[4])*r + a[3])*r + a[2])*r + a[1])*r + a[0])
                 / (((((((b[7]*r + b[6])*r + b[5])*r + b[4])*r + b[3])*r + b[2])*r + b[1])*r + 1);
    }
    else {                      /* p < 0.075 or p > 0.925 */
        r = (q < 0) ? p : 1-p;
        r = sqrt(-log(r));
        if(r <= 5.) {
            r -= 1.6;
            return s * (((((((c[7]*r + c[6])*r + c[5])*r + c[4])*r + c[3])*r + c[2])*r + c[1])*r + c[0])
                     / (((((((d[7]*r + d[6])*r + d[5])*r + d[4])*r + d[3])*r + d[2])*r + d[1])*r + 1);
        }
        else {
            r -= 5.;
            return s * (((((((e[7]*r + e[6])*r + e[5])*r + e[4])*r + e[3])*r + e[2])*r + e[1])*r + e[0])
                     / (((((((f[7]*r + f[6])*r + f[5])*r + f[4])*r + f[3])*r + f[2])*r + f[1])*r + 1);
        }
    }
}

double rng_normal_r(rng_state *rng) {
    double u = rng_uniform_r(rng);
    return qnorm(u);
}

double rng_normal() {
    return rng_normal_r(&default_rng);
}

double rng_exponential_r(rng_state *rng, double lambda) {
    double u = rng_uniform_r(rng);
    return -log(1-u)/lambda;
}

double rng_exponential(double lambda) {
    return rng_exponential_r(&default_rng, lambda);
}

static void ProcedureF(double k, double mu, double s, double *px, double *py, double *fx, double *fy) {
    /* k! for 0 <= k <= 9 */
    static double fact[10] = { 1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880 };
    /* Table I coefficients */
    static double a[] = { -0.5000000002,
                           0.3333333343,
                          -0.2499998565,
                           0.1999997049,
                          -0.1666848753,
                           0.1428833286,
                          -0.1241963125,
                           0.1101687109,
                          -0.1142650302,
                           0.1055093006 };

    double omega, b1, b2, c0, c1, c2, c3, X;
    double mu_k = mu - k;
    if(k < 10) {
        *px = -mu;
        *py = pow(mu, k) / fact[(int)k];
    }
    else {
        double delta, v;
        delta = 1/(12.*k);
        delta = delta*(1 - 4.8*delta*delta);
        v = mu_k/k;
        if(fabs(v) <= 0.25)
            *px = mu_k*v*(((((((((a[9]*v + a[8])*v + a[7])*v + a[6])*v + a[5])*v + a[4])*v + a[3])*v + a[2])*v + a[1])*v + a[0]) - delta;
        else
            *px = k*log(1 + v) - mu_k - delta;
        *py = M_1_SQRT_2PI / sqrt(k);
    }

    omega = M_1_SQRT_2PI / s;
    b1 = 1/(24.*mu);
    b2 = 0.3*b1*b1;
    c3 = b1*b2/7.;
    c2 = b2 - 15.*c3;
    c1 = b1 - 6.*b2 + 45.*c3;
    c0 = 1. - b1 + 3.*b2 - 15.*c3;

    X = (0.5 - mu_k)/s;
    X *= X;
    *fx = -0.5*X;
    *fy = omega*(((c3*X + c2)*X + c1)*X + c0);
}

static double sgn(double x) {
    return (x > 0) - (x < 0);
}

/* A modified implementation of the algorithm of Ahrens and Dieter, "Computer
 * Generation of Poisson Deviates from Modified Normal Distributions" (1982).
 * Here we don't bother retaining state between consecutive calls with the same
 * value of mu, since this is not a usage case we imagine being important. */
double rng_poisson_r(rng_state *rng, double mu) {
    double k;   /* resultant Poisson deviate */

    if(mu >= 10.) {
        /* Case A */
        double E, G, T, u, mu_k, c, fx, fy, px, py;
        double s = sqrt(mu);

        /* N: Normal sample */
        T = rng_normal_r(rng);
        G = mu + s*T;
        if(G >= 0) {
            k = floor(G);       // k is integral
            if(k >= floor(mu - 1.1484)) {
                /* I: Immediate acceptance */
                return k;
            }
            /* S: Squeeze acceptance */
            u = rng_uniform_r(rng);
            mu_k = mu - k;
            if(6*mu*mu*u >= mu_k*mu_k*mu_k)
                return k;
        }

        /* P: Preparations for Q and H (mostly relegated to ProcedureF) */
        c = 0.1069/mu;

        if(G >= 0) {
            ProcedureF(k, mu, s, &px, &py, &fx, &fy);
            /* Q: Quotient acceptance */
            if(fy*(1 - u) <= py*exp(px - fx))
                return k;
        }

        /* E: Double exponential sample */
        for(;;) {
            do {
                E = rng_exponential_r(rng, 1);
                u = rng_uniform_r(rng);
                u = u + u - 1;
                T = 1.8 + E*sgn(u);
            } while(T <= -0.6744);
            k = floor(mu + s*T);
            ProcedureF(k, mu, s, &px, &py, &fx, &fy);

            /* H: Hat acceptance */
            if(c*fabs(u) <= py*exp(px + E) - fy*exp(fx + E))
                return k;
        }

    }
    else {
        /* Case B: mu < 10 */
        double p, q, u;
        for(;;) {       /* This loop should almost never exceed one iteration: P[k > 35] < 2e-10 for all mu < 10 */
            p = q = exp(-mu);
            u = rng_uniform_r(rng);

            for(k = 0; k <= 35; ) {
                if(u <= q)
                    return k;
                k += 1;
                p *= mu/k;
                q += p;
            }
        }
    }
}

double rng_poisson(double lambda) {
    return rng_poisson_r(&default_rng, lambda);
}


/******************************************************************************
 * Sobol quasi-random sequence generator
 ******************************************************************************/

void rng_sobol_init(Sobol *sobol, int ndim) {
    static int neval = (1 << 30) - 1;   // use all available bits
    static int ini[9*40] = {
          3,   1,   0,   0,   0,   0,   0,   0,   0,
          7,   1,   1,   0,   0,   0,   0,   0,   0,
         11,   1,   3,   7,   0,   0,   0,   0,   0,
         13,   1,   1,   5,   0,   0,   0,   0,   0,
         19,   1,   3,   1,   1,   0,   0,   0,   0,
         25,   1,   1,   3,   7,   0,   0,   0,   0,
         37,   1,   3,   3,   9,   9,   0,   0,   0,
         59,   1,   3,   7,  13,   3,   0,   0,   0,
         47,   1,   1,   5,  11,  27,   0,   0,   0,
         61,   1,   3,   5,   1,  15,   0,   0,   0,
         55,   1,   1,   7,   3,  29,   0,   0,   0,
         41,   1,   3,   7,   7,  21,   0,   0,   0,
         67,   1,   1,   1,   9,  23,  37,   0,   0,
         97,   1,   3,   3,   5,  19,  33,   0,   0,
         91,   1,   1,   3,  13,  11,   7,   0,   0,
        109,   1,   1,   7,  13,  25,   5,   0,   0,
        103,   1,   3,   5,  11,   7,  11,   0,   0,
        115,   1,   1,   1,   3,  13,  39,   0,   0,
        131,   1,   3,   1,  15,  17,  63,  13,   0,
        193,   1,   1,   5,   5,   1,  27,  33,   0,
        137,   1,   3,   3,   3,  25,  17, 115,   0,
        145,   1,   1,   3,  15,  29,  15,  41,   0,
        143,   1,   3,   1,   7,   3,  23,  79,   0,
        241,   1,   3,   7,   9,  31,  29,  17,   0,
        157,   1,   1,   5,  13,  11,   3,  29,   0,
        185,   1,   3,   1,   9,   5,  21, 119,   0,
        167,   1,   1,   3,   1,  23,  13,  75,   0,
        229,   1,   3,   3,  11,  27,  31,  73,   0,
        171,   1,   1,   7,   7,  19,  25, 105,   0,
        213,   1,   3,   5,   5,  21,   9,   7,   0,
        191,   1,   1,   1,  15,   5,  49,  59,   0,
        253,   1,   1,   1,   1,   1,  33,  65,   0,
        203,   1,   3,   5,  15,  17,  19,  21,   0,
        211,   1,   1,   7,  11,  13,  29,   3,   0,
        239,   1,   3,   7,   5,   7,  11, 113,   0,
        247,   1,   1,   5,   3,  15,  19,  61,   0,
        285,   1,   3,   1,   1,   9,  27,  89,   7,
        369,   1,   1,   3,   7,  31,  15,  45,  23,
        299,   1,   3,   3,   9,   9,  25, 107,  39
    };

    int dim, bit, nbits;
    int max, *pini = ini;

    assert(ndim >= 1 && ndim <= SOBOL_MAXDIM);
    sobol->ndim = ndim;
    for(nbits = 0, max = 1; max <= neval; max <<= 1)
        ++nbits;
    sobol->norm = 1./max;

    for(bit = 0; bit < nbits; ++bit)
        sobol->v[0][bit] = (max >>= 1);

    for(dim = 1; dim < ndim; ++dim) {
        int *pv = sobol->v[dim], *pvv = pv;
        int powers = *pini++, j;
        int inibits = -1, bit;
        for(j = powers; j; j >>= 1)
            ++inibits;

        memcpy(pv, pini, inibits*sizeof(*pini));
        pini += 8;

        for(bit = inibits; bit < nbits; ++bit) {
            int newv = *pvv, j = powers;
            int b;
            for(b = 0; b < inibits; ++b) {
                if(j & 1)
                    newv ^= pvv[b] << (inibits - b);
                j >>= 1;
            }
            pvv[inibits] = newv;
            ++pvv;
        }

        for(bit = 0; bit < nbits - 1; ++bit)
            pv[bit] <<= nbits - bit - 1;
    }

    sobol->seq = 0;
    memset(sobol->prev, 0, ndim*sizeof(int));
}

void rng_sobol_get(Sobol *sobol, double *x) {
    int seq = sobol->seq++;
    int zerobit = 0, dim;

    while(seq & 1) {
        ++zerobit;
        seq >>= 1;
    }

    for(dim = 0; dim < sobol->ndim; ++dim) {
        sobol->prev[dim] ^= sobol->v[dim][zerobit];
        x[dim] = sobol->prev[dim]*sobol->norm;
    }
}

void rng_sobol_skip(Sobol *sobol, int n) {
    assert(n >= 0);
    while(n--) {
        int seq = sobol->seq++;
        int zerobit = 0, dim;

        while(seq & 1) {
            ++zerobit;
            seq >>= 1;
        }

        for(dim = 0; dim < sobol->ndim; ++dim)
            sobol->prev[dim] ^= sobol->v[dim][zerobit];
    }
}
