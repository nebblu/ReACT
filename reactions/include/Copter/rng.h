#ifndef RNG_H
#define RNG_H

#ifdef __cplusplus
extern "C" {
#endif

/* Mersenne Twister pseudo-random number generator.  Based on the
 * implementation of src/common/Random.c of Cuba 1.4, itself an adaptation of
 * the original C code of T. Nishimura and M. Matsumoto from 
 *   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html */

/* 32 or 53 random bits */
#ifndef MERSENNE_NBITS
#define MERSENNE_NBITS 32
#endif

/* Length of state vector */
#define MERSENNE_N 624

/* Period parameter */
#define MERSENNE_M 397

/* Standard interface.  Uses static state vector.  Not thread-safe. */
void rng_init(unsigned int seed);
unsigned int rng_random();
double rng_uniform();   // return value lies in half-open interval [0,1)

typedef struct {
    unsigned int state[MERSENNE_N];
    int next;
} rng_state;

/* Thread-safe variants.  User allocates a 'rng_state' variable, initializes it
 * with a call to 'rng_init_r()', then draws randoms from it with either
 * 'rng_random_r()' or 'rng_uniform_r()'. */
void rng_init_r(rng_state *rng, unsigned int seed);
unsigned int rng_random_r(rng_state *rng);
double rng_uniform_r(rng_state *rng);



/* Sampling from non-uniform distributions. */

/* Standard normal distribution */
double rng_normal();
double rng_normal_r(rng_state *rng);

/* Exponential distribution */
double rng_exponential(double lambda);
double rng_exponential_r(rng_state *rng, double lambda);

/* Poisson distribution (note that the return value is integral, even though
 * it is returned as a double; cast to int as desired) */
double rng_poisson(double mu);
double rng_poisson_r(rng_state *rng, double mu);


/* Sobol quasi-random sequence generator.  Based on ACM TOMS algorithm 659,
 * adapted from src/common/Random.c of Cuba 1.4.
 *
 * Here's a usage example showing how to perform quasi-Monte Carlo integration
 * of the function f(x,y,z) over the unit cube:
 *
 *   int i, npoints = 1000000;
 *   double x[3], sum;
 *   Sobol seq;
 *   rng_sobol_init(&seq, 3);
 *   for(i = 0; i < npoints; i++) {
 *       rng_sobol_get(&seq, x);
 *       sum += f(x[0], x[1], x[2]);
 *   }
 *   sum /= npoints;
 */

/* Maximum dimension of Sobol sample space */
#ifndef SOBOL_MAXDIM
#define SOBOL_MAXDIM 8
#endif

typedef struct {
    int ndim;
    double norm;
    int v[SOBOL_MAXDIM][30], prev[SOBOL_MAXDIM];
    int seq;
} Sobol;

void rng_sobol_init(Sobol *sobol, int ndim);
void rng_sobol_get(Sobol *sobol, double *x);
void rng_sobol_skip(Sobol *sobol, int n);

#ifdef __cplusplus
}
#endif

#endif // RNG_H
