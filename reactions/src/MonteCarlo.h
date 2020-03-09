#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include <vector>

#define NBINS 128
#define DEFAULT_MAXEVAL 100000
#define DEFAULT_BATCHSIZE 1000
#define DEFAULT_EPSREL 1e-3
#define DEFAULT_EPSABS 1e-12
#define DEFAULT_NSTART 1000
#define DEFAULT_NINCREASE 1000


/* MonteCarloIntegral
 *
 * n-dimensional Monte-Carlo integration, based on the implementation of the
 * VEGAS importance-sampling algorithm found in the Cuba library.  Uses GSL for
 * quasi-random number generation. */

class MonteCarloIntegral {
public:
    struct Range {
        double a;
        double b;
    };

    /* Initialize integral.
     *   n is the dimension of the integration region
     *   m is the number of components of the integrand (m > 1 for vector integrands)
     *   V is an array of n Range objects specifying the domain of integration */
    MonteCarloIntegral(int n, int m = 1, const Range* V = 0);

    virtual ~MonteCarloIntegral();

    /* Set domain of integration */
    void SetDomain(const Range* domain);

    /* Evaluate the integrand f(x).  Here x is an n-dimensional vector, and f
     * is an m-dimensional array. */
    virtual void Integrand(const double x[], double* f, double* param) const = 0;

    /* Evaluate the integral $I = \int f(x) d^nx$ */
    virtual void Integrate(double* I, double* param = 0, double* error = 0, int* neval = 0, double epsrel = 0, double epsabs = 0) const;

public:
    int n;                      // dimension of integral
    int m;                      // number of components of f(x)
    std::vector<Range> V;       // domain of integration

    /* Parameters */
    int maxeval;                // maximum number of evaluations of integrand
    int nstart;                 // number of samples in first iteration
    int nincrease;              // number of extra samples in each successive iteration
    int nbatch;                 // number of samples per batch (so everything fits in CPU cache)

protected:
    struct Cumulants {
        double sum;
        double sqsum;
        double weightsum;
        double avgsum;
        double chisum;
        double chisqsum;
        double guess;
        double avg;
        double err;
        double chisq;
    };

    struct Grid {
        double g[NBINS];
        double& operator[](int k) { return g[k]; }
        operator double*() { return &g[0]; }
    };

    struct State {
        int niter;
        int nsamples;
        int neval;
        std::vector<Cumulants> cumul;
        std::vector<Grid> grid;
    };

    static double Weight(double sum, double sqsum, int nsamples);
    static void RefineGrid(Grid& grid, Grid& margsum);
};

#if 0
typedef void (*Integrand)(const double x[], double* f, double* param);

template<int n, int m = 1>
void Integrate(Integrand f, const Range[] domain, double* F, double* param = 0, double* error = 0, int* neval = 0) {
    MonteCarloIntegral(n, m, domain);
}
#endif

#endif // MONTE_CARLO_H
