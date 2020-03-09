#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

//#include <gsl/gsl_qrng.h>

#include "Common.h"
#include "MonteCarlo.h"
#include "rng.h"


#if 0
struct QuasiRNG {
    gsl_qrng* q;

    QuasiRNG(int n) {
        q = gsl_qrng_alloc(gsl_qrng_sobol, n);
    }

    ~QuasiRNG() {
        gsl_qrng_free(q);
    }

    void Sample(double* x) {
        gsl_qrng_get(q, x);
    }
};
#endif


/******************************************************************************
 * MonteCarloIntegral
 ******************************************************************************/

MonteCarloIntegral::MonteCarloIntegral(int n_, int m_, const Range* domain) {
    n = n_;
    m = m_;

    /* Initialize domain */
    V.resize(n);
    SetDomain(domain);

    /* Set internal parameters */
    maxeval = DEFAULT_MAXEVAL;
    nbatch = DEFAULT_BATCHSIZE;
    nstart = DEFAULT_NSTART;
    nincrease = DEFAULT_NINCREASE;
}

MonteCarloIntegral::~MonteCarloIntegral() {
}

void MonteCarloIntegral::SetDomain(const Range* domain) {
    for(int i = 0; i < n; i++) {
        V[i].a = (domain == 0) ? 0 : domain[i].a;
        V[i].b = (domain == 0) ? 1 : domain[i].b;
    }
}

void MonteCarloIntegral::Integrate(double* I, double* param, double* errorout, int* nevalout, double epsrel, double epsabs) const {
    /* Make sure error thresholds are sensible */
    if(epsrel <= 0)
        epsrel = DEFAULT_EPSREL;
    if(epsabs <= 0)
        epsabs = DEFAULT_EPSABS;

    int neval = 0;

    /* Compute Jacobian */
    double jacobian = 1;
    for(int i = 0; i < n; i++)
        jacobian *= (V[i].b - V[i].a);

    /* Initialize state */
    State state;
    state.niter = 0;
    state.nsamples = nstart;
    state.cumul.resize(m);
    state.grid.resize(n);
    for(int i = 0; i < n; i++)
        for(int bin = 0; bin < NBINS; bin++)
            state.grid[i][bin] = (bin + 1)/(double)NBINS;

    /* Initialize quasi-random number generator */
//    QuasiRNG qrng(n);
    Sobol qrng;
    rng_sobol_init(&qrng, n);

    /* Initialize sample buffer */
    int samplebufsize = nbatch*((1 + n + m)*sizeof(double) + n*sizeof(short));
    char* samplebuf = (char*)malloc(samplebufsize);

    /* Main iteration loop */
    int notdone = 1;
    while(notdone) {
        int nsamples = state.nsamples;
        Grid margsum[m][n];
        memset(margsum, 0, sizeof(margsum));
        double base_weight = 1./nsamples;

        /* Sample integrand one batch at a time */
        while(nsamples > 0) {
            int size = (nbatch < nsamples) ? nbatch : nsamples;   // size of current batch
            double* w = (double*)samplebuf;
            double* x = w + size;
            double* f = x + size*n;
            double* lastf = f + size*m;
            short* bin = (short*)lastf;

            /* Prepare positions and weights for sampling */
            while(x < f) {
                double weight = base_weight;
//                qrng.Sample(x);
                rng_sobol_get(&qrng, x);
                for(int i = 0; i < n; i++) {
                    double pos = (*x)*NBINS;
                    int ipos = (int)pos;
                    double prev = (ipos == 0) ? 0 : state.grid[i][ipos-1];
                    double diff = state.grid[i][ipos] - prev;
                    *x++ = prev + (pos - ipos)*diff;
                    *bin++ = ipos;
                    weight *= diff*NBINS;
                }
                *w++ = weight;
            }

            /* Sample integrand (why does Cuba's DoSample() pass w to the integrand?) */
            f = x;
            x = w;
            double y[n];
            for(int s = 0; s < size; s++) {
                for(int i = 0; i < n; i++)
                    y[i] = V[i].a + x[i]*(V[i].b - V[i].a);
                Integrand(y, f, param);
                x += n;
                f += m;
            }
            neval += size;

            /* Adjust weights */
            f = x;
            w = (double*)samplebuf;
            bin = (short*)lastf;
            while(f < lastf) {
                double weight = *w++;
                for(int j = 0; j < m; j++) {
                    double wfun = weight * (*f++);
                    if(wfun != 0) {
                        Cumulants* c = &state.cumul[j];
                        Grid* ms = margsum[j];
                        c->sum += wfun;
                        c->sqsum += wfun*wfun;
                        for(int i = 0; i < n; i++)
                            ms[i][bin[i]] += wfun*wfun;
                    }
                }
                bin += n;
            }

            nsamples -= nbatch;
        }

        notdone = 0;

        /* Compute the integral and error values */
        for(int j = 0; j < m; j++) {
            Cumulants* c = &state.cumul[j];
            double w = Weight(c->sum, c->sqsum, state.nsamples);

            c->weightsum += w;
            c->avgsum += w*c->sum;
            double sigsq = 1/c->weightsum;
            c->avg = sigsq*c->avgsum;
            c->err = sqrt(sigsq);
            notdone |= (c->err > fmax(epsrel*fabs(c->avg), epsabs/jacobian));

            if(state.niter == 0)
                c->guess = c->sum;
            else {
                w *= (c->sum - c->guess);
                c->chisum += w;
                c->chisqsum += w*c->sum;
            }
            c->chisq = c->chisqsum - c->avg*c->chisum;

            c->sum = c->sqsum = 0;
        }

//#ifdef VERBOSE
//        {
//            char s[128 + 128*m], *p = s;
//            p += sprintf(p, "Iteration %d: %d integrand evaluations\n", state.niter+1, neval);
//            for(int j = 0; j < m; j++) {
//                Cumulants *c = &state.cumul[j];
//                p += sprintf(p, "  I_%d =  %g +/- %g  \tchisq %g (%d df)\n", j, jacobian*c->avg, jacobian*c->err, c->chisq, state.niter);
//            }
//            verbose("%s\n", s);
//        }
//#endif

        /* Finish if we're below error thresholds */
        if(notdone == 0)
            break;

        /* Abort (with a warning) if we've reached the maximum number of evaluations */
        if(neval >= maxeval) {
            warning("MonteCarlo: maximum number of evaluations reached\n");
            break;
        }

        /* Refine the grid for the next iteration */
        if(m == 1) {
            for(int i = 0; i < n; i++)
                RefineGrid(state.grid[i], margsum[0][i]);
        }
        else {
            for(int i = 0; i < n; i++) {
                Grid wmargsum;
                memset(wmargsum, 0, sizeof(wmargsum));
                for(int j = 0; j < m; j++) {
                    double w = state.cumul[j].avg;
                    if(w != 0) {
                        double* ms = margsum[j][i];
                        for(int bin = 0; bin < NBINS; bin++)
                            wmargsum[bin] += ms[bin]/(w*w);
                    }
                }
                RefineGrid(state.grid[i], wmargsum);
            }
        }

        /* Proceed to the next iteration */
        state.niter++;
        state.nsamples += nincrease;
    }

    /* Prepare results */
    for(int j = 0; j < m; j++) {
        Cumulants* c = &state.cumul[j];
        I[j] = jacobian * c->avg;
        if(errorout)
            errorout[j] = jacobian * c->err;
    }
    if(nevalout)
        *nevalout = neval;
}


double MonteCarloIntegral::Weight(double sum, double sqsum, int nsamples) {
    double w = sqrt(sqsum*nsamples);
    return (nsamples - 1)/fmax((w + sum)*(w - sum), DBL_MIN);
}

void MonteCarloIntegral::RefineGrid(Grid& grid, Grid& margsum) {
    double avgperbin, thisbin, newcur, delta;
    Grid imp, newgrid;
    int bin, newbin;

    /* Smooth the f^2 value stored for each bin */
    double prev = margsum[0];
    double cur = margsum[1];
    double norm = margsum[0] = 0.5*(prev + cur);
    for(bin = 1; bin < NBINS-1; bin++) {
        double s = prev + cur;
        prev = cur;
        cur = margsum[bin+1];
        norm += margsum[bin] = (s + cur)/3;
    }
    norm += margsum[NBINS-1] = 0.5*(prev + cur);

    if(norm == 0)
        return;
    norm = 1/norm;

    /* Compute the importance function for each bin */
    avgperbin = 0;
    for(bin = 0; bin < NBINS; bin++) {
        double impfun = 0;
        if(margsum[bin] > 0) {
            double r = margsum[bin]*norm;
            avgperbin += impfun = pow((r-1)/log(r), 1.5);
        }
        imp[bin] = impfun;
    }
    avgperbin /= NBINS;

    /* Redefine the current size of each bin */
    cur = newcur = 0;
    thisbin = 0;
    bin = -1;
    for(newbin = 0; newbin < NBINS-1; newbin++) {
        while(thisbin < avgperbin) {
            thisbin += imp[++bin];
            prev = cur;
            cur = grid[bin];
        }
        thisbin -= avgperbin;
        delta = (cur - prev)*thisbin;
        newgrid[newbin] = newcur = fmax(newcur, cur - 2*delta/(imp[bin] + imp[bin == 0 ? 0 : bin-1]));
    }

//    memcpy(grid, newgrid, (NBINS-1)*sizeof(double));
    for(bin = 0; bin < NBINS-1; bin++)
        grid[bin] = newgrid[bin];
    grid[NBINS-1] = 1;
}
