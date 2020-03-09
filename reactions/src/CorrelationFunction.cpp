#if HAVE_CONFIG_H
# include <config.h>
#endif
#include "CorrelationFunction.h"
#include "PowerSpectrum.h"
#include "Quadrature.h"
#include "SpecialFunctions.h"
#include "array.h"
#include "Spline.h"
#include "SPT.h"

#include <cuba.h>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cerrno>
#include <cstdlib>
#include <stdlib.h>
#include <sys/stat.h>
#include <math.h>
#include <vector>
#include <functional>

using std::cref;
using std::bind;


//#include<omp.h>

// LOAD LINEAR POWER SPECTRUM ///
CorrelationFunction::CorrelationFunction(const Cosmology& C, const PowerSpectrum& P, real kmin_, real kmax_)
    :C(C), P(P), kmin(kmin_), kmax(kmax_)
{ }

void ComputeXiLM(int l, int m, const PowerSpectrum& P,
                 int Nr, const double r[], double xi[],
                 int Nk, double kmin, double kmax)
 {
    assert(Nk > 0 && (Nk % 2) == 0);
    const double dk = (kmax - kmin)/Nk;

    /* Choose appropriate spherical Bessel function */
    double (*sj)(double x);
    if(l == 0)      sj = SphericalBesselJ0;
    else if(l == 1) sj = SphericalBesselJ1;
    else if(l == 2) sj = SphericalBesselJ2;
    else if(l == 3) sj = SphericalBesselJ3;
    else if(l == 4) sj = SphericalBesselJ4;
    else {
        fprintf(stderr, "ComputeXiLM: l = %d not supported\n", l);
        return;
    }

    array k = array::linspace(kmin, kmax, Nk+1);
    array mult(Nk+1);

    #pragma omp parallel for
    for(int j = 0; j <= Nk; j++) {
        /* Multiplicative factor for Simpson's rule: either 1, 2, or 4 */
        mult[j] = 2 + 2*(j % 2) - (j == 0) - (j == Nk);
        /* All other purely k-dependent factors */
        mult[j] *= P(k[j]) * pow(k[j], m) * (dk/3) / (2*M_PI*M_PI);
    }

    /* Integrate $P(k) k^m j_l(kr) dk$ over the interval $[kmin,kmax]$ using Simpson's rule */
    #pragma omp parallel for
    for(int i = 0; i < Nr; i++) {
        xi[i] = 0;
        for(int j = 0; j <= Nk; j++)
            xi[i] += mult[j] * sj(k[j]*r[i]);
    }
}
static real f(const PowerSpectrum& P, real r, real k) {
    return k*sin(k*r)*P(k);
}
real CorrelationFunction::Evaluate(real r) const {
    return 1./(2*M_PI*M_PI*r) * Integrate(bind(f, cref(P), r, std::placeholders::_1), kmin, kmax);
}
array CorrelationFunction::EvaluateMany(const array& r) const {
    int Nr = (int) r.size();
    array xi(Nr);
    #pragma omp parallel for
    for(int i = 0; i < Nr; i++)
        xi[i] = Evaluate(r[i]);

    return xi;
}
