#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <cassert>
#include <cstdio>
#include <ctime>
#include <functional>

#include "Quadrature.h"
#include "PowerSpectrum.h"

using std::cref;
using std::bind;


PowerSpectrum::PowerSpectrum() {
}

PowerSpectrum::~PowerSpectrum() {
}

array PowerSpectrum::EvaluateMany(const array& k) const {
    int n = (int)k.size();
    array pk(n);
    #pragma omp parallel for
    for(int i = 0; i < n; i++)
        pk[i] = Evaluate(k[i]);
    return pk;
}

/* Top-hat window function */
static real W(real x) {
    if(x < 1e-5)
        return 1 - (1./30.)*x*x;
    else
        return 3/pow3(x) * (sin(x) - x*cos(x));
}

static real f(const PowerSpectrum& P, real R, real k) {
    return k*k/(2*M_PI*M_PI) * P(k) * pow2(W(k*R));
}
real PowerSpectrum::Sigma(real R) const {
    real sigma2 = Integrate<ExpSub>(bind(f, cref(*this), R, std::placeholders::_1), 1e-5, 1e2, 1e-5, 1e-12);
    return sqrt(sigma2);
}

static real g(const PowerSpectrum& P, real k) {
    return P(k);
}
real PowerSpectrum::VelocityDispersion() const {
    return 1/(6*M_PI*M_PI) * Integrate<ExpSub>(bind(g, cref(*this), std::placeholders::_1), 1e-5, 1e2, 1e-5, 1e-12);
}

real PowerSpectrum::NonlinearScale() const {
    return 1/sqrt(VelocityDispersion());
}

void PowerSpectrum::Save(const char* filename, real kmin, real kmax, int Nk, bool logscale) {
    FILE* fp = fopen(filename, "w+");
    if(!fp) {
        fprintf(stderr, "PowerSpectrum::Save(): could not write to %s\n", filename);
        return;
    }

    time_t t = time(0);
    fprintf(fp, "# Power spectrum saved at %s\n", ctime(&t));
    fprintf(fp, "# k -- P(k)\n");
    for(int i = 0; i < Nk; i++) {
        real k = logscale ? kmin*exp(i*log(kmax/kmin)/(Nk-1)) : kmin + i*(kmax - kmin)/(Nk-1);
        real p = Evaluate(k);
        fprintf(fp, "%e %e\n", k, p);
    }
}


#if 0
RedshiftPowerSpectrum::RedshiftPowerSpectrum() {
}

RedshiftPowerSpectrum::~RedshiftPowerSpectrum() {
}
#endif
