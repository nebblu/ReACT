#if HAVE_CONFIG_H
# include <config.h>
#endif

#include "Kaiser.h"

#include "CorrelationFunction.h"
#include "Cosmology.h"
#include "GrowthFunction.h"
#include "LinearPS.h"
#include "Spline.h"


Kaiser::Kaiser(const Cosmology& C, const PowerSpectrum& P_, real z)
    : P(P_)
{
    GrowthFunction D(C);
    f = D.f(z);
}


Kaiser::Kaiser(const PowerSpectrum& P_, real f_)
    : P(P_), f(f_)
{
}

real Kaiser::P_s(real k, real mu) const {
    return pow2(1 + f*mu*mu) * P(k);
}

real Kaiser::PMultipole(int ell, real k) const {
    real prefactor;
    if(ell == 0)
        prefactor = 1. + (2./3.)*f + (1./5.)*f*f;
    else if(ell == 2)
        prefactor = (4./3.)*f + (4./7.)*f*f;
    else if(ell == 4)
        prefactor = (8./35.)*f*f;
    else
        prefactor = 0.;

    return prefactor * P(k);
}

array Kaiser::PMultipole(int ell, const array& k) const {
    int N = (int) k.size();
    array pell(N);
    #pragma omp parallel for
    for(int j = 0; j < N; j++)
        pell[j] = PMultipole(ell, k[j]);
    return pell;
}

real Kaiser::XiMultipole(int ell, real r, int Nk, real kmin, real kmax) const {
    double xi;
    ComputeXiLM(ell, 2, P, 1, &r, &xi, Nk, kmin, kmax);
    return xi;
}

array Kaiser::XiMultipole(int ell, const array& r, int Nk, real kmin, real kmax) const {
    int Nr = (int) r.size();
    array xi(Nr);
    ComputeXiLM(ell, 2, P, Nr, &r[0], &xi[0], Nk, kmin, kmax);
    return xi;
}

#if 0
Spline Kaiser::ComputeP_ell(int ell, real kmin, real kmax, int Nk) {
    real prefactor;
    if(ell == 0)
        prefactor = 1. + (2./3.)*f + (1./5.)*f*f;
    else if(ell == 2)
        prefactor = (4./3.)*f + (4./7.)*f*f;
    else if(ell == 4)
        prefactor = (8./35.)*f*f;
    else
        prefactor = 0.;

    array k = array::logspace(kmin, kmax, Nk);
    array pk = prefactor * P(k);
    return CubicSpline(k, pk);
}

Spline Kaiser::ComputeXi_ell(int ell, real kmin, real kmax, int N) {
    real sign = (ell % 2 == 1) ? 0. : ((ell % 4 == 0) ? +1. : -1.);

    real prefactor = sign;
    if(ell == 0)
        prefactor *= 1. + (2./3.)*f + (1./5.)*f*f;
    else if(ell == 2)
        prefactor *= (4./3.)*f + (4./7.)*f*f;
    else if(ell == 4)
        prefactor *= (8./35.)*f*f;

    array k = array::logspace(kmin, kmax, N);
    array pk = prefactor * P(k);
    array r(N), xi(N);
    ComputeXiLM(ell, 2, N, k, pk, r, xi);
    return CubicSpline(r, xi);
}

void Kaiser::ComputeXi_ell(int ell, real f, const array& k, const array& pk, array& r, array& xi) {
    int N = (int)k.size();
    real sign = (ell % 2 == 1) ? 0. : ((ell % 4 == 0) ? +1. : -1.);
    real prefactor = sign;
    if(ell == 0)
        prefactor *= 1. + (2./3.)*f + (1./5.)*f*f;
    else if(ell == 2)
        prefactor *= (4./3.)*f + (4./7.)*f*f;
    else if(ell == 4)
        prefactor *= (8./35.)*f*f;

    array pks = prefactor * pk;
    r.resize(N);
    xi.resize(N);
    ComputeXiLM(ell, 2, N, k, pks, r, xi);
}
#endif
