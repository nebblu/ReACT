#if HAVE_CONFIG_H
# include <config.h>
#endif

#include "ODE.h"
#include "SPT.h"
#include "SimpleRG.h"
#include "Quadrature.h"
#include <functional>
using std::cref;
using std::bind;

SimpleRG::SimpleRG(const Cosmology& C_, const PowerSpectrum& P_L_)
    : C(C_), P_L(P_L_)
{
}

InterpolatedPS SimpleRG::CalculateP(int Nk, real kmin, real kmax, real epsrel) {
    /* Set initial conditions: Q(k, A = 0) = P_L(k) */
    k.resize(Nk);
    vector<real> Q(Nk);
    for(int i = 0; i < Nk; i++) {
        k[i] = kmin * exp(i*log(kmax/kmin)/(Nk-1));
        Q[i] = P_L(k[i]);
    }

    /* Use order adaptive Dormand-Prince solver to find Q(k, A = 1) */
    int n = RKDP(0, 1, Q, bind(&SimpleRG::dQdA, *this, std::placeholders::_1, std::placeholders::_2), Q, epsrel, 0.01);
    if(n < 0)
        error("SimpleRG: numerical integration failed\n");
    return InterpolatedPS(C, k, Q);
}

vector<real> SimpleRG::dQdA(real A, const vector<real>& Q) {
    int Nk = k.size();
    InterpolatedPS P(C, k, Q);
    SPT spt(C, P);
    vector<real> dQdA(Nk);
    int i;
    #pragma omp parallel default(shared) private(i)
    #pragma omp for schedule(dynamic)
    for(i = 0; i < Nk; i++)
        dQdA[i] = spt.P22(k[i]) + spt.P13(k[i]);
    return dQdA;
}
