/* HigherSPT
 *
 * The SPT power spectrum to 2 loops and beyond. */

#ifndef HIGHER_SPT_H
#define HIGHER_SPT_H

#include "Common.h"
#include "MonteCarlo.h"

class HigherSPT {
public:
    HigherSPT(const PowerSpectrum& P_L, real qmin = 1e-4, real qmax = 1e2, int order = 3);

    real P1(real k);
    real P2(real k, real* P13 = 0, real* P22 = 0);
    real P3(real k, real* P15 = 0, real* P24 = 0, real* P33a = 0, real* P33b = 0);

    /* $P(k) = \sum_{n=1}^{order} P^{(n)}(k)$ */
    real P(real k);

protected:
    const PowerSpectrum& P_L;
    int order;

    struct P2_Integral : public MonteCarloIntegral {
        P2_Integral(const PowerSpectrum& P_L, real qmin, real qmax);
        void Integrand(const double x[], double* f, double* param) const;

        const PowerSpectrum& P_L;
    } p2int;

    struct P3_Integral : public MonteCarloIntegral {
        P3_Integral(const PowerSpectrum& P_L, real qmin, real qmax);
        void Integrand(const double x[], double* f, double* param) const;

        const PowerSpectrum& P_L;
    } p3int;
};

#endif // HIGHER_SPT_H
