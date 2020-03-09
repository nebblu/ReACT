#ifndef CORRELATION_FUNCTION_H
#define CORRELATION_FUNCTION_H
#include "Common.h"
#include "array.h"
#include "SPT.h"
#include "Cosmology.h"


/* Compute the quantity
 *   \xi_l^m(r) = \int_0^\infty \frac{dk}{2\pi^2} k^m j_l(kr) P(k)
 * using Simpson's rule.  The parameters Nk, kmin, and kmax determine the
 * accuracy of the integration.  Note that Nk must be even. */
void ComputeXiLM(int l, int m, const PowerSpectrum& P,
                 int Nr, const double r[], double xi[],
                 int Nk = 32768, double kmin = 0., double kmax = 30.);


class CorrelationFunction {
public:
    CorrelationFunction(const Cosmology& C, const PowerSpectrum& P, real kmin = 1e-3, real kmax = 30.);

    real Evaluate(real r) const;
    real operator()(real r) const { return Evaluate(r); }

    array EvaluateMany(const array& r) const;
    array operator()(const array& r) const { return EvaluateMany(r); }

    const PowerSpectrum& GetPowerSpectrum() const { return P; }

protected:
    const PowerSpectrum& P;
    const Cosmology& C;
    real kmin, kmax;
};

#endif // CORRELATION_FUNCTION_H
