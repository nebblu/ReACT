#ifndef LINEAR_PS_H
#define LINEAR_PS_H

#include "Common.h"
#include "NoWigglePS.h"
#include "PowerSpectrum.h"
#include "Spline.h"


/**********************************************************************
 * LinearPS
 *
 * Linear power spectrum interpolated from a transfer function, which is
 * presumably calculated from a Boltzmann code.  Uses the analytic
 * $\ln^2(k)/k^3$ formula at high $k$.
 *********************************************************************/
class LinearPS : public PowerSpectrum {
public:
    LinearPS(const Cosmology& C, real z = 0);

    const Cosmology& GetCosmology() const { return C; }
    real Evaluate(real k) const;

private:
//    const Cosmology& C;
    Cosmology C;
    real z;
    real k0, p0;        // k and P(k) of left-most data point
    real k1, p1;        // k and P(k) of right-most data point
    Spline pk;          // spline P(k) based on transfer function
    NoWigglePS Pnw;     // no-wiggle P(k) for high k
};


// Takes as input the matter power spectrum
class myLinearPS :
public PowerSpectrum {
public:
    myLinearPS(const Cosmology& C, real z = 0);

    const Cosmology& GetCosmology() const { return C; }
    real Evaluate(real k) const;

private:
//    const Cosmology& C;
    Cosmology C;
    real z;
    real k0, p0;        // k and P(k) of left-most data point
    real k1, p1;        // k and P(k) of right-most data point
    Spline pk;          // spline P(k) based on transfer function
    NoWigglePS Pnw;     // no-wiggle P(k) for high k
};



// Takes as input transfer function at z and A_s (not sigma8)
class LinearPS_as :
public PowerSpectrum {
public:
    LinearPS_as(const Cosmology& C, real z = 0);

    const Cosmology& GetCosmology() const { return C; }
    real Evaluate(real k) const;

private:
//    const Cosmology& C;
    Cosmology C;
    real z;
    real k0, p0;        // k and P(k) of left-most data point
    real k1, p1;        // k and P(k) of right-most data point
    Spline pk;          // spline P(k) based on transfer function
    NoWigglePS Pnw;     // no-wiggle P(k) for high k
};



#endif // LINEAR_PS_H
