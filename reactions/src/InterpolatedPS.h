#ifndef INTERPOLATED_PS_H
#define INTERPOLATED_PS_H

#include "Cosmology.h"
#include "PowerSpectrum.h"
#include "Spline.h"


/**********************************************************************
 * InterpolatedPS
 *
 * Power spectrum interpolated a set of data points.
 *********************************************************************/

class InterpolatedPS : public PowerSpectrum {
public:
    /* Empty constructor */
    InterpolatedPS();

    /* Guess whether points are logarithmic or linearly spaced */
    InterpolatedPS(const Cosmology& C, const char* filename, int kcol = 1, int pcol = 2);
    InterpolatedPS(const Cosmology& C, const array& k, const array& p);
    InterpolatedPS(const Cosmology& C, int N, const real* k, const real* p);

    /* Specify logarithmic or linear spacing explicitly */
    InterpolatedPS(const Cosmology& C, const array& k, const array& p, bool logscale);
    InterpolatedPS(const Cosmology& C, int N, const real* k, const real* p, bool logscale);

    real Evaluate(real k) const;
    const Cosmology& GetCosmology() const;

private:
    const Cosmology* C;
    real n;             // spectral index
    bool logscale;      // guess about whether the k-values are logarithmic
    real kmin, kmax;    // min/max k-values covered by interpolating function
    real pmin, pmax;    // corresponding P(k) values
    Spline interp;      // interpolation of data points
};

#endif // INTERPOLATED_PS_H
