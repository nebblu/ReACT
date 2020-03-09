#ifndef POWER_SPECTRUM_H
#define POWER_SPECTRUM_H

/* Conventions for power spectra:
 *   - k is measured in h/Mpc
 *   - P(k) is measured in (Mpc/h)^3
 *   - My Fourier conventions are
 *       \tilde f(k) = \int dx e^{-ikx} f(x)
 *       f(x) = \int \frac{dk}{2\pi} e^{ikx} \tilde f(k)
 *     which means P(k) is defined by
 *       (2\pi)^3 \delta_D(k+k') P(k) = < \delta(k) \delta(k') >
 *     This differs by (2\pi)^3 from many PT references. */

#include "Common.h"
#include "array.h"


/**********************************************************************
 * PowerSpectrum
 *
 * Base class for real-space power spectra.
 *********************************************************************/

class PowerSpectrum {
public:
    PowerSpectrum();
    virtual ~PowerSpectrum();

    virtual const Cosmology& GetCosmology() const = 0;

    /* Evaluate power spectrum at a single k */
    virtual real Evaluate(real k) const = 0;
    real operator()(real k) const { return Evaluate(k); }

    /* Evaluate power spectrum at many k values (parallelized for speed) */
    virtual array EvaluateMany(const array& k) const;
    array operator()(const array& k) const { return EvaluateMany(k); }

    /* Calculate the variance with top-hat smoothing radius R Mpc/h */
    virtual real Sigma(real R) const;

    /* Calculate the 1-D velocity dispersion $\sigma_v^2 = \frac{1}{6\pi^2} \int_0^\infty P(k) ~dk$ */
    virtual real VelocityDispersion() const;

    /* Calculate the non-linear scale $k_\text{nl} = 1/\sigma_v$ */
    virtual real NonlinearScale() const;

    /* Write power spectrum to file */
    virtual void Save(const char* filename, real kmin = 1e-3, real kmax = 1, int Nk = 1000, bool log = false);
};


#if 0
/**********************************************************************
 * RedshiftPowerSpectrum
 *
 * Base class for redshift-space power spectra.
 *********************************************************************/

class RedshiftPowerSpectrum {
public:
    RedshiftPowerSpectrum();
    virtual ~RedshiftPowerSpectrum();

    virtual real Evaluate(real k, real mu) const = 0;
    real operator()(real k, real mu) const { return Evaluate(k, mu); }

    virtual real EvaluateMoment(real k, int ell) const = 0;

};
#endif

#endif // POWER_SPECTRUM_H
