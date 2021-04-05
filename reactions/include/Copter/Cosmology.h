/* TODO:
 *  - generalize description of cosmology (maybe make Cosmology abstract,
 *    subclassed by LCDMCosmology, etc.) */

#ifndef COSMOLOGY_H
#define COSMOLOGY_H

#include "Common.h"
#include "array.h"

/******************************************************************************
 * Cosmology
 *
 * Contains all information needed to specify a cosmology, for the purposes of
 * perturbation theory.  The necessary parameters (and WMAP7+ values) are
 *   h, the Hubble parameter today (0.704)
 *   n, the scalar spectral index (0.963)
 *   Omega_m, the matter density parameter (0.272)
 *   Omega_b, the baryon density parameter (0.0456)
 *   sigma8, the linear power spectrum normalization factor (0.809)
 * In addition, in order to compute a realistic linear power spectrum, you must
 * specify a transfer function.  This is typically provided by a Boltzmann
 * code such as CMBFast or CAMB.
 *
 * In almost all cases you will want to instantiate a Cosmology from file.  The
 * file format is largely self-explanatory.  See data/example.ini for a sample.
 ******************************************************************************/
struct Cosmology {
    Cosmology();
    Cosmology(real h, real n, real Omega_m, real Omega_b);
    Cosmology(real h, real n, real Omega_m, real Omega_b, real sigma8, const array& ki, const array& Ti);
    Cosmology(real h, real n, real Omega_m, real Omega_b, real As, real k0, const array& ki, const array& Ti);
    Cosmology(const char* name);

    void Initialize(real h, real n, real Omega_m, real Omega_b);
    void Initialize_as(real h_, real n_, real Omega_m_, real Omega_b_, real As_, real k0_);
    void Initialize(const char* name);
    void SetTransferFunction(const array& ki, const array& Ti);
    void LoadTransferFunction(const char* tkfile, int kcol = 1, int tcol = 2);
    void NormalizeTransferFunction(real sigma8);
    void CalculateAdditionalParameters();

    /* Input parameters */
    real h;                     // Hubble parameter today
    real Tcmb;                  // CMB temperature today
    real n;                     // scalar spectral index
    real Omega_m;               // matter density parameter today
    real Omega_b;               // baryon density parameter today
    real sigma8;                // power spectrum variance smoothed at 8 Mpc/h
    real As;                // Initial amplitude
    real k0;                // Pivot scale
    array ki;                   // wavenumbers for transfer function (h/Mpc)
    array Ti;                   // transfer function

    /* Calculated parameters */
    real H0;                    // Hubble parameter today (H0 = 100*h km/sec/Mpc)
    real rho_crit;              // critical density today
    real Omega_gamma;           // photon density parameter
    real Omega_nu;              // neutrino density parameter
    real Omega_r;               // radiation (photon + neutrino) density parameter
    real Omega_Lambda;          // dark energy density parameter
    real a_eq;                  // scale factor at equality
    real z_eq;                  // redshift of equality
    real k_eq;
    real delta_H;               // normalization of linear power spectrum at z = 0

    /* Hubble parameter at scale factor a */
    real H(real a) const {
        return H0 * sqrt(Omega_m/pow3(a) + Omega_Lambda);
    }

    /* dH/da */
    real dHda(real a) const {
        return -3*pow2(H0)*Omega_m/(2*H(a)*pow4(a));
    }

    /* Print cosmological parameters to stdout */
    void Print(const char* prefix = "") const;
};

#endif // COSMOLOGY_H
