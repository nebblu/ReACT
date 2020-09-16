#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <cstdio>
#include <unistd.h>

#include "Cosmology.h"
#include "Datafile.h"
#include "SpecialFunctions.h"
#include "cfg.h"
#include "pstring.h"


Cosmology::Cosmology() {
    /* No initialization */
}

Cosmology::Cosmology(real h, real n, real Omega_m, real Omega_b) {
    Initialize(h, n, Omega_m, Omega_b);
}

Cosmology::Cosmology(real h, real n, real Omega_m, real Omega_b, real sigma8, const array& ki, const array& Ti) {
    Initialize(h, n, Omega_m, Omega_b);
    SetTransferFunction(ki, Ti);
    NormalizeTransferFunction(sigma8);
}

Cosmology::Cosmology(real h, real n, real Omega_m, real Omega_b, real As, real k0, const array& ki, const array& Ti) {
    Initialize_as(h, n, Omega_m, Omega_b, As, k0);
    SetTransferFunction(ki, Ti);
}

Cosmology::Cosmology(const char* cosmo) {
    sigma8 = 0;
    delta_H = 1;
    Initialize(cosmo);
}

void Cosmology::Initialize_as(real h_, real n_, real Omega_m_, real Omega_b_, real As_, real k0_) {
    h = h_;
    n = n_;
    Omega_m = Omega_m_;
    Omega_b = Omega_b_;
    As = As_;
    k0 = k0_;

    sigma8 = 0;
    delta_H = 1;

    CalculateAdditionalParameters();
}


void Cosmology::Initialize(real h_, real n_, real Omega_m_, real Omega_b_) {
    h = h_;
    n = n_;
    Omega_m = Omega_m_;
    Omega_b = Omega_b_;

    sigma8 = 0;
    delta_H = 1;

    CalculateAdditionalParameters();
}

void Cosmology::Initialize(const char* cosmo) {
    pstring path(cosmo);

    /* First check the current directory for a file named cosmo */
    bool found = (access(path.c_str(), R_OK) == 0);

    /* Next try adding a file extension to cosmo */
    if(!found) {
        path += ".ini";
        found = (access(path.c_str(), R_OK) == 0);
    }

    /* Finally search in the default data directory */
    if(!found) {
        path = DATADIR "/" + path;
        found = (access(path.c_str(), R_OK) == 0);
    }

    if(!found) {
        warning("Cosmology: could not find cosmology '%s'\n", cosmo);
        warning("  tried ./%s, ./%s.ini, and %s/%s.ini\n", cosmo, cosmo, DATADIR, cosmo);
        return;
    }

    verbose("Reading cosmology from %s\n", path.c_str());

    /* Read in parameters from cosmology file */
    Config cfg = cfg_new_from_file(path.c_str());
    if(!cfg_has_keys(cfg, "h,n,Omega_m,Omega_b,sigma8,tkfile", ","))
        error("Cosmology: not all parameters defined in %s\n", path.c_str());
    h = cfg_get_double(cfg, "h");
    n = cfg_get_double(cfg, "n");
    Omega_m = cfg_get_double(cfg, "Omega_m");
    Omega_b = cfg_get_double(cfg, "Omega_b");
    sigma8 = cfg_get_double(cfg, "sigma8");

    /* Load transfer function */
    const char* tkfile = cfg_get(cfg, "tkfile");
    int kcol = cfg_has_key(cfg, "tkfile:kcol") ? cfg_get_int(cfg, "tkfile:kcol") : 1;
    int tcol = cfg_has_key(cfg, "tkfile:tcol") ? cfg_get_int(cfg, "tkfile:tcol") : 2;
    LoadTransferFunction(tkfile, kcol, tcol);

    CalculateAdditionalParameters();
    if(!Ti.empty() && sigma8 != 0)
        NormalizeTransferFunction(sigma8);

    cfg_destroy(cfg);
}

void Cosmology::SetTransferFunction(const array& ki_, const array& Ti_) {
    ki = ki_;
    Ti = Ti_;
}

void Cosmology::LoadTransferFunction(const char* tkfile, int kcol, int tcol) {
    /* First check the current directory for a file named 'tkfile' */
    pstring tkpath(tkfile);
    FILE* fp = fopen(tkpath.c_str(), "r");

    if(!fp) {
        /* Next search in the default data directory */
        tkpath = DATADIR "/" + tkpath;
        fp = fopen(tkpath.c_str(), "r");
    }

    if(!fp) {
        fprintf(stderr, "Cosmology: could not find transfer function '%s'\n", tkfile);
        fprintf(stderr, "  tried ./%s and %s/%s\n", tkfile, DATADIR, tkfile);
        return;
    }

    verbose("Reading transfer function from %s\n", tkpath.c_str());

    Datafile data(fp);
    ki = data.GetColumn(kcol);
    Ti = data.GetColumn(tcol);
}

void Cosmology::NormalizeTransferFunction(real sigma8) {
    /* Calculate sigma8 directly from the transfer function (with delta_H = 1) */
    const real r = 8;
    int N = ki.size() - 1;      // number of subintervals
    array fi(N+1);
    for(int i = 0; i <= N; i++) {
        real k = ki[i];
        real T = Ti[i];
        real kr = k*r;
        real j1_kr = (kr < 1e-4) ? 1/3. - kr*kr/30. : SphericalBesselJ1(kr)/kr;
        fi[i] = pow(2997.925*k, 3+n) * pow2(T) * pow2(3*j1_kr);
    }
    real S = 0;
    for(int i = 0; i <= N-1; i++)
        S += log(ki[i+1]/ki[i]) * 0.5*(fi[i] + fi[i+1]);        // trapezoid rule in log(k)

    /* Set delta_H to get the desired sigma8 */
    delta_H = sigma8/sqrt(S);
    verbose("Normalizing linear P(k) to sigma8 = %g: delta_H = %e\n", sigma8, delta_H);
}

void Cosmology::CalculateAdditionalParameters() {
    /* These are just textbook definitions for a few miscellaneous parameters.
     * They will probably never be used. */

    using namespace Constants;
    Tcmb = 2.725;
    H0 = 100*h * km/second/Mpc;
    rho_crit = (3*H0*H0)/(8*M_PI*G);
    Omega_gamma = (M_PI*M_PI/15) * pow4(k*Tcmb)/(pow3(hbar)*pow5(c)) / rho_crit;
    Omega_nu = 3.*(7./8.) * pow(4./11., 4./3.) * Omega_gamma;
    Omega_r = Omega_gamma + Omega_nu;
    Omega_Lambda = 1 - Omega_m;
    a_eq = Omega_r/Omega_m;
    z_eq = 1/a_eq - 1;
    k_eq = (H0/c) * sqrt(2*Omega_m/a_eq) / (h/Mpc);
}

void Cosmology::Print(const char* prefix) const {
    printf("%sh = %g\n", prefix, h);
    printf("%sn = %g\n", prefix, n);
    printf("%sOmega_m = %g\n", prefix, Omega_m);
    printf("%sOmega_b = %g\n", prefix, Omega_b);
    printf("%ssigma8 = %g\n", prefix, sigma8);
    printf("%sTcmb = %g\n", prefix, Tcmb);
    printf("%sH0 = %g\n", prefix, H0);
    printf("%srho_crit = %g\n", prefix, rho_crit);
    printf("%sOmega_gamma = %g\n", prefix, Omega_gamma);
    printf("%sOmega_nu = %g\n", prefix, Omega_nu);
    printf("%sOmega_r = %g\n", prefix, Omega_r);
    printf("%sOmega_Lambda = %g\n", prefix, Omega_Lambda);
    printf("%sa_eq = %g\n", prefix, a_eq);
    printf("%sz_eq = %g\n", prefix, z_eq);
    printf("%sk_eq = %g\n", prefix, k_eq);
    printf("%sdelta_H = %g\n", prefix, delta_H);
}
