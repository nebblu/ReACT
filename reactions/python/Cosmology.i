%{
#include "Cosmology.h"
%}

struct Cosmology {
    Cosmology();
    Cosmology(real h, real n, real Omega_m, real Omega_b);
    Cosmology(real h, real n, real Omega_m, real Omega_b, real sigma8, const array& ki, const array& Ti);
    Cosmology(const char* name);
    ~Cosmology();

    void Initialize(real h, real n, real Omega_m, real Omega_b);
    void Initialize(const char* name);
    void SetTransferFunction(const array& ki, const array& Ti);
    void LoadTransferFunction(const char* tkfile, int kcol = 1, int tcol = 2);
    void NormalizeTransferFunction(real sigma8);
    void CalculateAdditionalParameters();

%immutable;
    /* Input parameters */
    real h;                     // Hubble parameter today
    real Tcmb;                  // CMB temperature today
    real n;                     // scalar spectral index
    real Omega_m;               // matter density parameter today
    real Omega_b;               // baryon density parameter today
    real sigma8;                // power spectrum variance smoothed at 8 Mpc/h scales
    array ki;                   // wavenumbers for transfer function (h/Mpc)
    array Ti;                   // transfer function

    /* Calculated parameters */
    real H0;                    // Hubble parameter today
    real rho_crit;              // critical density today
    real Omega_gamma;
    real Omega_nu;
    real Omega_r;
    real Omega_Lambda;
    real a_eq;                  // scale factor at equality
    real z_eq;                  // redshift of equality
    real k_eq;
    real delta_H;               // normalization of linear power spectrum at z = 0

%mutable;

    real H(real a) const;
    real dHda(real a) const;
    void Print(const char* prefix = "");
};
