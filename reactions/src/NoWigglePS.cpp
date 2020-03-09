#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <cstdio>

#include "GrowthFunction.h"
#include "NoWigglePS.h"


void NoWigglePrintGnuplotFormula(const Cosmology& C, real z, int formula) {
    GrowthFunction D(C);

    if(formula == EisensteinHu) {
        real h = C.h;
        real Omega_0 = C.Omega_m;
        real Omega_b = C.Omega_b;
        real Theta = C.Tcmb/2.7;
        real alpha_Gamma = 1 - 0.328*log(431*Omega_0*h*h) * Omega_b/Omega_0 + 0.38*log(22.3*Omega_0*h*h) * pow2(Omega_b/Omega_0);
        real s = 44.5 * log(9.83/(Omega_0*h*h)) / sqrt(1 + 10*pow(Omega_b*h*h, 0.75)) * h;

        /* Print reduced formula */
        printf("D = %g\n", D(z)/D(0));
        printf("Gamma_eff(k) = %g*(%g + %g/(1 + %g*k**4))\n", C.Omega_m*C.h, alpha_Gamma, 1-alpha_Gamma, pow4(0.43*s));
        printf("q(k) = %g*k/Gamma_eff(k)\n", pow2(Theta));
        printf("P_nw(k) = %g * D**2 * k**(%g) * (log(%g + 1.8*q(k))/(log(%g + 1.8*q(k)) + (14.2 + 731/(1 + 62.5*q(k))) * q(k)**2))**2\n", 2*M_PI*M_PI*pow2(C.delta_H)*pow(2997.92, 3+C.n), C.n, 2*M_E, 2*M_E);
    }
    else if(formula == BBKS) {
        printf("D = %g\n", D(z)/D(0));
    }
    else if(formula == EfstathiouBondWhite) {
        printf("D = %g\n", D(z)/D(0));
    }
}


NoWigglePS::NoWigglePS(const Cosmology& C_, real z_, int formula_)
    : C(C_)
{
    z = z_;
    formula = formula_;
    GrowthFunction D(C);
    Dz = D(z)/D(0);
}

real NoWigglePS::Evaluate(real k) const {
    if(k <= 0)
        return 0;

    real h = C.h;
    real n = C.n;
    real Tcmb = C.Tcmb;
    real Omega_m = C.Omega_m;
    real Omega_b = C.Omega_b;
    real delta_H = C.delta_H;

    if(formula == EisensteinHu) {
        real Theta = Tcmb/2.7;
        real alpha_Gamma = 1 - 0.328*log(431*Omega_m*h*h) * Omega_b/Omega_m + 0.38*log(22.3*Omega_m*h*h) * pow2(Omega_b/Omega_m);
        real s = 44.5 * log(9.83/(Omega_m*h*h)) / sqrt(1 + 10*pow(Omega_b*h*h, 0.75)) * h;
        real Gamma_eff = Omega_m*h * (alpha_Gamma + (1-alpha_Gamma)/(1 + pow4(0.43*k*s)));
        real q = k * pow2(Theta) / Gamma_eff;
        real L0 = log(2*M_E + 1.8*q);
        real C0 = 14.2 + 731/(1 + 62.5*q);
        real T0 = L0/(L0 + C0*pow2(q));
        return 2*M_PI*M_PI/pow3(k) * pow2(delta_H) * pow(2997.92*k, 3+n) * pow2(T0) * pow2(Dz);
    }
    else if(formula == BBKS) {
        real q = k / (h*(Omega_m - Omega_b));
        real T = log(1 + 2.34*q)/(2.34*q) * pow(1 + 3.89*q + pow2(16.1*q) + pow3(5.46*q) + pow4(6.71*q), -0.25);
        return 2*M_PI*M_PI/pow3(k) * pow2(delta_H) * pow(2997.92*k, 3+n) * pow2(T) * pow2(Dz);
    }
    else if(formula == EfstathiouBondWhite) {
        real B = pow2(delta_H) * pow4(2997.92);
        real Gamma = Omega_m*h;
        real a = 6.4/Gamma;
        real b = 3.0/Gamma;
        real c = 1.7/Gamma;
        real nu = 1.13;
        return B*k / pow(1 + pow(a*k + pow(b*k, 1.5) + pow2(c*k), nu), 2/nu) * pow2(Dz);
    }
    else {
        return 0.;
    }
}
