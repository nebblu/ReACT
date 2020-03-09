#ifndef KAISER_H
#define KAISER_H

#include "Common.h"
#include "array.h"

/* Formulas for Kaiser's linear theory redshift-space power spectrum, and
 * Hamilton's corresponding formulas for the redshift-space correlation
 * function. */
class Kaiser {
public:
    Kaiser(const Cosmology& C, const PowerSpectrum& P_L, real z);
    Kaiser(const PowerSpectrum& P_L, real f);

    /* Redshift-space power spectrum */
    real P_s(real k, real mu) const;

    /* Multipoles of redshift-space power spectrum */
    real PMultipole(int ell, real k) const;
    array PMultipole(int ell, const array& k) const;

    /* Multipoles of redshift-space correlation function */
    real XiMultipole(int ell, real r, int Nk = 32768, real kmin = 0, real kmax = 100) const;
    array XiMultipole(int ell, const array& r, int Nk = 32768, real kmin = 0, real kmax = 100) const;

//    Spline ComputeP_ell(int ell, real kmin = 1e-4, real kmax = 1e1, int Nk = 1024);
//    Spline ComputeXi_ell(int ell, real kmin = 1e-4, real kmax = 1e1, int Nk = 1024);

//    static void ComputeXi_ell(int ell, real f, const array& k, const array& pk, array& r, array& xi);


protected:
    const PowerSpectrum& P;
    real f;
};

#endif // KAISER_H
