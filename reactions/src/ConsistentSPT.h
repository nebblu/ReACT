#ifndef CONSISTENT_SPT_H
#define CONSISTENT_SPT_H

#include "Common.h"
#include "Cosmology.h"
#include "PowerSpectrum.h"
#include "Spline.h"
#include "array.h"

/* TODO:
 * - use eta = log D for the time variable instead of sigma */


/*******************************************************************************
 * ConsistentSPT
 *
 * Power spectrum from 1-loop standard perturbation theory with consistent time
 * dependence of higher order terms (i.e. without making the Einstein-de Sitter
 * approximation).
 ******************************************************************************/
class ConsistentSPT {
public:
    ConsistentSPT(const Cosmology& C, const PowerSpectrum& P0, real z0, real zdefault = 0, real epsrel = 1e-5);

    /* 1-loop power spectrum */
    real P(real k, int a = 1, int b = 1) const;
    real P_dd(real k, real z = -1) const;
    real P_dt(real k, real z = -1) const;
    real P_tt(real k, real z = -1) const;

    /* 1-loop contributions to PS */
    real P22(real k, int a = 1, int b = 1) const;
    real P13(real k, int a = 1, int b = 1) const;

    real P22_dd(real k, real z = -1) const;
    real P22_dt(real k, real z = -1) const;
    real P22_tt(real k, real z = -1) const;
    real P13_dd(real k, real z = -1) const;
    real P13_dt(real k, real z = -1) const;
    real P13_tt(real k, real z = -1) const;

    void ComputeAll(real k, real& p11, real& p22dd, real& p13dd, real& p22tt, real& p13tt, real& p22dt, real& p13dt, real z) const;

    /* Mode-coupling integrals */
    real P_AA(real k) const;
    real P_AB(real k) const;
    real P_BB(real k) const;
    real P_A0A(real k) const;
    real P_AA0(real k) const;
    real P_A0B(real k) const;
    real P_AB0(real k) const;
    real P_B0A(real k) const;
    real P_B0B(real k) const;

private:
    static const int Nsteps = 10000;

    const Cosmology& C;
    const PowerSpectrum& P0;
    real a0, r, c0;
    real zdefault;
    real A;
    real epsrel, epsabs;

    /* Time-dependent functions */
    Spline U, VA, VB, WA0A, WAA0, WA0B, WAB0, WB0A, WB0B;
    Spline X, YA, YB, ZA0A, ZAA0, ZA0B, ZAB0, ZB0A, ZB0B;

    real sigma(real a) const;
    real epsilon(real s) const;
    vector<real> dfds(real s, const vector<real>& f);
};

#endif // CONSISTENT_SPT_H
