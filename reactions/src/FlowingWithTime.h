#ifndef FLOWINGWITHTIME_H
#define FLOWINGWITHTIME_H

#include "Common.h"
#include "InterpolatedPS.h"
#include "PowerSpectrum.h"
#include "array.h"

/* Pietroni's time-RG theory
 *
 * Based on M. Pietroni, "Flowing with time: a new approach to non-linear
 * cosmological perturbations", 2008.  See Appendix B of that paper for an
 * explanation of the code below. */


/* Interpolated power spectrum $P_{ab}(k)$ */
struct InterpolatedP_ab {
    int N;
    real ka, kb, p11a, p11b, p12a, p12b, p22a, p22b;
    InterpolatedPS P_11, P_12, P_22;

    InterpolatedP_ab(const Cosmology& C, int N, const real* k, const real* p11, const real* p12, const real* p22);
    real operator()(int a, int b, real k) const;
};


class FlowingWithTime {
public:
    FlowingWithTime(const Cosmology& C, real z_i, const PowerSpectrum& P_i, const array& k);
    ~FlowingWithTime();

    /* Return a list of power spectra, evaluated at z[i] */
    vector<InterpolatedP_ab> CalculateP_ab(int Nz, const real z[]) const;

    real A(const InterpolatedP_ab& P, int a, int c, int d, int b, int e, int f, real k) const;

    vector<real> dXdeta(real eta, vector<real>& X) const;

protected:
    const Cosmology& C;
    real a_i;
    const PowerSpectrum& P_i;
    int Nk;
    array k;
};

#endif // FLOWINGWITHTIME_H
