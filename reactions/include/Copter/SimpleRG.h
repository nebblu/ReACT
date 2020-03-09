#ifndef SIMPLE_RG_H
#define SIMPLE_RG_H

#include "Common.h"
#include "Cosmology.h"
#include "InterpolatedPS.h"
#include "array.h"


/**********************************************************************
 * SimpleRG
 *
 * McDonald's "simple renormalization group" power spectrum.
 *********************************************************************/
class SimpleRG {
public:
    SimpleRG(const Cosmology& C, const PowerSpectrum& P_L);
    
    InterpolatedPS CalculateP(int Nk, real kmin, real kmax, real epsrel = 1e-3);

private:
    const Cosmology& C;
    const PowerSpectrum& P_L;
    int Nk;
    vector<real> k;

    vector<real> dQdA(real A, const vector<real>& Q);
};

#endif // SIMPLE_RG_H
