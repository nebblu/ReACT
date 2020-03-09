%{
#include "SimpleRG.h"
%}

class SimpleRG {
public:
    SimpleRG(const Cosmology& C, const PowerSpectrum& P_L);
    
    InterpolatedPS CalculateP(int Nk, real kmin, real kmax, real epsrel = 1e-3);
};
