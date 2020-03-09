%{
#include "InterpolatedPS.h"
%}

class InterpolatedPS : public PowerSpectrum {
public:
    InterpolatedPS(const Cosmology& C, const char* filename, int kcol = 1, int pcol = 2);
    InterpolatedPS(const Cosmology& C, const array& k, const array& p);
    InterpolatedPS(const Cosmology& C, const array& k, const array& p, bool logscale);

    const Cosmology& GetCosmology() const;
};
