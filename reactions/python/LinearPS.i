%{
#include "LinearPS.h"
%}

class LinearPS : public PowerSpectrum {
public:
    LinearPS(const Cosmology& C, real z = 0);

    const Cosmology& GetCosmology() const;
};
