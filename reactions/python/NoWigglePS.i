%{
#include "NoWigglePS.h"
%}

enum {
    EisensteinHu,
    BBKS,
    EfstathiouBondWhite
};

class NoWigglePS : public PowerSpectrum {
public:
    NoWigglePS(const Cosmology& C, real z = 0, int formula = EisensteinHu);

    const Cosmology& GetCosmology() const;
};

void NoWigglePrintGnuplotFormula(const Cosmology& C, real z = 0, int formula = EisensteinHu);
