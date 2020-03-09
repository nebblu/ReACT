#ifndef NOWIGGLE_PS_H
#define NOWIGGLE_PS_H

#include "Common.h"
#include "Cosmology.h"
#include "GrowthFunction.h"
#include "PowerSpectrum.h"

enum {
    EisensteinHu = 0,           // Eisenstein & Hu 1998 (astro-ph/9709112)
    BBKS,                       // Bardeen, Bond, Kaiser, Szalay 1986
    EfstathiouBondWhite         // Efstathiou, Bond, White 1992 (requires n = 1)
};

/**********************************************************************
 * NoWigglePS
 *
 * Parametric power spectrum based on various fits to data/simulations.
 *********************************************************************/
class NoWigglePS : public PowerSpectrum {
public:
    NoWigglePS(const Cosmology& C, real z = 0, int formula = EisensteinHu);

    const Cosmology& GetCosmology() const { return C; }
    real Evaluate(real k) const;

private:
//    const Cosmology& C;
    Cosmology C;
    real z;
    real Dz;    // normalized growth function, D(z)/D(0)
    int formula;
};

void NoWigglePrintGnuplotFormula(const Cosmology& C, real z = 0, int formula = EisensteinHu);

#endif // NOWIGGLE_PS_H
