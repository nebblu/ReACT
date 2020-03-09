/* Test of Matsubara's Lagrangian resummation scheme. */

#include <cmath>

#include "array.h"
#include "Cosmology.h"
#include "LagrangianResummation.h"
#include "LinearPS.h"
#include "Timer.h"

int main(int argc, char* argv[]) {
    const char* cosmology = "M000";
    real z = 0.;

    Cosmology C(cosmology);
    LinearPS P_L(C, z);
    LagrangianResummation lr(C, z);

    int Nk = 100;
    real kmin = 0.001, kmax = 0.5;
    array k = array::logspace(kmin, kmax, Nk);
    array plin = array::zeros(Nk);
    array pmat = array::zeros(Nk);

    Timer T;
    for(int i = 0; i < Nk; i++) {
        plin[i] = P_L(k[i]);
        pmat[i] = lr.P(k[i]);
        info("%e %e %e\n", k[i], plin[i], pmat[i]);
    }
    info("Elapsed time: %d seconds\n", T.WallTimeElapsed());

    return 0;
}
