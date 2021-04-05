#include <cmath>

#include "Common.h"
#include "Cosmology.h"
#include "LinearPS.h"
#include "SPT.h"
#include "Timer.h"

int main(int argc, char* argv[]) {
    const char* cosmology = "M000";
    real z = 0.;

    Cosmology C(cosmology);
    LinearPS P_L(C, z);
    SPT spt(C, P_L);

    int Nout = 100;
    real kminout = 0.001;
    real kmaxout = 0.5;
    real k[Nout], p_11[Nout], p_12[Nout], p_22[Nout];
    Timer T;
    for(int i = 0; i < Nout; i++) {
        k[i] = kminout * exp(i*log(kmaxout/kminout)/(Nout-1));
        p_11[i] = spt.P(k[i], 1, 1);
        p_12[i] = spt.P(k[i], 1, 2);
        p_22[i] = spt.P(k[i], 2, 2);
        info("%e %e %e %e\n", k[i], p_11[i], p_12[i], p_22[i]);
    }
    info("Elapsed time: %d seconds\n", T.WallTimeElapsed());

//    printf("# k -- P1 -- P2 -- P3\n");
//    for(int i = 0; i < Nout; i++)
//        printf("%e %e %e %e\n", k[i], p1[i], p2[i], p3[i]);

    return 0;
}
