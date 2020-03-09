#include <cmath>
#include <cstdio>

#include "Cosmology.h"
#include "LinearPS.h"
#include "RPT.h"

int main(int argc, char* argv[]) {
    
    real z0 = 100.;
    real z = 0.;
    real kcut = 1.;
    int Neta = 20;
    int Nk = 50;

    Cosmology C("rpt");
    LinearPS P0(C, z0);
    RPT rpt(C, P0, z0, z, Neta, Nk, kcut);

    int Nout = 100;
    real kminout = 0.001;
    real kmaxout = 0.5;
    real k[Nout], p1[Nout], p2[Nout], p3[Nout];
    for(int i = 0; i < Nout; i++) {
        k[i] = kminout * exp(i*log(kmaxout/kminout)/(Nout-1));
        p1[i] = rpt.P1_11(k[i]);
        p2[i] = rpt.P2_11(k[i]);
        p3[i] = rpt.P3_11(k[i]);
        printf("%e %e %e %e\n", k[i], p1[i], p2[i], p3[i]); fflush(stdout);
    }

//    printf("# k -- P1 -- P2 -- P3\n");
//    for(int i = 0; i < Nout; i++)
//        printf("%e %e %e %e\n", k[i], p1[i], p2[i], p3[i]);

    return 0;
}
