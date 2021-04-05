#include "Common.h"
#include "Cosmology.h"
#include "FlowingWithTime.h"
#include "LinearPS.h"
#include "Timer.h"

int main(int argc, char* argv[]) {
    const char* cosmology = "M000";
    real z = 0.;

    Cosmology C(cosmology);
    LinearPS P_L(C, z);
    vector<InterpolatedP_ab> P;
    for(int i = 0; i < 5; i++)
        P.push_back(InterpolatedP_ab(C, C.ki.size(), &C.ki[0], &C.Ti[0], &C.Ti[0], &C.Ti[0]));
    info("%e\n", P[0](1,1, 0.2));

//    printf("# k -- P1 -- P2 -- P3\n");
//    for(int i = 0; i < Nout; i++)
//        printf("%e %e %e %e\n", k[i], p1[i], p2[i], p3[i]);

    return 0;
}
