#include <cmath>
#include <cstdio>

#include "Common.h"
#include "Quadrature.h"
#include "Timer.h"

DECLARE_SUB(TanSub, tan(u), atan(x), 1/pow2(cos(u)))

double f1(double x) {
//    return 2*pow2(sin(M_PI*x));
    return x*exp(-x);
}

double g1(double x) {
    return exp(1-x);
}

double f2(double x, double y) {
//    return sin(10*(x*x + y*y));
    return x*y*y;
}

int main(int argc, char* argv[]) {
    double eps = 1e-5;

    {
        double a = 0.;
        double b = 1.;

        double I1 = Integrate(f1, a, b, 1e-12);
        printf("I1 = %15.15f\n", I1);

        I1 = Integrate(f1, a, b, eps);
        printf("I1 = %15.15f\n", I1);

        I1 = Integrate<TanSub>(f1, a, b, eps);
        printf("I1 = %15.15f\n", I1);

        I1 = Integrate<1>(f1, &a, &b, eps);
        printf("I1 = %15.15f\n", I1);
    }

    {
        double a = 1.;
        double b = 1e20;

        double J1 = Integrate<ExpSub>(g1, a, b, eps);
        printf("J1 = %15.15f\n", J1);

        J1 = Integrate<InverseSub>(g1, a, b, eps);
        printf("J1 = %15.15f\n", J1);
    }


    {
        double I2;
        Timer T;
        real a[2] = { 0, 0 };
        real b[2] = { 1, 1 };
        for(int i = 0; i < 10000; i++)
            I2 = Integrate<2>(f2, a, b, eps);
        printf("I2 = %15.15f\n", I2);
        printf("t = %d seconds\n", T.WallTimeElapsed());
    }

    return 0;
}
