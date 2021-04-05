#include "Common.h"
#include "ODE.h"
#include "array.h"


vector<real> dydt(real t, const vector<real>& y) {
    vector<real> a(2);
    a[0] = -y[1];
    a[1] = y[0];
    return a;
}

int main(int argc, char* argv[]) {
    real ti = 0;
    real tf = 2*M_PI;
    vector<real> y(2);
    y[0] = 1;
    y[1] = 0;
    real epsrel = 1e-10;

    int N = RKDP(ti, tf, y, dydt, y, epsrel);
    info("cos(2*pi) = %e\n", y[0]);
    info("sin(2*pi) = %e\n", y[1]);

    return 0;
}
