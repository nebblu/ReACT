#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <vector>

#include "GrowthFunction.h"
#include "Quadrature.h"
#include "ODE.h"
#include "array.h"

#include <functional> 

using std::cref;
using std::bind;
using std::vector;



GrowthFunction::GrowthFunction() {
}

GrowthFunction::GrowthFunction(const Cosmology& C_) {
    C = C_;
    const int N = 10000;
    const real a_i = 1.0/N;
    const real a_f = 1.0 + 1.0/N;;

    /* Solve for D(a) between a_i and a_f, assuming D(a) ~ a for small a */
    vector<real> A(N+1);
    for(int i = 0; i <= N; i++)
        A[i] = a_i + i*(a_f - a_i)/N;
    vector<real> x0(2);
    x0[0] = a_i;
    x0[1] = 1;
    vector<vector<real> > X = RungeKutta4(a_i, a_f, x0, bind(&GrowthFunction::dXda, this, std::placeholders::_1, std::placeholders::_2), N);

    /* Interpolate to a = 0 */
    A.insert(A.begin(), 0.);
    X[0].insert(X[0].begin(), 0.);
    X[1].insert(X[1].begin(), 0.);

    /* Normalize to D=1 at a=1 */
    D = CubicSpline(A, X[0]);
    real D0 = D(1);
    for(int i = 0; i <= N+1; i++) {
        X[0][i] /= D0;
        X[1][i] /= D0;
    }

    D = CubicSpline(A, X[0]);
    dDda = CubicSpline(A, X[1]);
}

GrowthFunction::GrowthFunction(const Cosmology& C_, real a_i, real dDda_i) {
    C = C_;
    const int N = 10000;
    const real a_f = 1.01;

    vector<real> A(N+1);
    for(int i = 0; i <= N; i++)
        A[i] = a_i + i*(a_f - a_i)/N;
    vector<real> x0(2);
    x0[0] = a_i;        // D_i = a_i
    x0[1] = dDda_i;
    vector<vector<real> > X = RungeKutta4(a_i, a_f, x0, bind(&GrowthFunction::dXda, this, std::placeholders::_1, std::placeholders::_2), N);

    /* Normalize to D=1 at a=1 */
    D = CubicSpline(A, X[0]);
    real D0 = D(1);
    for(int i = 0; i <= N; i++) {
        X[0][i] /= D0;
        X[1][i] /= D0;
    }

    D = CubicSpline(A, X[0]);
    dDda = CubicSpline(A, X[1]);
}

real GrowthFunction::Evaluate(real z) const {
    real a = 1/(1+z);
    return D(a);
}

real GrowthFunction::f(real z) const {
    real a = 1/(1+z);
    return a/D(a) * dDda(a);
}

vector<real> GrowthFunction::dXda(real a, const vector<real>& X) {
    vector<real> dX(2);
    dX[0] = X[1];
    dX[1] = -(3/a + C.dHda(a)/C.H(a)) * X[1] + 1.5*C.Omega_m*pow2(C.H0)/(pow5(a)*pow2(C.H(a))) * X[0];
    return dX;
}
