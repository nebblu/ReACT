#ifndef ODE_H
#define ODE_H

#include <vector>
using std::vector;

#include "Common.h"

/******************************************************************************
 * RK4: 4th order Runge-Kutta method.
 *   Integrates y(t) over the interval [t0,t1] using N fixed time steps of
 *   the standard Runge-Kutta 4 method.
 * - In the first form (for a single dependent variable), dydt must have the signature
 *     real dydt(real t, real y);
 *   the returned array gives the solution y(t) evaluated at the fixed times
 *     t[i] = t0 + i*(t1 - t0)/N.
 * - In the second form (for m dependent variables0, dydt must have the signature
 *     vector<real> dydt(real t, const vector<real>& y);
 *   the 2-dimensional returned array y[j][i] gives the solutions y_j(t) (1 <= j <= m)
 *   evaluated at the fixed times
 *     t[i] = t0 + i*(t1 - t0)/N.
 */

template<typename Function>
vector<real> RungeKutta4(real t0, real t1, real y0, Function dydx, int N = 1000);

template<typename Function>
vector<vector<real> > RungeKutta4(real t0, real t1, const vector<real>& y0, Function dydx, int N = 1000);

/******************************************************************************
 * RKDP: Adaptive Dormand-Prince method.
 *   Integrates y over the interval [t0,t1], starting with initial value y0.
 * - In the first form (for a single dependent variable), dydt must have the signature
 *      real dydt(real t, real y);
 *   the output arrays t and y give the solution y(t) over the entire interval [t0,t1].
 * - In the second form (for multiple dependent variables), the signature must be
 *      vector<real> dydt(real t, const vector<real>& y);
 *   the output array y1 gives the solution only at the final time t1.
 * - Both forms return the number of time steps taken.
 * - epsrel is the local error tolerance (global error may be larger, but typically
 *   of the same order of magnitude).
 * - h0 is the initial step size; if set to 0, the initial step size defaults to t1-t0. */

template<typename Function>
int RKDP(real t0, real t1, real y0, Function dydt, vector<real>& t, vector<real>& y, real epsrel = 1e-3, real h0 = 0);

template<typename Function>
int RKDP(real t0, real t1, const vector<real>& y0, Function dydt, vector<real>& y1, real epsrel = 1e-3, real h0 = 0);

#include "ODE.inl"

#endif // ODE_H
