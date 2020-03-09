#ifndef SPLINE_H
#define SPLINE_H

/* TODO:
 * - Allow user-defined behavior for out-of-range requests (i.e. x < xmin or
 *   x > xmax).
 * - Also allow configurable index lookup, in cases where the domain points
 *   are known to be spaced linearly, logarithmically, etc.
 */

#include <vector>
using std::vector;

#include "Common.h"

class Spline;
class SplineImpl;

/***** Spline factory functions *****/

/* Linear interpolation */
Spline LinearSpline(const vector<real>& X, const vector<real>& Y);
Spline LinearSpline(int N, const real* X, const real* Y);

/* Shifted linear interpolation (see Blu, Thevenaz, and Unser, 2004) */
Spline ShiftedLinearSpline(const vector<real>& X, const vector<real>& Y, real tau = 0.2);
Spline ShiftedLinearSpline(int N, const real* X, const real* Y, real tau = 0.2);

/* Natural cubic spline */
Spline CubicSpline(const vector<real>& X, const vector<real>& Y);
Spline CubicSpline(int N, const real* X, const real* Y);


/***************************************************************
 * Spline
 *
 * Generic spline wrapper class.
 ***************************************************************/
class Spline {
public:
    /* Default to cubic spline */
//    Spline(const vector<real>& X, const vector<real>& Y);
//    Spline(int N, const real* X, const real* Y);

    Spline(SplineImpl* impl = NULL);
    Spline(const Spline& F);
    Spline& operator=(const Spline& F);
    ~Spline();

    real Evaluate(real x) const;
    real EvaluateDerivative(real x) const;

    real operator()(real x) const { return Evaluate(x); }

    /* Find a local maximum or minimum of the interpolated function */
//    real FindMaximum(real xguess, real* ymax = 0);
//    real FindMinimum(real xguess, real& ymin = 0);

protected:
    SplineImpl* impl;   // internal spline implementation
};


/**************************************************
 * SplineImpl
 *
 * Base class for internal spline implementations.
 **************************************************/
struct SplineImpl {
    int refcount;

    SplineImpl() { refcount = 0; }
    virtual ~SplineImpl() {}

    virtual real y(real x) const = 0;
    virtual real dydx(real x) const = 0;
    virtual SplineImpl* clone() const = 0;

    real xmin, xmax;    // domain
};

#endif // SPLINE_H
