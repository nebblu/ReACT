#ifndef GROWTH_FUNCTION_H
#define GROWTH_FUNCTION_H

#include <vector>

#include "Common.h"
#include "Cosmology.h"
#include "Spline.h"

/* Linear growth function, normalized to D=1 at z=0.  (Set Dp_i = -1.5 to
 * get the decaying mode solution. */
class GrowthFunction {
public:
    GrowthFunction();
    GrowthFunction(const Cosmology& C);
    GrowthFunction(const Cosmology& C, real a_i, real Dp_i);

    real Evaluate(real z) const;
    real operator()(real z) const { return Evaluate(z); }

    /* Logarithmic growth rate: f = dlogD/dloga */
    real f(real z) const;

private:
    Cosmology C;

    Spline D;           // D(a)
    Spline dDda;        // dD/da

    /* X = (D, dD/da) */
    std::vector<real> dXda(real a, const std::vector<real>& X);
};

#endif // GROWTH_FUNCTION_H
