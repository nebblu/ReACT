%{
#include "GrowthFunction.h"
%}

class GrowthFunction {
public:
    GrowthFunction();
    GrowthFunction(const Cosmology& C);
    GrowthFunction(const Cosmology& C, real a_i, real Dp_i);

    real Evaluate(real z) const;
    %extend {
        real __call__(real z) const { return $self->Evaluate(z); }
    }

    /* f = dlogD/dloga */
    real f(real z) const;
};
