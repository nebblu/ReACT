%{
#include "CorrelationFunction.h"
%}

class CorrelationFunction {
public:
    CorrelationFunction(const PowerSpectrum& P, real kmin = 1e-4, real kmax = 1e1);

    real Evaluate(real r) const;
    array EvaluateMany(const array& r) const;

    %extend {
        real __call__(real r) const { return $self->Evaluate(r); }
        array __call__(const array& r) const { return $self->EvaluateMany(r); }
    };

    const PowerSpectrum& GetPowerSpectrum() const { return P; }

protected:
    const PowerSpectrum& P;
    real kmin, kmax;
};
