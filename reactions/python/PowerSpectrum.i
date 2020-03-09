%{
#include "PowerSpectrum.h"
%}

class PowerSpectrum {
public:
    PowerSpectrum();
    virtual ~PowerSpectrum();

    %extend {
        real __call__(real k) const { return $self->Evaluate(k); }
        array __call__(const array& k) const { return $self->EvaluateMany(k); }
    };

    const Cosmology& GetCosmology() const = 0;

    real Sigma(real R) const;
    real VelocityDispersion() const;
    real NonlinearScale() const;
    void Save(const char* filename, real kmin = 1e-3, real kmax = 1, int Nk = 1001, bool log = false);
};

#if 0
class RedshiftPowerSpectrum {
public:
    RedshiftPowerSpectrum();

    real Evaluate(real k, real mu) const = 0;
    real EvaluateMoment(real k, int ell) const = 0;
    %extend {
        real __call__(real k, real mu) { return $self->Evaluate(k, mu); }
    };

};
#endif
