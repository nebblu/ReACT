%{
#include "Spline.h"
%}

class Spline {
public:
    real Evaluate(real x);
    real EvaluateDerivative(real x);

    %extend {
        real __call__(real x) {
            return $self->Evaluate(x);
        }

        array __call__(const array& X) {
            int N = (int)X.size();
            array Y(N);
            for(int i = 0; i < N; i++)
                Y[i] = $self->Evaluate(X[i]);
            return Y;
        }
    };
};

Spline LinearSpline(const array& X, const array& Y);
Spline ShiftedLinearSpline(const array& X, const array& Y, real tau = 0.2);
Spline CubicSpline(const array& X, const array& Y);
