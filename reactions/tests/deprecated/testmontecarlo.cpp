#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "Common.h"
#include "MonteCarlo.h"


class TestIntegral : public MonteCarloIntegral {
public:
    TestIntegral()
        : MonteCarloIntegral(3, 1)
    {
        maxeval = 100000;
    }

    void Integrand(const double x[], double* f, double*) const {
        f[0] = sin(x[0]) * cos(x[1]) * exp(x[2]);
    }
};

class TestIntegral2 : public MonteCarloIntegral {
public:
    TestIntegral2()
        : MonteCarloIntegral(5, 1)
    {
        double kmin = 1e-4;
        double kmax = 1e2;
        Range domain[5] = {
            { log(kmin), log(kmax) },
            { log(kmin), log(kmax) },
            { -1, 1 },
            { -1, 1 },
            { 0, 2*M_PI }
        };
        SetDomain(domain);

//        epsrel = 1e-4;
//        maxeval = 100000000;
    }

    void Integrand(const double x[], double* f, double* param) const {
        double k = *param;
        double p = exp(x[0]);
        double q = exp(x[1]);
        double mu_p = x[2];
        double mu_q = x[3];
        double phi_q = x[4];
        double kp2 = k*k + p*p - 2*k*p*mu_p;
        double pq2 = p*p + q*q - 2*p*q*(sqrt(1-mu_p*mu_p)*sqrt(1-mu_q*mu_q)*cos(phi_q) + mu_p*mu_q);
        f[0] = 1/pow3(M_PI) * 2*M_PI*pow3(p)*pow3(q) * exp(-kp2 - pq2);
    }
};


int main(int argc, char* argv[]) {
    double integral, error;
    int neval;

    TestIntegral I;
    I.Integrate(&integral, 0, &error, &neval, 0, 0.001);
    printf("I = %.15f +/- %.15f\n", integral, error);
    printf("[expecting (1-cos(1)) * sin(1) * (exp(1)-1) = 0.664669679978137714]\n");
    printf("neval = %d\n", neval);

    double k = 0.2;
    TestIntegral2 I2;
    I2.Integrate(&integral, &k, &error, &neval, 0, 0.01);
    printf("I2 = %.15f +/- %.15f\n", integral, error);
    printf("[expecting 1]\n");
    printf("neval = %d\n", neval);

    return 0;
}
