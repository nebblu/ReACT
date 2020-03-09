#ifndef HIGHER_STAT_H
#define HIGHER_STAT_H
#include "Common.h"
#include "MonteCarlo.h"

class HigherStat {
public:
  HigherStat(const Cosmology& C, const PowerSpectrum& P_l, real epsrel = 1e-4);
      ~HigherStat();

/*SKEWNESS */
// a = 1 : tophat
// a = 2 : exponential
// a = 3 : no filter
// b=1 analytic
// b=2 numerical
// R : smoothing scale

  double sig_sqr(double R, int a, int b, double vars[]) const ;

  double skewness(double R, int a, int b, double vars[]) const ;

  double kurtosis(double R,int a, int b, double vars[])const ;

  void sigsqr_init(int a) const;

// private:
//     const Cosmology& C;
//     const PowerSpectrum& P_l;
//     real epsrel;

  protected:
    const Cosmology& C;
    const PowerSpectrum& P_l;
    real epsrel;
    struct kurtosis_integral : public MonteCarloIntegral {
        kurtosis_integral(const HigherStat& hs);
        virtual void Integrand(const double x[], double* f, double* param) const;

        const HigherStat& hs;
    } kintegral;

};

#endif // HIGHER_STAT_H
