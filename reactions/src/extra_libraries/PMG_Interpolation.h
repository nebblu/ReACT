#ifndef PMG_H
#define PMG_H
#include "Common.h"

class PMGI {
public:
  PMGI(const Cosmology& C, const PowerSpectrum& P_l, real epsrel = 1e-4);

  void spline_init(double scalef, double omega0, double mgvals[], const int mgsize, int a)const;
  double PMG_TNS(double k, int a, double bl, double sigma_v, double mg1)const;
  void PMG_FREE()const;


private:
    const Cosmology& C;
    const PowerSpectrum& P_l;
    real epsrel;

};

#endif // PMG
