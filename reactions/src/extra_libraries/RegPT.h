#ifndef REG_PT_H
#define REG_PT_H
#include "Common.h"

class RegPT {
public:
  RegPT(const Cosmology& C, const PowerSpectrum& P_l, real epsrel = 1e-4);

  // used to fix log growth
  void rempreg(real f) const;
  // used to fix linear growth
  void Greg(real dfid) const;

  // Initialization of sigma_d damping term
  // a = 0 :analytic
  // a = 1 numeric
  void sigmad_init() const;

// All functions similar to those found in SPT.h (with reordering of arguments perhaps - laziness!)

  real  PLOOPr(int a, real k) const ;

  real  PLOOPnr(real kmin, real kmax, int a,  real k) const ;

  real PTNSMnDGPr(real k, real bl, real sigma_v, int a) const ;

  real PTNSMmgr(real k, double kmin, double kmax, int a,real bl, real sigma_v ) const;

  real ABr(real k, real bl, real U, int a) const;

  real ABnr(double kmin, double kmax, real bl, real k, real U, int e) const;

// used for interpolation of scale dep models (P0P2)
  void ABr_selec(double myarray[], real k) const ;

private:
    const Cosmology& C;
    const PowerSpectrum& P_l;
    real epsrel;

};

#endif // REG_PT_H
