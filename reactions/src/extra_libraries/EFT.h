#ifndef EFT_H
#define EFT_H
#include "Common.h"

class EFT {
public:
  EFT(const Cosmology& C, const PowerSpectrum& P_l, real epsrel=1e-3);

   void NWinit(double scalef, double omega0, double mg1, double mg2, double mg3) const;
   void sigmav_init(double scalef, double omega0, double mg1, double mg2, double mg3) const;


   void eft_init(double s8, double scalef, double omega0, double mg1, double mg2, double mg3)const;

// Resummed 1-loop power spectrum
   double Presum(double k, int a) const;

// Individual loop terms with no wiggle spectrum instead of linear
   double Pres22_dd(double k) const;
   double Pres13_dd(double k) const;
   double Pres22_dt(double k) const;
   double Pres13_dt(double k) const;
   double Pres22_tt(double k) const;
   double Pres13_tt(double k) const;

//Normalise linear and logarithmic growth
  void Geft(real dfid) const;
  void rempeft(real ffid) const;



// Total 1-loop power spectra
   double PNWloop(double k, int a) const;

// EFT redshift space resummed spectrum

real PTNSeft(real k, double barr[], real U, int a, double ds1) const;
real PTNSspt(real k, double barr[], real u1, int e) const;
real PTNSspt_mat(real k, real bl, real u1, int e) const ;
real ABCnow(real k, real bl, real U,  real anw, int a, int c, int e) const;
real NL_mat(real k, real bl, real U, int c, int e) const;
real Lag_bias_nw(int a, real k, real bias[]) const;

// decomposition for GC-EFT project
real Lag_bias_decomp(real k, int a, int c) const;
real TNS_decomp(real k, int a, int c) const;

// Splined model multipoles
real pkmu_eft(double k, int a, int c, double params[]) const;
real pkmu_tns(double k, int a, int c, double params[]) const;
real Model_init(int wig) const;

private:
    const Cosmology& C;
    const PowerSpectrum& P_l;
    real epsrel;

};

#endif // EFT_H
