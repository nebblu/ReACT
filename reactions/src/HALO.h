#ifndef HALO_H
#define HALO_H
#include "Common.h"


// limits on mass exponents
const double Mmax = 18.;
const double Mmin = 5.;

class HALO {
public:
  HALO(const Cosmology& C, const PowerSpectrum& P_l, real epsrel = 1e-4);


// initialise spherical collapse quantities
      int scol_init(double vars[]) const;
      int scol_initp(double vars[]) const;

// halo model components
      double rvirial(double Mvir, double vars[]) const;
      double cvirial(double Mvir, double acol) const;
      double nvirial(double Mvir, double omega0) const;

// pseudo components
      double rvirialp(double Mvir, double vars[]) const;
      double cvirialp(double Mvir, double acol) const;
      double nvirialp(double Mvir, double omega0) const;

// standard FT of NFW
      double halo_profileK2(double k,double Rvir, double mycvir) const; // standard FT of NFW

// halo model power spectra terms
      double one_halo(double k, double vars[]) const;
      double one_halop(double k, double vars[]) const;

// reactions
      void react_init(double vars[]) const;
      void react_init2(double vars[],Spline ploopr, Spline ploopp) const;
      double reaction(double k, double vars[]) const;

// Linear spectrum for CosmoSIS
      double plinear_cosmosis(double k) const;
// Linear growth 
      double Lin_Grow(double k) const;
  
// halofit pseudo spectrum
      double PHALO_pseudo(double k) const;
      // initialiser for halofit quantities - vars is as in all other functions, only call once for all k but at fixed scale factor(a = vars[0])
      void phinit_pseudo(double vars[])const;


// extras
      // density profile in real and fourier space (see expression in .cpp to edit to desired profile )
            double halo_profileR(double Mvir, double Rvir,  double mycvir, double r) const;
      // FT of above profile
            double halo_profileK(double k, double Mvir, double Rvir, double mycvir) const; // used for general profiles rho(r)


private:
    const Cosmology& C;
    const PowerSpectrum& P_l;
    real epsrel;


} ;

#endif // HALO
