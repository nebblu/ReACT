#ifndef HALO_H
#define HALO_H
#include "Common.h"


// limits on mass exponents
// The maximum mass is chosen sufficiently high so that we get convergence of the 1-halo term
// The minimum mass is chosen so that we always find a solution for nu = delta_nl/sigma (M) = 1 used in the virial concentration (see  eq 46 of 1812.05594)
// The current values were chosen to satisfy these properties for a wide range of test cases in f(R) gravity where nu is mass dependent.
// Ideally you want as small a range as possible so that you are able to solve spherical collapse less often (mass_loop parameter in cosmosis pipeline)
const double Mmax = 20.;
const double Mmin = 5.;


// For massless neutrinos or no neutrinos, user should set P_cb = P_nu = P_l (total matter spectrum)
class HALO {
public:
  HALO(const Cosmology& C, const PowerSpectrum& P_l, const PowerSpectrum& P_cb,const PowerSpectrum& P_nu,const PowerSpectrum& P_cbl, real epsrel = 1e-4);

// initialise spherical collapse quantities
      int scol_init(double vars[], bool mgcamb = false , int model = 1) const;
      int scol_initp(double vars[], bool mgcamb = false) const;

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
      void react_init(double vars[], bool modg = true, int model = 1) const;
      void react_init2(double vars[],Spline ploopr, Spline ploopp, bool modg = true) const;
      double reaction(double k, double vars[]) const;

// reactions with massive neutrinos
      void react_init_nu(double vars[], bool mgcamb = false, bool modg = true, int model = 1) const;
      double reaction_nu(double k, double vars[]) const;

// Multiple redshift intialisation for cosmosis - work in progress
//  void react_init_nu2(double vars[], Spline ploopr, Spline ploopp, bool mgcamb = false, bool modg = true) const;


// 1-loop SPT reaction
      double reaction_spt(double k0, double vars[], bool mgcamb = false, int model = 1) const;

// Initialise everything for single redshift
      void initialise(double vars[], bool  mgcamb = false, bool modg = true, int model = 1) const;

// Linear spectrum for CosmoSIS
      double plinear_cosmosis(double k) const;

// Linear growth
      double Lin_Grow(double k) const;

// halofit pseudo spectrum
      double PHALO_pseudo(double k, bool mgcamb = false) const;
      // initialiser for halofit quantities - vars is as in all other functions, only call once for all k but at fixed scale factor(a = vars[0])
      void phinit_pseudo(double vars[], bool mgcamb = false)const;

// extras
      // density profile in real and fourier space (see expression in .cpp to edit to desired profile )
            double halo_profileR(double Mvir, double Rvir,  double mycvir, double r) const;
      // FT of above profile
            double halo_profileK(double k, double Mvir, double Rvir, double mycvir) const; // used for general profiles rho(r)


private:
    const Cosmology& C;
    const PowerSpectrum& P_l;
    const PowerSpectrum& P_cb;
    const PowerSpectrum& P_nu;
    const PowerSpectrum& P_cbl;
    real epsrel;


} ;



#endif // HALO
