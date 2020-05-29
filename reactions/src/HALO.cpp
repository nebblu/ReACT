#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <cfloat>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include <math.h>       /* pow */
#include <chrono>


#include "Common.h"
#include "Quadrature.h"
#include "HALO.h"
#include "SCOL.h"
#include "Spline.h"
#include "SpecialFunctions.h"
#include "PowerSpectrum.h"
#include "SPT.h"
#include "BSPT.h"
#include "BSPTN.h"


#include <stdio.h>
#include <gsl/gsl_sf.h>       /* pow */
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_math.h>

#include <random>
#include <functional>
#include<omp.h>

using std::cref;
using std::bind;

extern const int ERROR_MESSAGE_LEN_HALO = 512;
char error_message_halo[ERROR_MESSAGE_LEN_HALO];

void react_error_halo(const char* msg) {
    const int truncate = ERROR_MESSAGE_LEN_HALO-10;

    #pragma omp critical
    {
        snprintf(error_message_halo, truncate, "HALO.cpp error: %s", msg);
        printf("%s\n", error_message_halo);
    }
}

/* Code to calculate halo terms */

/*Load the linear power spectrum */
HALO::HALO(const Cosmology& C_, const PowerSpectrum& P_l_, real epsrel_)
: C(C_), P_l(P_l_)
{
    epsrel = epsrel_;
}


// Initial time
const double ainitial = 0.00003;

// Window function (top hat)
static real W(real x) {
    if(x < 1e-5)
        return 1. - (1./30.)*x*x;
    else
        return 3./pow3(x) * (sin(x) - x*cos(x));
}

// Window function derivative wrt R (top hat)
static real Wd(double k, double R) {
  double kR = k*R;
  double skr = sin(kR);
  double ckr = cos(kR);
        return (3.*skr - 9.*(skr-kR*ckr)/pow2(kR))/kR/R ;
}


// sigma^2 (EQ 42 of 1812.05594) in terms of log_10(mass)
// --- LCDM @ z=0 value --
// include scale factor later
static real sigma_integrand(const PowerSpectrum& P, real R, real k) {
    return k*k/(2.*M_PI*M_PI) * P(k) * pow2(W(k*R)*Dl_spt/dnorm_spt);
}


// d(sigma_8^2)/dR 
static real sigma8d_integrand(const PowerSpectrum& P, real R, real k) {
    return k*k/(2.*M_PI*M_PI) * P(k) * pow2(Dl_spt/dnorm_spt) * 2.*W(k*R)*Wd(k,R) ;
}

///// Modified Spherical Collapse  ////////

/* Initialise splines for delta_col , a_vir, delta_avir  and nu = delta_col/sigma as functions of mass index (10^mass index) */
Spline delta_col;
Spline a_vir;
Spline delta_avir;
Spline mysig;
Spline mysigp;
int HALO::scol_init(double vars[]) const{
  SCOL scol;

  int status = 0;
  // number of points in mass loop (default is 30)
  int loop_N = (int)vars[5];

  // selects whether we need a guess for the SC initial delta based on y_env
  // 0 means no guess is made and default value 3e-5 is used, 1 sets the guess to 10% larger than y_env initial.
  int yenvf = 1;

// arrays to store values
double lgmass[loop_N],sigar[loop_N],scolar0[loop_N],scolar1[loop_N],scolar2[loop_N];

// calculate sigma_8^2 and it's smoothing scale derivative @ R=8.
double sig1 = Integrate<ExpSub>(bind(sigma_integrand, cref(P_l), 8., std::placeholders::_1), 1e-4, 50., 1e-3);
double sig2 = Integrate<ExpSub>(bind(sigma8d_integrand, cref(P_l), 8., std::placeholders::_1), 1e-4, 50., 1e-3);

 if (!gsl_finite(sig1)) {
  react_error_halo("sigma_8 evaluated to non-numerical value");
}

if (!gsl_finite(sig2)) {
  react_error_halo("sigma_8 derivative evaluated to non-numerical value");
}

// theory parameters
double pars[3];
pars[0] = vars[2];
pars[1] = vars[3];
pars[2] = vars[4];

// log(mass) loop for calculating the collapse quantities
//#pragma omp parallel for
for(int i = 0; i< loop_N; i++){
   double myscolparams[3];

// log mass
   lgmass[i] = Mmin + i*(Mmax-Mmin)/(loop_N-1.); // linear spacing

// top hat radius
   double Rth = 0.1*pow((Gnewton*pow(10, lgmass[i]))/(5.*vars[1]),ONE/THREE);

// initialise delta_c, a_vir, delta_vir
   scol.myscol(myscolparams, vars[0], vars[1], Rth, sig1, sig2, pars, 2, yenvf);
   if(scol.error.errorno != 0) {
    status = scol.error.errorno;
   }

// calculate sigma^2
   sigar[i] = sqrt(Integrate<ExpSub>(bind(sigma_integrand, cref(P_l), Rth, std::placeholders::_1), 1e-4, 50., 1e-5, 1e-5));

// store values in arrays
   scolar0[i] = myscolparams[0];
   scolar1[i] = myscolparams[1];
   scolar2[i] = myscolparams[2];

}


// Spline over arrays
// sigma^2
  mysig = CubicSpline(loop_N,lgmass,sigar);
// delta_c
  delta_col = CubicSpline(loop_N,lgmass,scolar0);
// a_virial
  a_vir = CubicSpline(loop_N,lgmass,scolar1);
// delta_virial
  delta_avir = CubicSpline(loop_N,lgmass,scolar2);

  return status;

}



///// Pseudo spherical collapse  ////////

Spline delta_colp;
Spline a_virp;
Spline delta_avirp;
Spline linear_growth;

// sigma^2 for pseudo cosmology --- modified
static real sigma_integrandp(const PowerSpectrum& P, real R, real k) {
    return k*k/(2.*M_PI*M_PI) * P(k) * pow2(linear_growth(k)*W(k*R));
}

// d(sigma_8^2)/dR
static real sigma8d_integrandp(const PowerSpectrum& P, real R, real k) {
    return k*k/(2.*M_PI*M_PI) * P(k) * pow2(linear_growth(k)) * 2.*W(k*R)*Wd(k,R) ;
}


int HALO::scol_initp(double vars[]) const{
  SCOL scol;
  IOW iow;

  int status = 0;
  int loop_N = 30;
  int loop_Nk = 300;

    // arrays to store values
  double myscolparams[3],lgmass[loop_N],sigar[loop_N],scolar0[loop_N],scolar1[loop_N],scolar2[loop_N];
  double kval_tab[loop_Nk],ling_tab[loop_Nk];

  // initialise linear growth spline for the theory of gravity/dark energy. See specialfunctions.cpp mu,gamma2,gamma3 functions
  double kmin = 1e-5;
  double kmax = 100.;

  for(int i = 0; i< loop_Nk; i++){
    kval_tab[i] =  kmin * exp(i*log(kmax/kmin)/(loop_Nk-1.));
    // numerical initialisation of linear growth in general theory of gravity/dark energy - see SpecialFunctions.cpp
    iow.initn_lin(vars[0], kval_tab[i], vars[1],vars[2], vars[3], vars[4]);
    ling_tab[i] = F1_nk/dnorm_spt;
      }

// spline linear growth factor
   linear_growth = LinearSpline(loop_Nk,kval_tab,ling_tab);


// calculate sigma_8^2 and it's smoothing scale derivative @ R=8.
double sig1 = Integrate<ExpSub>(bind(sigma_integrandp, cref(P_l), 8., std::placeholders::_1), 1e-4, 50., 1e-3);
double sig2 = Integrate<ExpSub>(bind(sigma8d_integrandp, cref(P_l), 8., std::placeholders::_1), 1e-4, 50., 1e-3);

    
   // theory params
   double pars[4];
   pars[0] = 1e-15;
   pars[1] = 1.;
   pars[2] = 1.;
   // y env flag
   int yenvf = 0;

// initialise spherical collapse quantities in GR (independent of R (or M))
  scol.myscol(myscolparams, vars[0], vars[1], 1., sig1, sig2, pars, 1, 0);
  if(scol.error.errorno != 0) {
    status = scol.error.errorno;
  }

//#pragma omp parallel for
for(int i = 0; i< loop_N; i++){
   lgmass[i] = Mmin + i*(Mmax-Mmin)/(loop_N-1.); // linear spacing

   double Rth = 0.1*pow((Gnewton*pow(10, lgmass[i]))/(5.*vars[1]),ONE/THREE);

// calculate sigma^2
    sigar[i] = sqrt(Integrate<ExpSub>(bind(sigma_integrandp, cref(P_l), Rth, std::placeholders::_1), 1e-4, 50., 1e-5, 1e-5));

// store values in arrays
    scolar0[i] = myscolparams[0];
    scolar1[i] = myscolparams[1];
    scolar2[i] = myscolparams[2];
      }

// Spline over arrays
  // sigma^2
    mysigp = CubicSpline(loop_N,lgmass,sigar);
  // delta_c
    delta_colp = CubicSpline(loop_N,lgmass,scolar0);
  // a_virial
    a_virp = CubicSpline(loop_N,lgmass,scolar1);
  // delta_virial
    delta_avirp = CubicSpline(loop_N,lgmass,scolar2);

    return status;
}


/* HALO MODEL QUANTITIES */

//  nu and dnu/dlnM , a being the scale factor,  used in mass functions

// real
static double mynu(double lgmass, void * params){
  (void)(params);
  return  delta_col(lgmass)/mysig(lgmass);
}
static double nu_der(double lgmass){
  gsl_function F;
  double result, abserr;
  F.params = 0;
  F.function = &mynu;
  gsl_deriv_central(&F, lgmass, 1e-5, &result, &abserr);
return result;
}

// pseudo
static double mynup(double lgmass, void * params){
  (void)(params);
  return delta_colp(lgmass)/mysigp(lgmass);
}
static double nu_derp(double lgmass){
  gsl_function F;
  double result, abserr;
  F.params = 0;
  F.function = &mynup;
  gsl_deriv_central(&F, lgmass, 1e-5, &result, &abserr);
return result;
}



////////////////////////////////////////////////////////////////////////////////////////


/* Virial radius */

// vars:
// 0 = acol, 1= omega0, 2 = mg param

// real
double HALO::rvirial(double lgmass, double vars[]) const {
      //virial overdensity
      double Deltavir = (1.+ delta_avir(lgmass))*pow3(vars[0]/a_vir(lgmass));
      return 0.1*pow((pow(10.,lgmass)*Gnewton)/(5.*vars[1]*Deltavir), 1./3.);
    }

//pseudo
double HALO::rvirialp(double lgmass, double vars[]) const {
      //virial overdensity
      double Deltavir = (1.+ delta_avirp(lgmass))*pow3(vars[0]/a_virp(lgmass));
      return 0.1*pow((pow(10.,lgmass)*Gnewton)/(5.*vars[1]*Deltavir), 1./3.);
    }


////////////////////////////////////////////////////////////////////////////////////////


/* Virial concentration */

// real
double HALO::cvirial(double lgmass, double acol) const {

      double myc0 = 9.;
      double alpha = 0.13;

      double lgmstar;

      const double upper = Mmax;
      const double lower = Mmin;


      std::mt19937 gen(2);

      std::uniform_real_distribution<> dis(lower, upper);

      double pos_pt = dis(gen);
      double neg_pt = dis(gen);

      while ( delta_col(lgmass)/mysig(pos_pt)-1.  < 0.){
          pos_pt = dis(gen);}

      auto t1 = std::chrono::high_resolution_clock::now();
     
      while ( delta_col(lgmass)/mysig(neg_pt)-1. > 0.){
          neg_pt = dis(gen);
      
          auto t2 = std::chrono::high_resolution_clock::now();
          double duration1 = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
          duration1 *= 1e-9;
       
          if (duration1>1e-4) {
          return g_de*myc0*pow(10.,-alpha*(lgmass-Mmin))*acol;
                        } 
      
      }

       const double about_zero_mag = 1e-3;
      for (;;)
      {
          double mid_pt = (pos_pt + neg_pt)/2.0;
          double f_mid_pt = delta_col(lgmass)/mysig(mid_pt)-1.;

          if (fabs(f_mid_pt)  < about_zero_mag){
              lgmstar = mid_pt;
              break;
                }

          if (f_mid_pt >= 0.){
              pos_pt = mid_pt;}
          else{
              neg_pt = mid_pt;}
      }

    return g_de*myc0*pow(10.,-alpha*(lgmass-lgmstar))*acol;

  }


//pseudo

  double HALO::cvirialp(double lgmass, double acol) const {

        double myc0 = 9.;
        double alpha = 0.13;
        double lgmstar;

        const double upper = Mmax;
        const double lower = Mmin;


        std::mt19937 gen(2);

        std::uniform_real_distribution<> dis(lower, upper);

        double pos_pt = dis(gen);
        double neg_pt = dis(gen);

        while ( delta_colp(lgmass)/mysigp(pos_pt)-1.  < 0.){
            pos_pt = dis(gen);}

        while ( delta_colp(lgmass)/mysigp(neg_pt)-1. > 0.){
            neg_pt = dis(gen);}

         const double about_zero_mag = 1e-3;
        for (;;)
        {
            double mid_pt = (pos_pt + neg_pt)/2.0;
            double f_mid_pt = delta_colp(lgmass)/mysigp(mid_pt)-1.;

            if (fabs(f_mid_pt)  < about_zero_mag){
                lgmstar = mid_pt;
                break;
                  }

            if (f_mid_pt >= 0.){
                pos_pt = mid_pt;}
            else{
                neg_pt = mid_pt;}
        }

      return myc0*pow(10.,-alpha*(lgmass-lgmstar))*acol;

    }

    ////////////////////////////////////////////////////////////////////////////////////////

/* Virial number density */

// real
  double HALO::nvirial(double lgmass, double omega0) const {
       void *params;
       double q = 0.75;
       double p = 0.3;
       double aconst = 1./3.056;


       double v2 = pow2(mynu(lgmass,params));
       double vfv = aconst*sqrt(2./M_PI * q * v2 ) * (1.+ pow(v2*q,-p)) * exp(-q*v2/2.);
       double myrho =  15.* omega0/(4.*M_PI*Gnewton*pow3(0.1));
       double Mvir = pow(10.,lgmass);

       double nuprime = 0.434294*nu_der(lgmass); // log10 (Mass) derivative of  nu, constant is conversion from ln to log10 derivative (=log_10(e))

       return myrho * vfv/Mvir *  1./ mynu(lgmass,params) * nuprime;    /// Halo mass function (ST for now)
  }

// pseudo
  double HALO::nvirialp(double lgmass, double omega0) const {
       void *params;
       double q = 0.75;
       double p = 0.3;
       double aconst = 1./3.056;


       double v2 = pow2(mynup(lgmass,params));
       double vfv = aconst*sqrt(2./M_PI * q * v2 ) * (1.+ pow(v2*q,-p)) * exp(-q*v2/2.);
       double myrho =  15.* omega0/(4.*M_PI*Gnewton*pow3(0.1));
       double Mvir = pow(10.,lgmass);

       double nuprime = 0.434294*nu_derp(lgmass); // log10 (Mass) derivative of  nu, constant is conversion from ln to log10 derivative (=log_10(e))

       return myrho * vfv/Mvir *  1./ mynup(lgmass,params) * nuprime;    /// Halo mass function (ST for now)

  }


  ////////////////////////////////////////////////////////////////////////////////////////

// STANDARD FT OF NFW PROFILE (SEE EQ. 81 OF 0206508)
double HALO::halo_profileK2(double k, double Rvir, double mycvir) const {
     double rs = Rvir/mycvir;
     double krs = k*rs;
     double csp = 1.+mycvir;
     double term1 = 1./(log(csp) - mycvir/csp);

     return term1  * (sin(krs) * (gsl_sf_Si(csp*krs) - gsl_sf_Si(krs))
                    - sin(krs*mycvir)/(csp*krs)
                    + cos(krs) * (gsl_sf_Ci(csp*krs) - gsl_sf_Ci(krs)));
}



// 1 - halo term integrand
//vars:
// 0 = acol, 1= omega0, 2 = mg param
static double halo_integrand(const Cosmology& C, const PowerSpectrum& P_l, double vars[], double k, double lgmvir){
        HALO halo(C, P_l, 1e-3);

        // background matter density

        double myrho = 15.*vars[1]/(4.*M_PI*Gnewton*pow3(0.1));

        // Halo density profile components //

        // virial radius
        double Rvir = halo.rvirial(lgmvir, vars); // run first to initialise delta_collapse used in c_vir

        // virial concentration
        double mycvir = halo.cvirial(lgmvir,vars[0]);

        // virial mass
        double Mvir = pow(10,lgmvir);

        // initialise normalisation of fourier transform of halo density profile (only for general,for NFW not needed)
        //profn_init(cref(C), cref(P_l),  Mvir,  Rvir,  mycvir);

        // halo density profile
        double hpk = fabs(halo.halo_profileK2(k,Rvir,mycvir));
        double halof2 = pow2(hpk);

        // virial number density
        double mynvir = halo.nvirial(lgmvir, vars[1]);

        double prefac = pow2(Mvir/myrho);

        return prefac * mynvir * halof2;

}


// real one halo term

// vars:
// 0 = acol, 1= omega0, 2 = mg param
double HALO::one_halo(double k, double vars[]) const {
    // constant converts from an integral over ln M to log10 M
    return  Integrate(bind(halo_integrand, cref(C), cref(P_l), vars, k, std::placeholders::_1),  1.01*Mmin , 0.99*Mmax  , epsrel)/0.434294;
  }



// 1 - halo term integrand for pseudo
// vars:
// 0 = acol, 1= omega0, 2 = mg param
  static double halo_integrandp(const Cosmology& C, const PowerSpectrum& P_l, double vars[], double k, double lgmvir){
          HALO halo(C, P_l, 1e-3);

          // background matter density
          double myrho = 15.*vars[1]/(4.*M_PI*Gnewton*pow3(0.1));

          // Halo density profile components //

          // virial radius
          double Rvir = halo.rvirialp(lgmvir, vars); // run first to initialise delta_collapse used in c_vir

          // virial concentration
          double mycvir = halo.cvirialp(lgmvir,vars[0]);

          // virial mass
          double Mvir = pow(10,lgmvir);

          // initialise normalisation of fourier transform of halo density profile (only for general,for NFW not needed)
          //profn_init(cref(C), cref(P_l),  Mvir,  Rvir,  mycvir);

          // halo density profile
          double hpk = halo.halo_profileK2(k, Rvir,mycvir);
          double halof2 = pow2(hpk);


          // virial number density
          double mynvir = halo.nvirialp(lgmvir, vars[1]);

          double prefac = pow2(Mvir/myrho);

          return prefac * mynvir * halof2;

  }

// pseudo one halo term
  double HALO::one_halop(double k, double vars[]) const {
      // constant converts from an integral over ln M to log10 M
      return  Integrate(bind(halo_integrandp, cref(C), cref(P_l), vars, k, std::placeholders::_1), 1.01*Mmin ,0.99*Mmax , epsrel)/0.434294;
    }



////////////// REACTIONS //////////////

/// reaction parameter intiialisation for a single redshift

// kstar, mathcalE and linear power spectrum spline
double kstar, bigE;

// vars:
// 0 = acol, 1= omega0, 2 = mg param
void HALO::react_init(double vars[]) const{
  SPT spt(C, P_l, epsrel);

  double pspt, psptp, p1h, p1hp, plreal,argument;

/// mathcal E///

  bigE = one_halo(0.01,vars)/one_halop(0.01,vars);

  /// Calculate kstar///

// 1halo terms

  p1h = one_halo(0.06,vars); // real
  p1hp = one_halop(0.06,vars); // pseudo

// spt terms
// Real PT spectrum
  pspt = spt.PLOOPn2(1, vars, 0.06, 1e-3);

  // GR-1-loop spectrum with linear growth replaced by modified gravity growth (unscreened)
  psptp = spt.PLOOPn2(8, vars, 0.06, 1e-3);

  // real linear spectrum
  plreal = pow2(linear_growth(0.06))*P_l(0.06);

  argument = ( (((pspt + p1h)/(psptp + p1hp)) * (plreal + p1hp) - p1h ) / plreal  - bigE )/(1.-bigE);

  if(argument<0.001 || argument>1.){
    kstar = 1e-6;
    bigE = 1.;
  }
  else{
  kstar = (-0.06/log(argument));
  }

}


/// reaction parameter initialisation for multiple redshifts --- useful for lensing --  much faster than running react_init multiple times

// Reaction using 1-loop splines (ploopr and ploopp) as a function of redshift -- these should be initialised of course, one can use the ploop_init function in SPT.cpp
void HALO::react_init2(double vars[], Spline ploopr, Spline ploopp) const{
  SPT spt(C, P_l, epsrel);

  double pspt, psptp, p1h, p1hp, plreal,argument;

/// mathcal E///

  bigE = one_halo(0.01,vars)/one_halop(0.01,vars);

  /// Calculate kstar///

// 1halo terms
  p1h = one_halo(0.06,vars); // real
  p1hp = one_halop(0.06,vars); // pseudo

// spt terms
// Real PT spectrum
  pspt = ploopr(vars[0]);
// GR-1-loop spectrum with linear growth replaced by modified gravity growth (unscreened)
  psptp = ploopp(vars[0]);

  // real linear spectrum
  plreal = pow2(linear_growth(0.06))*P_l(0.06);

  argument = ( (((pspt + p1h)/(psptp + p1hp)) * (plreal + p1hp) - p1h ) / plreal  - bigE )/(1.-bigE);

  if(argument<0.001 || argument>1.){
    kstar = 1e-6;
    bigE = 1.;
  }
  else{
  kstar = (-0.06/log(argument));
  }

}

// Compute the reaction (Eq.53 of 1812.05594 )

// vars is 0: scale factor, 1: omega_matter today , 2: modified gravity param
double HALO::reaction(double k, double vars[]) const {
  // real linear spectrum
 double plreal = pow2(linear_growth(k))*P_l(k);
 return (((1.-bigE)*exp(-k/kstar) + bigE)*plreal+  one_halo(k,vars))/(plreal+ one_halop(k,vars));
}


// Compute linear power spectrum to be handed to CosmoSIS -  linear_growth must be initialised with scol_initp (needed for pseudo halo model quantities)
double HALO::plinear_cosmosis(double k) const {
  // real linear spectrum
  return pow2(linear_growth(k))*P_l(k);
}


// Pseudo halofit prescription
////////////////////////////////////////////////
///////HALO-FIT FOR NON-SCALE DEP MODELS : takashi et al 2012 fits  1208.2701///////////
///////////////////////////////////////////////

double wintcambs_pseudo[3];
static double wintcamb_pseudo(const PowerSpectrum& P_L, double r){
 int n1 = 3000;
 double sum1 = 0.;
 double sum2 = 0.;
 double sum3 = 0.;
 double anorm = 1./(2.*pow2(M_PI));
 double t,k, pkterm,x,x2,w1,w2,w3,fac;
 for(int i = 1; i<=n1 ; i++){
   t = (i-0.5)/n1;
   k = -1. + 1./t;
   pkterm = pow2(linear_growth(k))*P_L(k)*pow3(k)*anorm;
   x = k*r;
   x2 = pow2(x);
   w1 = exp(-x2);
   w2 = 2.*x2*w1;
   w3 = 4.*x2*(1.-x2)*w1;
   fac = pkterm/k/t/t;
   sum1 += w1*fac;
   sum2 += w2*fac;
   sum3 += w3*fac;
 }
sum1 /= n1;
sum2 /= n1;
sum3 /= n1;
wintcambs_pseudo[0] = sqrt(sum1);
wintcambs_pseudo[1] = -sum2/sum1;
wintcambs_pseudo[2] = -sum2*sum2/sum1/sum1 - sum3/sum1;
}

double phpars_pseudo[13];
void parscamb_pseudo(const PowerSpectrum& P_L)
{
  double xlogr1 = -2.;
  double xlogr2 = 3.5;
  double rmid = 0.;
  double diff = 1.;
  for (;;){
    rmid = (xlogr2 + xlogr1)/2.;
    rmid = pow(10,rmid);
    wintcamb_pseudo(cref(P_L),rmid);
    diff = wintcambs_pseudo[0]-1.;
    if (fabs(diff)<=0.0001) {
      phpars_pseudo[0] = 1./rmid;
      phpars_pseudo[11] = -3.-wintcambs_pseudo[1];
      phpars_pseudo[12] = -wintcambs_pseudo[2];
      break;
    }
   else if(diff>0.0001){
     xlogr1 = log10(rmid);
   }
   else if (diff < -0.0001){
     xlogr2 = log10(rmid);
   }
    if (xlogr2<-1.9999){
      phpars_pseudo[0] = 1000.;
      break;
    }
    else if(xlogr2>3.4999){
      phpars_pseudo[0] = 1000.;
      break;
    }
  }
}

// Initialises components
void HALO::phinit_pseudo(double vars[])const{
    double scalef = vars[0];
    double omega0 = vars[1];
    parscamb_pseudo(cref(P_l));

    double neff = phpars_pseudo[11];
    double curv = phpars_pseudo[12];

   double neff2 = pow2(neff);
   double neff3 = neff2*neff;
   double neff4 = pow2(neff2);


   double ane = 1.5222 + 2.8553*neff + 2.3706*neff2+ 0.9903*neff3 + 0.2250*neff4- 0.6038*curv;
   double bne = -0.5642 + 0.5864*neff + 0.5716*neff2- 1.5474*curv;
   double cne = 0.3698 + 2.0404*neff + 0.8161*neff2 + 0.5869*curv;

   double nue = 5.2105 + 3.6902*neff;
   double acub = pow3(scalef);
   double omgm = (omega0/acub)/(omega0/acub + (1.-omega0));
   double omgl = (1.-omega0)/(omega0/acub + (1.-omega0));

     phpars_pseudo[1] = pow(10.,ane); // an
     phpars_pseudo[2] = pow(10.,bne); // bn
     phpars_pseudo[3] = pow(10.,cne); //cn

     phpars_pseudo[4] = 0.1971 - 0.0843*neff + 0.8460*curv; //gan
     phpars_pseudo[5] = fabs(6.0835 + 1.3373*neff - 0.1959*neff2- 5.5274*curv); // alpha
     phpars_pseudo[6] = 2.0379 - 0.7354*neff + 0.3157*neff2 + 1.2490*neff3 + 0.3980*neff4- 0.1682*curv; //beta

     phpars_pseudo[7] = pow(10.,nue); // nun

     if (fabs(1-omgm) > 0.01) {
       double f1a = pow(omgm,-0.0732);
       double f2a = pow(omgm,-0.1423);
       double f3a = pow(omgm,0.0725);

       double f1b = pow(omgm,-0.0307);
       double f2b = pow(omgm,-0.0585);
       double f3b = pow(omgm,0.0743);

        double frac = omgl/(1. - omgm);
       phpars_pseudo[8] = frac*f1b + (1.-frac)*f1a;
       phpars_pseudo[9] = frac*f2b + (1.-frac)*f2a;
       phpars_pseudo[10] = frac*f3b + (1.-frac)*f3a;

        }
      else{
        phpars_pseudo[8]=1.;
        phpars_pseudo[9]=1.;
        phpars_pseudo[10]=1.;
      }

  }


double HALO::PHALO_pseudo(double k) const{
  IOW iow;
  if (phpars_pseudo[0] == 1000. || k<=0.005) {
    return pow2(linear_growth(k))*P_l(k);
  }
  else{
  double kcub = pow3(k)/2./pow2(M_PI);
  double deltahp = phpars_pseudo[1]*pow(k/phpars_pseudo[0],phpars_pseudo[8]*3.)/(1.+phpars_pseudo[2]*pow(k/phpars_pseudo[0],phpars_pseudo[9])+pow(phpars_pseudo[3]*phpars_pseudo[10]*k/phpars_pseudo[0],3.-phpars_pseudo[4]));
  double deltah = deltahp/(1.+ phpars_pseudo[7]/pow2(k/phpars_pseudo[0]));
  double deltaplin = kcub*pow2(linear_growth(k))*P_l(k);
  double deltaq = deltaplin*pow(1.+ deltaplin,phpars_pseudo[6])/(1.+phpars_pseudo[5]*deltaplin)*exp(-(k/phpars_pseudo[0])/4.-pow2(k/phpars_pseudo[0])/8.);
  double deltanl = deltaq + deltah;

  return  deltanl/kcub;
}
}



/* Extra stuff */

/// Halo density profile in real space - NFW for now  but  can edit this easily
double HALO::halo_profileR(double Mvir, double Rvir,  double mycvir, double r) const {
  double rs = Rvir/mycvir;
  double rhos = Mvir/(4.*M_PI*pow3(rs))/(log(1.+mycvir) - mycvir/(1.+mycvir));
  // NFW
  return  rhos/(r/rs * pow2(1.+ r/rs));
}

// profile in real space
static double halo_profile_integrand(const Cosmology& C, const PowerSpectrum& P_l, double Mvir, double Rvir, double mycvir, double k, double r){
      HALO halo(C, P_l, 1e-3);
      return  sin(k*r)*r/k*halo.halo_profileR(Mvir, Rvir, mycvir,r);
  }

/// Halo density profile in k space

// normalise it
double profile_norm = 1.;

void profn_init(const Cosmology& C, const PowerSpectrum& P_l, double Mvir, double Rvir, double mycvir){
  HALO halo(C, P_l, 1e-3);
  profile_norm =  4.*M_PI*Integrate(bind(halo_profile_integrand, cref(C), cref(P_l), Mvir, Rvir, mycvir, 1e-4, std::placeholders::_1), QMINp, Rvir , 1e-3);
}

// ONE CAN USE THIS FUNCTION IF YOU WANT TO EDIT NFW PROFILE USED IN COMPUTATIONS
double HALO::halo_profileK(double k, double Mvir, double Rvir, double mycvir) const {
     profn_init(cref(C), cref(P_l),  Mvir,  Rvir,  mycvir);
     return 4.*M_PI/profile_norm*Integrate(bind(halo_profile_integrand, cref(C), cref(P_l), Mvir, Rvir, mycvir, k, std::placeholders::_1), QMINp, Rvir , epsrel);

}
