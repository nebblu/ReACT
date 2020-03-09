#if HAVE_CONFIG_H
# include <config.h>
#endif


#include <boost/bind.hpp>
using boost::cref;


#include "Cosmology.h"
#include "PowerSpectrum.h"
#include "Quadrature.h"
#include "BSPT.h"
#include "BSPTN.h"
#include "SPT.h"
#include "RegPT.h"

#include "Spline.h"
#include "SpecialFunctions.h"
#include "LinearPS.h"
#include "CMBc.h"


#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include <math.h>       /* pow */
#include <boost/math/tools/roots.hpp>
#include <random>


#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>


// bottom limit of comoving distance integral
const double xim = 10.;

CMBc::CMBc(const Cosmology& C_, const PowerSpectrum& P_L_, real epsrel_)
: C(C_), P_L(P_L_)
{
    epsrel = epsrel_;
}

//////////////////////////////////////////
 /* Functions used in initialisation */
//////////////////////////////////////////

// Hubble function
double HAc(double omega0, double z){
	double omegaL= 1.-omega0;
	return  1./sqrt(omega0*pow(1.+z,3)+omegaL);}

// Functions for GM or SC formula
// neff(k)
  static double myfuncnw(double k,void * params){
    (void)(params);
    return log(neffnw(k));
  }

  static double myfuncehu(double k,void * params){
    (void)(params);
    return log(neffehu(k));
  }

  static double myfunclin(double k,void * params){
    (void)(params);
    return log(mylinearps(k));
  }

static double myneff(double k){
  gsl_function F;
  double result, abserr;
    F.params = 0;
    F.function = &myfuncnw;
    gsl_deriv_central(&F, k, 1e-6, &result, &abserr);
return k*result;
}


// Find K_nl
  double myknls(const PowerSpectrum& P_L)
  {
      std::mt19937 gen(2);
       const double lower_bound =  1e-4;
       const double upper_bound =  1000.;
      std::uniform_real_distribution<> dis(lower_bound, upper_bound);

      double pos_pt = dis(gen);
      double neg_pt = dis(gen);

      while ( (pow2(D_spt/dnorm_spt)*P_L(pos_pt)*pow3(pos_pt)/2./M_PI/M_PI  - 1.) < 0.){
          pos_pt = dis(gen);}

      while ( (pow2(D_spt/dnorm_spt)*P_L(neg_pt)*pow3(neg_pt)/2./M_PI/M_PI- 1.) > 0.){
          neg_pt = dis(gen);}

       const double about_zero_mag = 1e-5;
      for (;;)
      {
          double mid_pt = (pos_pt + neg_pt)/2.0;
          double f_mid_pt =  pow2(D_spt/dnorm_spt)*P_L(mid_pt)*pow3(mid_pt)/2./M_PI/M_PI -1.;

          if (fabs(f_mid_pt)  < about_zero_mag){
              return mid_pt;
                }

          if (f_mid_pt >= 0.){
              pos_pt = mid_pt;}
          else{
              neg_pt = mid_pt;}

      }
  }

// regpt damp term for post born prediction
static real regpt_exp(const PowerSpectrum& P_l, double q){
    return  P_l(q)/(6.*pow2(M_PI));
  }

// Calculate sigma_8 at some redshift
static double sigma8int(const PowerSpectrum& P_L, double k){
  return pow2(3.*((sin(k*8.) - k*8.*cos(k*8.))/pow2(k*8.))/k/8.) * pow2(k) * pow2(D_spt/dnorm_spt) * P_L(k)/2./M_PI/M_PI;
}

static double mysigma8(const PowerSpectrum& P_L){
  return sqrt(Integrate(bind(sigma8int,cref(P_L),_1),0.0001,30.,1e-4));
}

// Spline function z->comoving distance
Spline zofc;
// Spline linear growth factor + F2_dgp
Spline Dlin, Dlingr, F2DGP;
// Splines for HALOFIT
Spline phpars0,phpars1,phpars2,phpars3,phpars4,phpars5,phpars6,phpars7,phpars8,phpars9,phpars10,phpars11,phpars12;
// Splines for dgp 1-loop matter power Spectrum
Spline fdgps, cdgps, idgps, kldgps, p22spl,p13spl, damp_termcmb;
// Comoving distance to z=1000 (~CMB)
double xis;
// comoving distance to z=40 as in TN's computation
double xis40;

// Splines for GM formula
Spline myneffspl, mys8spl, myknlspl;

// Initialize z(comoving_dist) and xi_star = comoving distance to CMB
// a = 1 linear
// a = 2 1-loop (+ initialises regpt components)
// a = 3 Gil-Marin or Scoccimaro [+ initialises all halofit params, n_eff(k), k_nl(z) and sigma_8(z)]
 void CMBc::comdist_init(int a, double omega0, double mg1){
  SPT spt(C, P_L, 1e-3);
  RegPT rpt(C, P_L, 1e-3);
  BSPT bspt(C, P_L, 1e-3);
  IOW iow;

  // initialise damping factor for regpt
  if(a==2){
  rpt.sigmad_init();
    }

   // initialise EHU PS and Fonseca NW
    if(a==3){
      bspt.mypspline(2);
    }

  double h0inv = 2997.92;
  const double zmin = 1e-4;
  const double zmax = 1100.;
  const double kmin = 1e-4;
  const double kmax = 100.;

  vector<double> ztable, ktable, ctable, dtable, dgrtable, f2table;
  // Halofit tables
  vector<double> phtable0,phtable1,phtable2,phtable3,phtable4,phtable5,phtable6,phtable7,phtable8,phtable9,phtable10,phtable11,phtable12;
  // 1-LOOP SPT nDGP GROWTH TABLES
  vector<double> ftab,ctab,itab,kltab;
  // 1-Loop k-components
  vector<double> p13tab,p22tab, dampt;
  // GM formula tables
  vector<double> nefftab, knltab, sig8tab;


  xis =  h0inv*Integrate(bind(HAc,omega0,_1),0.,zmax,1e-4);
  xis40 = h0inv*Integrate(bind(HAc,omega0,_1),0.,40.,1e-4);

  const int n1 = 1000;
  double cval,zval,scalef,kval,sigd;

  for(int i = 0; i<n1; i++ ){
    zval = zmin*exp(i*log(zmax/zmin)/(n1-1.));
    kval = kmin*exp(i*log(kmax/kmin)/(n1-1.));

    scalef = 1./(1.+zval);
    iow.inite(scalef,omega0,mg1,1.,1.);

    cval = h0inv*Integrate(bind(HAc,omega0,_1),zmin,zval,1e-4);

    ztable.push_back(zval);
    ktable.push_back(kval);
    ctable.push_back(cval);

    dtable.push_back(D_spt/dnorm_spt);
    dgrtable.push_back(Dl_spt/dnorm_spt);
    f2table.push_back(F_spt/pow2(dnorm_spt));

    if (a==2) {
        p22tab.push_back(pow4(dnorm_spt/Dl_spt)*spt.P22_dd(kval));
        p13tab.push_back(pow4(dnorm_spt/Dl_spt)*spt.P13_dd(kval));
        sigd = pow2(kval)*Integrate<ExpSub>(bind(regpt_exp, cref(P_L), _1), kmin, kval/2., 1e-4);
        dampt.push_back(sigd);
    }

    if (a ==3 ) {
      spt.phinit(scalef,omega0);
      phtable0.push_back(phpars[0]);
      phtable1.push_back(phpars[1]);
      phtable2.push_back(phpars[2]);
      phtable3.push_back(phpars[3]);
      phtable4.push_back(phpars[4]);
      phtable5.push_back(phpars[5]);
      phtable6.push_back(phpars[6]);
      phtable7.push_back(phpars[7]);
      phtable8.push_back(phpars[8]);
      phtable9.push_back(phpars[9]);
      phtable10.push_back(phpars[10]);
      phtable11.push_back(phpars[11]);
      phtable12.push_back(phpars[12]);

        if (i==0) {
          knltab.push_back(myknls(cref(P_L)));
          }
        else if( knltab[i-1] > 200. ){
           knltab.push_back(knltab[i-1]);
        }
        else{
           knltab.push_back(myknls(cref(P_L)));
              }
             printf("%d %e \n", i, knltab[i]);

      sig8tab.push_back(mysigma8(cref(P_L)));
      nefftab.push_back(myneff(kval));

    }
      }

      zofc = LinearSpline(ctable,ztable);
      Dlin = LinearSpline(ztable,dtable);
      Dlingr = LinearSpline(ztable,dgrtable);
      F2DGP = LinearSpline(ztable,f2table);

      if (a==2) {
        p22spl=LinearSpline(ktable,p22tab);
        p13spl=LinearSpline(ktable,p13tab);
        damp_termcmb = LinearSpline(ktable,dampt);
            }
      if (a==3) {
        phpars0 = LinearSpline(ztable,phtable0);
        phpars1 = LinearSpline(ztable,phtable1);
        phpars2 = LinearSpline(ztable,phtable2);
        phpars3 = LinearSpline(ztable,phtable3);
        phpars4 = LinearSpline(ztable,phtable4);
        phpars5 = LinearSpline(ztable,phtable5);
        phpars6 = LinearSpline(ztable,phtable6);
        phpars7 = LinearSpline(ztable,phtable7);
        phpars8 = LinearSpline(ztable,phtable8);
        phpars9 = LinearSpline(ztable,phtable9);
        phpars10 = LinearSpline(ztable,phtable10);
        phpars11 = LinearSpline(ztable,phtable11);
        phpars12 = LinearSpline(ztable,phtable12);

        myknlspl = LinearSpline(ztable,knltab);
        mys8spl = LinearSpline(ztable,sig8tab);
        myneffspl = LinearSpline(ktable,nefftab);

      }
}

// B_LSS integrands
// tree analytic
static double btlssa_integrand(const Cosmology& C, const PowerSpectrum& P_L, int a, double params[], double l1, double l2, double x, double xi ){
   BSPT bspt(C, P_L, 1e-4);
   IOW iow;
   double bis, kernel, F2A, F2B, F2C, growth;
   double h0 = 1./2997.92;
   double lensingk = (xis-xi)/xi/xis;
   double myz = zofc(xi);
   double scalef = 1./(1.+myz);
   double prefactor = pow3(3.*params[0]*pow2(h0)/(2.*scalef) * lensingk);
   double k1 = l1/xi;
   double k2 = l2/xi;

   D_spt = Dlin(myz) * dnorm_spt;
   Dl_spt = Dlin(myz) * dnorm_spt; //Dlingr(myz) * dnorm_spt;
   F_spt = F2DGP(myz) * dnorm_spt;
   bis = bspt.Btree(a, k1, k2, x);

   return prefactor * bis * pow2(xi);
}

// 1-loop analytic
static double bllssa_integrand(const Cosmology& C,const PowerSpectrum& P_L, int a, double params[], double l1, double l2, double x, double xi ){
   BSPT bspt(C, P_L, 1e-2);
   IOW iow;
   double bis, kernel;
   double h0 = 1./2997.92;
   double lensingk = (xis-xi)/xi/xis;
   double myz = zofc(xi);
   double scalef = 1./(1.+myz);
   double prefactor = pow3(3.*params[0]*pow2(h0)/(2.*scalef) * lensingk);
   switch (a) {
     case 1:
      Dl_spt = Dlin(myz) * dnorm_spt;
      break;
      case 2:
   iow.inite2(scalef,params[0],params[1],params[2],params[3]);
      break;
 }
   double k1 = l1/xi;
   double k2 = l2/xi;
   bis = bspt.Bloop(a, k1, k2, x);
   return prefactor * bis * pow2(xi);
}

// Gil-Marin formula
static double bgmlssa_integrand(const Cosmology& C,const PowerSpectrum& P_L,  double params[], double l1, double l2, double x, double xi ){
   SPT spt(C, P_L, 1e-3);
   double vars[6],F2A,F2B,F2C,p1,p2,p3;

   double h0 = 1./2997.92;
   double lensingk = (xis-xi)/xi/xis;
   double myz = zofc(xi);
   double scalef = 1./(1.+myz);
   double prefactor = pow3(3.*params[0]*pow2(h0)/(2.*scalef) * lensingk);


   double k1 = l1/xi;
   double k2 = l2/xi;
   double k3 = sqrt(pow2(k1) + pow2(k2) + 2.*k2*k1*x);


    D_spt =  Dlin(myz) * dnorm_spt;
    F_spt =  F2DGP(myz) * pow2(dnorm_spt);
    vars[0] = mys8spl(myz);
    vars[1] = myknlspl(myz);
    vars[2] = 1.;
    vars[3] = (1.- F_spt * 7./2./pow2(D_spt));

     vars[4] = myneffspl(k1);
     vars[5] = myneffspl(k2);
     F2A = F2fit(vars,k1,k2,x);

     vars[4] = myneffspl(k1);
     vars[5] = myneffspl(k3);
     F2B = F2fit(vars,k1,k3,-(k2*x+k1)/k3);

     vars[4] = myneffspl(k2);
     vars[5] = myneffspl(k3);
     F2C = F2fit(vars,k2,k3,-(k1*x+k2)/k3);

    phpars[0] = phpars0(myz);
    phpars[1] = phpars1(myz);
    phpars[2] = phpars2(myz);
    phpars[3] = phpars3(myz);
    phpars[4] = phpars4(myz);
    phpars[5] = phpars5(myz);
    phpars[6] = phpars6(myz);
    phpars[7] = phpars7(myz);
    phpars[8] = phpars8(myz);
    phpars[9] = phpars9(myz);
    phpars[10] = phpars10(myz);
    phpars[11] = phpars11(myz);
    phpars[12] = phpars12(myz);

     if (k1>100.) {
      p1 = spt.PHALO(100.);
     }
     else{
     p1 = spt.PHALO(k1);
        }
        if (k2>100.) {
          p2 = spt.PHALO(100.);
        }
        else{
        p2 = spt.PHALO(k2);
           }

           if (k3>100.) {
              p3 = spt.PHALO(100.);
           }
           else{
           p3 = spt.PHALO(k3);
              }

    double bis =2.*(p1*p2*F2A + p3*p1*F2B + p2*p3*F2C);

    return prefactor * bis * pow2(xi);
}


// SC fitting formula
static double bsclssa_integrand(const Cosmology& C,const PowerSpectrum& P_L,  double params[], double l1, double l2, double x, double xi ){
   SPT spt(C, P_L, 1e-3);
   double vars[6],F2A,F2B,F2C,p1,p2,p3;

   double h0 = 1./2997.92;
   double lensingk = (xis-xi)/xi/xis;
   double myz = zofc(xi);
   double scalef = 1./(1.+myz);
   double prefactor = pow3(3.*params[0]*pow2(h0)/(2.*scalef) * lensingk);


   double k1 = l1/xi;
   double k2 = l2/xi;
   double k3 = sqrt(pow2(k1) + pow2(k2) + 2.*k2*k1*x);

    D_spt =  Dlin(myz) * dnorm_spt;
    F_spt =  F2DGP(myz) * pow2(dnorm_spt);
    vars[0] = mys8spl(myz);
    vars[1] = myknlspl(myz);
    vars[2] = 1.;
    vars[3] = (1.- F_spt * 7./2./pow2(D_spt));

     vars[4] = myneffspl(k1);
     vars[5] = myneffspl(k2);
     F2A = F2fitsc(vars,k1,k2,x);

     vars[4] = myneffspl(k1);
     vars[5] = myneffspl(k3);
     F2B = F2fitsc(vars,k1,k3,-(k2*x+k1)/k3);

     vars[4] = myneffspl(k2);
     vars[5] = myneffspl(k3);
     F2C = F2fitsc(vars,k2,k3,-(k1*x+k2)/k3);

    phpars[0] = phpars0(myz);
    phpars[1] = phpars1(myz);
    phpars[2] = phpars2(myz);
    phpars[3] = phpars3(myz);
    phpars[4] = phpars4(myz);
    phpars[5] = phpars5(myz);
    phpars[6] = phpars6(myz);
    phpars[7] = phpars7(myz);
    phpars[8] = phpars8(myz);
    phpars[9] = phpars9(myz);
    phpars[10] = phpars10(myz);
    phpars[11] = phpars11(myz);
    phpars[12] = phpars12(myz);

     if (k1>100.) {
      p1 = spt.PHALO(100.);
     }
     else{
     p1 = spt.PHALO(k1);
        }
        if (k2>100.) {
          p2 = spt.PHALO(100.);
        }
        else{
        p2 = spt.PHALO(k2);
           }

           if (k3>100.) {
              p3 = spt.PHALO(100.);
           }
           else{
           p3 = spt.PHALO(k3);
              }

    double bis =2.*(p1*p2*F2A + p3*p1*F2B + p2*p3*F2C);

    return prefactor * bis * pow2(xi);
}

// Numerical tree level
static double btnlssa_integrand(const Cosmology& C, const PowerSpectrum& P_L, int a, double params[], double l1, double l2, double x, double xi ){
   BSPT bspt(C, P_L, 1e-3);
   IOW iow;
   double bis, kernel, F2A, F2B, F2C, growth;
   double h0 = 1./2997.92;
   double lensingk = (xis-xi)/xi/xis;
   double myz = zofc(xi);
   double scalef = 1./(1.+myz);
   double prefactor = pow3(3.*params[0]*pow2(h0)/(2.*scalef) * lensingk);
   double k1 = l1/xi;
   double k2 = l2/xi;
   double k3 = sqrt(pow2(k1) + pow2(k2) + 2.*k2*k1*x);

   double vars[6];
   vars[0] = params[0];
   vars[1] = params[1];
   vars[2] = params[2];
   vars[3] = params[3];
   vars[4] = scalef;
   vars[5] = 1.*a; // doubling as choice of dgp (1), ide (2) , fr(3)

   bis = bspt.Btreen(vars, k1, k2, x);
   return prefactor * bis * pow2(xi);
}


// numerical 1-loop level
static double blnlssa_integrand(const Cosmology& C, const PowerSpectrum& P_L, int a, double params[], double l1, double l2, double x, double xi ){
   BSPT bspt(C, P_L, 1e-4);
   IOW iow;
   double bis, kernel, F2A, F2B, F2C, growth;
   double h0 = 1./2997.92;
   double lensingk = (xis-xi)/xi/xis;
   double myz = zofc(xi);
   double scalef = 1./(1.+myz);
   double prefactor = pow3(3.*params[0]*pow2(h0)/(2.*scalef) * lensingk);
   double k1 = l1/xi;
   double k2 = l2/xi;
   double vars[9];
   vars[0] = params[0];
   vars[1] = params[1];
   vars[2] = params[2];
   vars[3] = params[3];
   vars[4] = scalef;
   vars[5] = 1.*a; // doubling as choice of dgp (1), ide (2) , fr(3)
   vars[6] = params[4];
   vars[7] = params[5];
   vars[8] = params[6];
   bis = bspt.Bloopn(vars, k1, k2, x);
   return prefactor * bis * pow2(xi);
}


// POST BORN STUFF

// C_l(xi1, xi2) integrands
// linear
static double cl_integrandt(const Cosmology& C, const PowerSpectrum& P_L, double omega0, double l2, double xi1, double xi2){
  double h0 = 1./2997.92;
  double k2 = l2/xi2;
  double lenker1 = (xi1 - xi2)/xi2/xi1 ;
  double lenker2 = (xis - xi2)/xi2/xis;
  double weylpot2 =  pow2(Dlin(zofc(xi2))) * P_L(k2) * pow2(3. * omega0 * pow2(h0) * (1.+zofc(xi2)) / 2. / pow2(k2));

  return  lenker1 * lenker2 / pow2(xi2) * weylpot2 ;
}

// 1-loop (regpt)
static double cl_integrandl(const Cosmology& C, const PowerSpectrum& P_L, double omega0, double l2, double xi1, double xi2){
  double D0, D0sqr, D2, expon, p13, pl,  ploop,ploopr,h0,k2,weylpot2,lenker1,lenker2;
   D0 =  pow2(Dlin(zofc(xi2)));
   h0 = 1./2997.92;
   k2 = l2/xi2;

  expon = D0 * damp_termcmb(k2) ;
  if(expon >= 80.){
   return 0.;
  }
  else{

   lenker1 = (xi1 - xi2)/xi2/xi1;
   lenker2 = (xis - xi2)/xi2/xis;

  D0sqr = pow2(D0);
  pl = D0 * P_L(k2);
  p13 =  D0sqr*p13spl(k2);
  ploop = pl+ D0sqr*p22spl(k2) + p13;
  D2 = pow2(p13/pl/2.)*D0;

  ploopr = exp(-expon)*(ploop + expon/2.*(2.*pl + p13) + pow2(expon)/4.*pl + P_L(k2)*D2);

  weylpot2 =   ploopr * pow2(3. * omega0 * pow2(h0) * (1.+zofc(xi2)) / 2. / pow2(k2));

  return  lenker1 * lenker2 / pow2(xi2) * weylpot2 ;
}
}


// 1-loop (spt)
static double cl_integrandspt(const Cosmology& C, const PowerSpectrum& P_L, double omega0, double l2, double xi1, double xi2){
  double D0, D0sqr, D2, expon, p13, pl,  ploop,ploopr,h0,k2,weylpot2,lenker1,lenker2;
   D0 =  pow2(Dlin(zofc(xi2)));
   h0 = 1./2997.92;
   k2 = l2/xi2;
   lenker1 = (xi1 - xi2)/xi2/xi1;
   lenker2 = (xis - xi2)/xi2/xis;

  D0sqr = pow2(D0);
  pl = D0 * P_L(k2);
  p13 =  D0sqr*p13spl(k2);
  ploop = pl+ D0sqr*p22spl(k2) + p13;

  weylpot2 =   ploop * pow2(3. * omega0 * pow2(h0) * (1.+zofc(xi2)) / 2. / pow2(k2));

  return  lenker1 * lenker2 / pow2(xi2) * weylpot2 ;
}



// halofit
static double cl_integrandh(const Cosmology& C, const PowerSpectrum& P_L, double omega0, double l2, double xi1, double xi2){
  SPT spt(C, P_L, 1e-3);

  double h0 = 1./2997.92;
  double k2 = l2/xi2;
  double lenker1 = (xi1 - xi2)/xi2/xi1 ;
  double lenker2 = (xis - xi2)/xi2/xis;
  double myz = zofc(xi2);

    D_spt =  Dlin(myz) * dnorm_spt;
    phpars[0] = phpars0(myz);
    phpars[1] = phpars1(myz);
    phpars[2] = phpars2(myz);
    phpars[3] = phpars3(myz);
    phpars[4] = phpars4(myz);
    phpars[5] = phpars5(myz);
    phpars[6] = phpars6(myz);
    phpars[7] = phpars7(myz);
    phpars[8] = phpars8(myz);
    phpars[9] = phpars9(myz);
    phpars[10] = phpars10(myz);
    phpars[11] = phpars11(myz);
    phpars[12] = phpars12(myz);

  double weylpot2 =   spt.PHALO(k2) * pow2(3. * omega0 * pow2(h0) * (1.+zofc(xi2)) / 2. / pow2(k2));

  return  lenker1 * lenker2 / pow2(xi2) * weylpot2 ;
}


// M(l1,l2) function integrands
// tree level
static double mfunc_integrandt(const Cosmology& C, const PowerSpectrum& P_L, double omega0, double l1, double l2, double xi1){
  double h0 = 1./2997.92;
  double mycl = pow4(l2)*Integrate(bind(cl_integrandt,cref(C),cref(P_L), omega0, l2, xi1, _1), xim,xi1,1e-4);  // C_l2 (xi1, xis)
  double k1 = l1/xi1;
  double weylpot1 =  pow2(Dlin(zofc(xi1))) * P_L(k1) * pow2(3. * omega0 * pow2(h0) * (1.+zofc(xi1)) / 2. / pow2(k1));

  return pow2((xis - xi1)/pow2(xi1)/xis) * weylpot1 * mycl;
}

// 1-loop level (RegPT implementation)
static double mfunc_integrandl(const Cosmology& C, const PowerSpectrum& P_L, double omega0, double l1, double l2, double xi1){
  double D0,D0sqr, D2, expon, p13, pl,  ploop,ploopr,h0,mycl,k1,weylpot1;
  h0 = 1./2997.92;
  D0 =  pow2(Dlin(zofc(xi1)));
  k1 = l1/xi1;
  expon = D0 * damp_termcmb(k1) ;
  if(expon >= 80. ){
   return 0.;
  }
  else{
     mycl = pow4(l2) * Integrate(bind(cl_integrandl, cref(C), cref(P_L), omega0, l2, xi1, _1), xim,xi1,1e-2);   // C_l2 (xi1, xis)
     D0sqr = pow2(D0);
     pl = D0 * P_L(k1);
     p13 =  D0sqr*p13spl(k1);
     ploop = pl+ D0sqr*p22spl(k1) + p13;
     D2 = pow2(p13/pl/2.)*D0;

     ploopr = exp(-expon)*(ploop + expon/2.*(2.*pl + p13) + pow2(expon)/4.*pl + P_L(k1)*D2);

     weylpot1 =  ploopr * pow2(3. * omega0 * pow2(h0) * (1.+zofc(xi1)) / 2. / pow2(k1));

  return pow2( (xis - xi1)/pow2(xi1)/xis ) * weylpot1 * mycl;
}
}

// pure spt
static double mfunc_integrandspt(const Cosmology& C, const PowerSpectrum& P_L, double omega0, double l1, double l2, double xi1){
  double D0,D0sqr, D2, expon, p13, pl,  ploop,ploopr,h0,mycl,k1,weylpot1;
  h0 = 1./2997.92;
  D0 =  pow2(Dlin(zofc(xi1)));
  k1 = l1/xi1;

     mycl = pow4(l2) * Integrate(bind(cl_integrandspt, cref(C), cref(P_L), omega0, l2, xi1, _1), xim,xi1,1e-2);   // C_l2 (xi1, xis)
     D0sqr = pow2(D0);
     pl = D0 * P_L(k1);
     p13 =  D0sqr*p13spl(k1);
     ploop = pl+ D0sqr*p22spl(k1) + p13;

     weylpot1 =  ploop * pow2(3. * omega0 * pow2(h0) * (1.+zofc(xi1)) / 2. / pow2(k1));

  return pow2( (xis - xi1)/pow2(xi1)/xis ) * weylpot1 * mycl;
}



// halofit
static double mfunc_integrandh(const Cosmology& C, const PowerSpectrum& P_L, double omega0, double l1, double l2, double xi1){
  SPT spt(C, P_L, 1e-2);

  double h0 = 1./2997.92;
  double mycl = pow4(l2)*Integrate(bind(cl_integrandh,cref(C), cref(P_L), omega0, l2, xi1, _1), xim,xi1,1e-4);   // C_l2 (xi1, xis)
  double k1 = l1/xi1;
  double myz = zofc(xi1);

    D_spt =  Dlin(myz) * dnorm_spt;
    phpars[0] = phpars0(myz);
    phpars[1] = phpars1(myz);
    phpars[2] = phpars2(myz);
    phpars[3] = phpars3(myz);
    phpars[4] = phpars4(myz);
    phpars[5] = phpars5(myz);
    phpars[6] = phpars6(myz);
    phpars[7] = phpars7(myz);
    phpars[8] = phpars8(myz);
    phpars[9] = phpars9(myz);
    phpars[10] = phpars10(myz);
    phpars[11] = phpars11(myz);
    phpars[12] = phpars12(myz);

  double weylpot1 =  spt.PHALO(k1) * pow2(3. * omega0 * pow2(h0) * (1.+zofc(xi1)) / 2. / pow2(k1));

  return pow2( (xis - xi1)/pow2(xi1)/xis ) * weylpot1 * mycl;
}

// M(l1,l2) as in Eq.4.5 of 1605.05662
static double Mfunc(int a, const Cosmology& C, const PowerSpectrum& P_L, double omega0, double l1, double l2){
    switch (a) {
      case 1:
      return pow4(l1) * Integrate(bind(mfunc_integrandt,cref(C),cref(P_L), omega0, l1, l2, _1), xim,xis40,1e-3); // tree level
      break;
      case 2:
      return pow4(l1) * Integrate(bind(mfunc_integrandl,cref(C),cref(P_L), omega0, l1, l2, _1), xim,xis40,1e-3); // 1-loop (regpt)
      break;
      case 3:
      return pow4(l1) * Integrate(bind(mfunc_integrandh,cref(C),cref(P_L), omega0, l1, l2, _1), xim,xis,1e-3); // halofit
      break;
      case 4:
      return pow4(l1) * Integrate(bind(mfunc_integrandspt,cref(C),cref(P_L), omega0, l1, l2, _1), xim,xis40,1e-3); // 1-loop (spt)
      break;
}}

// B_LSS contribution
// a chooses LCDM (1) or nDGP/MG analytic (2) (tree or 1-loop only) or 3(IDE) or 4 [f(R)] (only numerical)
// b chooses tree analytic (1) or 1-loop analytic (2) or Gil-Marin (3) or Scoccimaro (4) or tree numerical (5) or 1-loop numerical (6)
// params[0] = omega0, params[1]=mg1, params[2]=mg2, params[3] =mg3

// LSS bispectrum
double CMBc::BLSS(int a, int b,  double params[], double l1, double l2, double x) const{
switch (b) {
  case 1:
   return   Integrate(bind(btlssa_integrand,cref(C),cref(P_L), a, params, l1, l2, x, _1),xim,xis,epsrel); // tree level EdS
   break;
  case 2:
    return  Integrate(bind(bllssa_integrand,cref(C),cref(P_L),a, params, l1, l2, x, _1),xim,xis40,epsrel); // 1-loop EdS
    break;
  case 3:
    return  Integrate(bind(bgmlssa_integrand,cref(C),cref(P_L), params, l1, l2, x, _1),xim,xis,epsrel); // Gil-Marin fit
    break;
  case 4:
    return  Integrate(bind(bsclssa_integrand,cref(C),cref(P_L), params, l1, l2, x, _1),xim,xis40,epsrel); // Scoccimaro fit
    break;
  case 5:
  if (a==1) {
    return  Integrate(bind(btlssa_integrand,cref(C),cref(P_L), a, params, l1, l2, x, _1),xim,xis40,epsrel); // For GR we use EdS
    }
  else{
    return  Integrate(bind(btnlssa_integrand,cref(C),cref(P_L),a-1, params, l1, l2, x, _1),xim,xis40,epsrel); // tree level numerical
  }
  break;
  case 6:
  if (a==1) {
    return  Integrate(bind(bllssa_integrand,cref(C),cref(P_L), a, params, l1, l2, x, _1),xim,xis40,epsrel);  // For GR we use EdS
    }
  else{
    return  Integrate(bind(blnlssa_integrand,cref(C),cref(P_L),a-1, params, l1, l2, x, _1),xim,xis40,epsrel); // 1-loop level numerical
  }
    break;
  }
}

// BLSS lensing kernel
double CMBc::BLSS_kernel(int a, int b, double params[], double l1, double l2, double x, double xi) const{
  switch (b) {
    case 1:
     return   btlssa_integrand(cref(C),cref(P_L), a, params, l1, l2, x, xi); // tree level EdS
     break;
    case 2:
      return   bllssa_integrand(cref(C),cref(P_L), a, params, l1, l2, x, xi); // 1-loop EdS
      break;
    case 3:
      return   bgmlssa_integrand(cref(C),cref(P_L), params, l1, l2, x, xi); // Gil-Marin fit
      break;
    }
}


// post born correction
double CMBc::BPB(int a, double omega0, double l1, double l2, double x) const{
   double l1sqr = pow2(l1);
   double l2sqr = pow2(l2);
   double l3 = sqrt(l1sqr + l2sqr + 2.*l1*l2*x);
   double l3sqr = pow2(l3);
   double l1l3 = -(l1sqr + l1*l2*x);
   double l2l3 = -(l2sqr + l1*l2*x);
   return   2.*x/l1/l2 * (l1l3 * Mfunc(a, cref(C), cref(P_L), omega0, l1,l2) + l2l3 * Mfunc(a,cref(C), cref(P_L), omega0, l2,l1))
           +2.*l2l3/l2sqr/l3sqr * (l2*l1*x * Mfunc(a,cref(C), cref(P_L), omega0, l2,l3) + l1l3 * Mfunc(a,cref(C), cref(P_L), omega0, l3,l2))
           +2.*l1l3/l3sqr/l1sqr * (l2l3 * Mfunc(a,cref(C), cref(P_L), omega0, l3,l1) + l2*l1*x * Mfunc(a,cref(C), cref(P_L), omega0, l1,l3));
    }


// CMB lensing spectrum
// nDGP or LCDM are chosen via comdist_init(int a, double omega0, double mg1) function (Only linear growth dependence - 1-loop uses UsA + EdS)
double CMBc::CMBla(int a, int b, double params[],double l1, double l2, double x) const{
  double l1sqr = pow2(l1);
  double l2sqr = pow2(l2);
  double l3 = sqrt(l1sqr + l2sqr + 2.*l1*l2*x);
  int lw1,lw2,lw3;
  lw1 = (int)l1;
  lw2 = (int)l2;
  lw3 = (int)round(l3);
  return  wigner3j(lw1, lw2, lw3) * sqrt((2.*l1+1.)*(2.*l2+1.)*(2.*l3+1.)/4./M_PI) * (BLSS(a,b, params, l1,l2,x) + BPB(a, params[0], l1,l2,x));
}
