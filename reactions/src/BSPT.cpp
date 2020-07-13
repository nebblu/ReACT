
#if HAVE_CONFIG_H
# include <config.h>
#endif

#include "Cosmology.h"
#include "PowerSpectrum.h"
#include "Quadrature.h"
#include "BSPT.h"
#include "BSPTN.h"
#include "SPT.h"
#include "Spline.h"
#include "SpecialFunctions.h"
#include "LinearPS.h"


#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include <math.h>       /* pow */

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_deriv.h>

#include <random>
#include <functional>

using std::cref;
using std::bind;

BSPT::BSPT(const Cosmology& C_, const PowerSpectrum& P_L_, real epsrel_)
: C(C_), P_L(P_L_)
{
    epsrel = epsrel_;
}


static double neff_integrand(const PowerSpectrum& P_L, double lambda, double k, double q){
return P_L(q)/neffehu(q)*exp(-pow2(log(k/q)/lambda)/2.)/q;
}


static double myfunc(double k,void * params){
  (void)(params);
  return log(neffnw(k));
}

static double myneff(double k){
gsl_function F;
double result, abserr;
F.function = &myfunc;
F.params = 0;
abserr = 1e-5;
gsl_deriv_central(&F, k, 1e-5, &result, &abserr);
return k*result;
}

// Calculate sigma_8 at some redshift
static double sigma8int2(const PowerSpectrum& P_L, double k){
  return pow2(3.*((sin(k*8.) - k*8.*cos(k*8.))/pow2(k*8.))/k/8.) * pow2(k) * pow2(D_spt/dnorm_spt) * P_L(k)/2./M_PI/M_PI;
}

static double mysigma82(const PowerSpectrum& P_L){
  return sqrt(Integrate(bind(sigma8int2,cref(P_L),std::placeholders::_1),0.0001,30.,1e-4));
}

// Find K_nl
  double myknls2(const PowerSpectrum& P_L)
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

// No wiggle spectrum for neff used in Vlah et al
Spline neffnw;
Spline neffehu;
Spline mylinearps;
double mys8,myknl;

// Initialize the NW power spectrum, A^{nw,0l}, sigma8 and knl
// j initialises k_nl and sigma8 at the redshift set by inite (defined in GR)
void BSPT::mypspline(double scalef, double omega0, int j) const{
SPT spt(C, P_L, 1e-2);
NoWigglePS now(C, 0. , EisensteinHu);
vector<double> kval_table, now_table, ehu_table, linps_table;
double k, nowps, lambda, ehups, linps;
const double kmin = 1e-4;
const double kmax = 80.;
const int n1 = 2000;
// initialise halofit functions
spt.phinit(scalef, omega0);

for(int i = 0; i<n1; i++ ){
  k = kmin*exp(i*log(kmax/kmin)/(n1-1.));

  kval_table.push_back(k);
  ehups = now.Evaluate(k);
  linps = P_L(k);

  ehu_table.push_back(ehups);
  linps_table.push_back(linps);
    }
  neffehu = CubicSpline(kval_table,ehu_table);
  mylinearps = CubicSpline(kval_table,linps_table);

  for(int i = 0; i<n1; i++ ){
    k = kmin*exp(i*log(kmax/kmin)/(n1-1.));
    lambda = 0.25*pow(k/0.05,0.04);
    nowps = now.Evaluate(k)/sqrt(2.*M_PI*pow2(lambda))*Integrate<ExpSub>(bind(neff_integrand,cref(P_L),lambda, k,std::placeholders::_1),kmin,kmax,1e-5);
    now_table.push_back(nowps);
      }
    neffnw = CubicSpline(kval_table,now_table);
    if(j==1){
    mys8 = mysigma82(cref(P_L));
    myknl = myknls2(cref(P_L));
     }

  }

// Fitting formula prescription of 1805.10567

real BSPT::Bfit(double k1, double k2, double x) const {
  double k3 =  sqrt(pow2(k1) + pow2(k2) + 2.*k2*k1*x);
  double F2A,F2B,F2C, vars[6];
    vars[0] = mys8;
    vars[1] = myknl;
    vars[2] = 1.;
    vars[3] = (1.- F_spt * 7./2./pow2(D_spt));
    vars[4] = myneff(k1);
    vars[5] = myneff(k2);
    F2A = F2fit(vars,k1,k2,x);
    vars[4] = myneff(k1);
    vars[5] = myneff(k3);
    F2B = F2fit(vars,k1,k3,-(k2*x+k1)/k3);
    vars[4] = myneff(k2);
    vars[5] = myneff(k3);
    F2C = F2fit(vars,k2,k3,-(k1*x+k2)/k3);

    SPT spt(C, P_L, 1e-3);
    double p1,p2,p3;
    p1 = spt.PHALO(k1);
    p2 = spt.PHALO(k2);
    p3 = spt.PHALO(k3);

  return 2.*(p1*p2*F2A + p3*p1*F2B + p2*p3*F2C);
}



// GR/nDGP Tree level Bispectrum (Real space)
// a==1 (GR analytic depending if you run inite )
// a=2 (DGP analytic)
real BSPT::Btree(int a, double k1, double k2, double x) const {
  double k3 =  sqrt(pow2(k1) + pow2(k2) + 2.*k2*k1*x);
  double F2A,F2B,F2C,growth;
  switch (a) {
    case 1:
     growth = pow2(Dl_spt/dnorm_spt);
     F2A = growth*F2eds(k1, k2, x);
     F2B = growth*F2eds(k1,k3,-(k2*x+k1)/k3);
     F2C = growth*F2eds(k2,k3,-(k1*x+k2)/k3);
     break;
   case 2:
    growth = pow2(D_spt/dnorm_spt);
    F2A = F2ndgp(k1,k2,x);
    F2B = F2ndgp(k1,k3,-(k2*x+k1)/k3);
    F2C = F2ndgp(k2,k3,-(k1*x+k2)/k3);
    break;
   default:
    warning("BSPT: invalid indices, a = %d\n", a);
    return 0;
  }
//Tree
 return     2.*growth*(P_L(k1)*P_L(k2)*F2A + P_L(k3)*P_L(k1)*F2B + P_L(k2)*P_L(k3)*F2C);
}


// GR/DGP 1-loop Bispectrum (Real space)

// a=1 (GR analytic)
// a=2 (DGP analytic)
real BSPT::Bloop(int a, double k1, double k2, double x) const {
  SPT spt(C, P_L, 1e-2);
  double k3 =  sqrt(pow2(k1) + pow2(k2) + 2.*k2*k1*x);
  double F2A,F2B,F2C,P131,P132,P133,growth;
  switch (a) {
  case 1:
    growth = pow2(Dl_spt/dnorm_spt);
    F2A = growth*F2eds(k1, k2, x);
    F2B = growth*F2eds(k1,k3,-(k2*x+k1)/k3);
    F2C = growth*F2eds(k2,k3,-(k1*x+k2)/k3);
    P131 = spt.P13_dd(k1);
    P132 = spt.P13_dd(k2);
    P133 = spt.P13_dd(k3);
    break;
  case 2:
   growth = pow2(D_spt/dnorm_spt);
   F2A = F2ndgp(k1, k2, x);
   F2B = F2ndgp(k1,k3,-(k2*x+k1)/k3);
   F2C = F2ndgp(k2,k3,-(k1*x+k2)/k3);
   P131 = spt.P13D_dd(k1);
   P132 = spt.P13D_dd(k2);
   P133 = spt.P13D_dd(k3);
   break;
   default:
         warning("BSPT: invalid indices, a = %d\n", a);
         return 0;
}
//Tree
 return     2.*growth*(P_L(k1)*P_L(k2)*F2A + P_L(k3)*P_L(k1)*F2B + P_L(k2)*P_L(k3)*F2C)
// B222, B321-I and B411
            + Bloopterms(a,k1,k2,k3,x)
// B123-II
            + F2A*(P_L(k1)*P132 + P_L(k2)*P131) + F2B*(P_L(k3)*P131 + P_L(k1)*P133) + F2C*(P_L(k3)*P132 + P_L(k2)*P133);
 }

 static real B_222(const PowerSpectrum& P_L, real k1, real k2, real x, real k3, real r, real u, real o) {
   double k2s = pow2(k2);

   double p = k1*r;
   double ps = pow2(p);
   double d = sqrt(1. + r*r - 2.*r*u);

   real k1mp = k1*d;
   real k2p = u*x - sqrt((1.-pow2(u))*(1.-pow2(x)))*cos(o);
   real k2pp = sqrt(k2s + ps + 2*k2*p*k2p);

     if(d < 1e-5 )
         return 0;
     else{
 return   8.*pow2(r)*P_L(p)*(F2eds(r*k1,k1mp,(u-r)/d)*F2eds(p,k2pp,-(k2*k2p+p)/k2pp)*F2eds(k1mp,k2pp,(k2*x+p*u-k2*r*k2p-p*r)/d/k2pp)*P_L(k1mp)*P_L(k2pp));
 }
 }

 static real B_222D(const PowerSpectrum& P_L, real k1, real k2, real x, real k3, real r, real u, real o) {
   double k2s = pow2(k2);

   double p = k1*r;
   double ps = pow2(p);
   double d = sqrt(1. + r*r - 2.*r*u);

   real k1mp = k1*d;
   real k2p = u*x - sqrt((1.-pow2(u))*(1.-pow2(x)))*cos(o);
   real k2pp = sqrt(k2s + ps + 2*k2*p*k2p);

     if(d < 1e-5 )
         return 0;
     else{
 return   8.*pow2(r)*P_L(p)*(F2ndgp(r*k1,k1mp,(u-r)/d)*F2ndgp(p,k2pp,-(k2*k2p+p)/k2pp)*F2ndgp(k1mp,k2pp,(k2*x+p*u-k2*r*k2p-p*r)/d/k2pp)*P_L(k1mp)*P_L(k2pp));
 }
 }


 static real B_321I(const PowerSpectrum& P_L, real k1, real k2, real x, real k3, real r, real u, real o) {
     double k1s = pow2(k1);
     double k2s = pow2(k2);
     double p = k1*r;
     double ps = pow2(p);
     double d = sqrt(1. + r*r - 2.*r*u);
     double d1 = sqrt(1. + r*r + 2.*r*u);

     real k1mp = k1*d;
     real k2p = u*x - sqrt((1.-pow2(u))*(1.-pow2(x)))*cos(o);
     real k2mp = sqrt(k2s + ps - 2.*k2*p*k2p);
     real k3mp = sqrt(k2s + k1s + ps + 2.*k1*k2*x + 2.*k1s*r*u + 2.*k2*p*k2p);

     double kernel1= F2eds(r*k1,k1mp,(u-r)/d);
     double kernel2= F2eds(r*k1,k2mp,(k2*k2p-p)/k2mp);
     double kernel3 = F2eds(r*k1,k3mp,-(k1*u+k2*k2p+p)/k3mp);

     if(d < 1e-5 || d1<1e-5)
         return 0;
     else{
 return    6.*pow2(r)*P_L(p)*(
           kernel2*F3eds(k1,p,k2mp,(k2*k2p-p)/k2mp,u,(k2*x-p*u)/k2mp)*P_L(k2mp)*P_L(k1)
        +  kernel3*F3eds(k1,p,k3mp,-(k1*u+k2*k2p+p)/k3mp,u,-(k1+k2*x+p*u)/k3mp)*P_L(k3mp)*P_L(k1)

        +  kernel1*F3eds(k2,p,k1mp,(u-r)/d,k2p,(x-r*k2p)/d)*P_L(k1mp)*P_L(k2)
        +  kernel3*F3eds(k2,p,k3mp,-(k1*u+k2*k2p+p)/k3mp,k2p,-(k2+k1*x+p*k2p)/k3mp)*P_L(k3mp)*P_L(k2)

       +  kernel1*F3eds(k3,p,k1mp,(u-r)/d,-(k1*u+k2*k2p)/k3, (-k1+p*u-k2*x+k2*r*k2p)/d/k3)*P_L(k1mp)*P_L(k3)
        +  kernel2*F3eds(k3,p,k2mp,(k2*k2p-p)/k2mp,-(k1*u+k2*k2p)/k3,-(k1*k2*x-k1*p*u+k2s- k2*p*k2p)/k3/k2mp)*P_L(k2mp)*P_L(k3));

 }}


  static real B_321ID(const PowerSpectrum& P_L, real k1, real k2, real x, real k3, real r, real u, real o) {
      double k1s = pow2(k1);
      double k2s = pow2(k2);
      double p = k1*r;
      double ps = pow2(p);
      double d = sqrt(1. + r*r - 2.*r*u);
      double d1 = sqrt(1. + r*r + 2.*r*u);

      real k1mp = k1*d;
      real k2p = u*x - sqrt((1.-pow2(u))*(1.-pow2(x)))*cos(o);
      real k2mp = sqrt(k2s + ps - 2.*k2*p*k2p);
      real k3mp = sqrt(k2s + k1s + ps + 2.*k1*k2*x + 2.*k1s*r*u + 2.*k2*p*k2p);

      double kernel1= F2ndgp(p,k1mp,(u-r)/d);
      double kernel2= F2ndgp(p,k2mp,(k2*k2p-p)/k2mp);
      double kernel3 = F2ndgp(p,k3mp,-(k1*u+k2*k2p+p)/k3mp);

      if(d < 1e-5 || d1<1e-5)
          return 0;
      else{
  return    D_spt/dnorm_spt*6.*pow2(r)*P_L(p)*(
            kernel2*F3ndgp(k1,p,k2mp,(k2*k2p-p)/k2mp,u,(k2*x-p*u)/k2mp)*P_L(k2mp)*P_L(k1)
        +  kernel3*F3ndgp(k1,p,k3mp,-(k1*u+k2*k2p+p)/k3mp,u,-(k1+k2*x+p*u)/k3mp)*P_L(k3mp)*P_L(k1)

         +  kernel1*F3ndgp(k2,p,k1mp,(u-r)/d,k2p,(x-r*k2p)/d)*P_L(k1mp)*P_L(k2)
         +  kernel3*F3ndgp(k2,p,k3mp,-(k1*u+k2*k2p+p)/k3mp,k2p,-(k2+k1*x+p*k2p)/k3mp)*P_L(k3mp)*P_L(k2)

        +  kernel1*F3ndgp(k3,p,k1mp,(u-r)/d,-(k1*u+k2*k2p)/k3, (-k1+p*u-k2*x+k2*r*k2p)/d/k3)*P_L(k1mp)*P_L(k3)
        +  kernel2*F3ndgp(k3,p,k2mp,(k2*k2p-p)/k2mp,-(k1*u+k2*k2p)/k3,-(k1*k2*x-k1*p*u+k2s- k2*p*k2p)/k3/k2mp)*P_L(k2mp)*P_L(k3));

  }}



 static real B_411(const PowerSpectrum& P_L, real k1, real k2, real x, real k3, real r, real u, real o) {
     double p = k1*r;
     double k2p = u*x - sqrt((1.-pow2(u))*(1.-pow2(x)))*cos(o);
     return   12.*pow2(r)*P_L(p)*(F4edsb(k1,k2,p,p,-k2p,k2p,x,-u,u)*P_L(k1)*P_L(k2)
                               +  F4edsb(k2,k3,p,p,(k1*u+k2*k2p)/k3,-(k1*u+k2*k2p)/k3,-(k1*x+k2)/k3,-k2p,k2p)*P_L(k2)*P_L(k3)
                               +  F4edsb(k3,k1,p,p,-u,u,-(k1+k2*x)/k3,(k1*u+k2*k2p)/k3,-(k1*u+k2*k2p)/k3)*P_L(k3)*P_L(k1));
 }

 static real B_411D(const PowerSpectrum& P_L, real k1, real k2, real x, real k3, real r, real u, real o) {
     double p = k1*r;
     double k2p = u*x - sqrt((1.-pow2(u))*(1.-pow2(x)))*cos(o);
     double karg1[4],karg2[4],karg3[4];
     double xarg1[6],xarg2[6],xarg3[6];
     karg1[0]=k1;
     karg1[1]=k2;
     karg1[2]=p;
     karg1[3]=p;

     karg2[0]=k2;
     karg2[1]=k3;
     karg2[2]=p;
     karg2[3]=p;

     karg3[0]=k1;
     karg3[1]=k3;
     karg3[2]=p;
     karg3[3]=p;

     xarg1[0]=k2p;
     xarg1[1]=x;
     xarg1[2]=-u;
     xarg1[3]=-k2p;
     xarg1[4]= XMIN;
     xarg1[5]=u;

     xarg2[0]=-(k1*u+k2*k2p)/k3;
     xarg2[1]=-(k1*x+k2)/k3;
     xarg2[2]=-k2p;
     xarg2[3]=(k1*u+k2*k2p)/k3;
     xarg2[4]= XMIN;
     xarg2[5]=k2p;

     xarg3[0]=u;
     xarg3[1]=-(k1+k2*x)/k3;
     xarg3[2]=(k1*u+k2*k2p)/k3;
     xarg3[3]=-u;
     xarg3[4]= XMIN;
     xarg3[5]=-(k1*u+k2*k2p)/k3;

     return  pow2(D_spt/dnorm_spt)*12.*pow2(r)*P_L(p)*(F4dgp(karg1,xarg1)*P_L(k1)*P_L(k2)
                                        +F4dgp(karg2,xarg2)*P_L(k2)*P_L(k3)
                                        +F4dgp(karg3,xarg3)*P_L(k1)*P_L(k3));
 }


double BSPT::Bloopterms(int a, real k1, real k2, real k3, real x) const{
    	real KMAX = QMAXp/k1;
   	  real KMIN = QMINp/k1;
      real tp = 2.*M_PI;
      real c[3] = {KMIN,-1.,0.};
      real d[3] = {KMAX,0.99999999,tp};
      switch (a) {
        case 1:
        return pow6(Dl_spt/dnorm_spt)*pow3(k1/tp)*(Integrate<3>(bind(B_222, cref(P_L), k1, k2, x, k3, std::placeholders::_1,std::placeholders::_2,std::placeholders::_3), c, d,epsrel)
                                       +    Integrate<3>(bind(B_321I, cref(P_L), k1, k2, x, k3, std::placeholders::_1,std::placeholders::_2,std::placeholders::_3), c, d,epsrel)
                                       +    Integrate<3>(bind(B_411, cref(P_L), k1, k2, x, k3, std::placeholders::_1,std::placeholders::_2,std::placeholders::_3), c, d,epsrel));
        break;
        case 2:
        return   pow3(k1/tp)*(Integrate<3>(bind(B_222D, cref(P_L), k1, k2, x, k3, std::placeholders::_1,std::placeholders::_2,std::placeholders::_3), c, d,epsrel)
                          +    Integrate<3>(bind(B_321ID, cref(P_L), k1, k2, x, k3, std::placeholders::_1,std::placeholders::_2,std::placeholders::_3), c, d,epsrel)
                          +   Integrate<3>(bind(B_411D, cref(P_L), k1, k2, x, k3, std::placeholders::_1,std::placeholders::_2,std::placeholders::_3), c, d,epsrel));
          break;
        default:
            warning("BSPT: invalid indices, a = %d \n", a);
            return 0;
      }

}


///////////////////// NUMERICAL TREE LEVEL SPECTRUM ////////////////////////
// 0 = omega0, 1 = mg1, 2=mg2, 3=mg3, 4=a, 5=mg o
// vars[5] = 1 (DGP or GR)
// vars[5] = 2 (f(R) Hu-Sawicki)
real BSPT::Btreen(double vars[], double k1, double k2, double x) const {
BSPTN bsptn;
double ka[3],xa[3],kargs[3];
double k1s,k2s,k3;
k1s = pow2(k1);
k2s = pow2(k2);
ka[0]=k1;
ka[1]=k2;
ka[2]=sqrt(k1s+k2s+2.*x*k1*k2);
k3 = ka[2];
kargs[0] = k3;
kargs[1] = k2;
kargs[2] = k1;
xa[0] = -(k1*x+k2)/k3; // k2.k3
xa[1] = x; //k1.k2
xa[2] = -(k2*x+k1)/k3;// k1.k3


switch( (int)vars[5] ) {
  case 1:
  bsptn.initnb0_dgp(vars[4], ka, xa, kargs, vars[0], vars[1], vars[2],vars[3]);
  return 2./pow4(dnorm_spt)*(P_L(k1)*P_L(k2)*F1b_k1*F1b_k2*F2b_k12 + P_L(k3)*P_L(k1)*F1b_k1*F1b_k3*F2b_k13 + P_L(k2)*P_L(k3)*F1b_k2*F1b_k3*F2b_k23);

  case 2:
  bsptn.initnb0_fofr(vars[4], ka, xa, kargs, vars[0], vars[1], vars[2],vars[3]);
  return 2./pow4(dnorm_spt)*(P_L(k1)*P_L(k2)*F1b_k1*F1b_k2*F2b_k12 + P_L(k3)*P_L(k1)*F1b_k1*F1b_k3*F2b_k13 + P_L(k2)*P_L(k3)*F1b_k2*F1b_k3*F2b_k23);

  default:
   warning("BSPT: invalid indices, vars[5] = %e\n", vars[5]);
  return 0.;

}
}
///////////////////// NUMERICAL 1-LOOP SPECTRUM ///////////////////////////

// MG INITIALISATION
 static real bloopdgp(const PowerSpectrum& P_L, double vars[], double k1, double k2, double x, double k3, double p, double mu, double phi){
   BSPTN bsptn;
   double B411,B222,B321a,B321b,pkp,pk1,pk2,pk3,pka1,pka4,pka3;
   double scalef = vars[4];
   double omega0= vars[0];
   double mg1= vars[1];
   double mg2= vars[2];
   double mg3= vars[3];
   double epars[3];
   epars[0]= vars[7];
   epars[1]= vars[8];
   epars[2]= 1.; //vars[9];

    double k2p,k2pp,k3p,ps,ka[8],xa[8],kargs[25];
    double sx = (1.-pow2(x));
    double k1s = pow2(k1);
    double k2s = pow2(k2);
    double k3s = pow2(k3);
    kargs[5] = x;
    ka[5] = k1;
    ka[6] = k2;
    ka[7] = k3;
    kargs[9] = -(k1 + k2*x)/k3; // k1.k3
    kargs[13] = -(k2 + k1*x)/k3; // k2.k3
   ka[0] = p;
   ps = pow2(p);
   k2pp = sqrt(sx*(1.-pow2(mu)));
   ka[1] = sqrt(k1s+ps-2.*p*k1*mu); // k1-p
   xa[1] = (k1*mu-p)/ka[1]; // p. k1-p
   xa[5] = mu;
   kargs[1] = sqrt(k1s + ps + 2.*p*k1*mu);  // p + k1
   kargs[14] = -(k1*mu + p)/kargs[1]; // - p . k1+p
   k2p = mu*x - k2pp*cos(phi*M_PI+M_PI); // k2.p
   k3p = (-k1*mu-k2*k2p)/k3; // k3.p
   ka[2] = sqrt(k2s+ps+2.*p*k2*k2p); // k2+p
   ka[3] = sqrt(k2s + ps - 2.*k2*p*k2p); // k2 - p
   ka[4] = sqrt(k3s + ps - 2.*k3*p*k3p); // k3-p
   xa[0] = (k1*k2*x + k1*p*mu - p*k2*k2p - ps)/ka[1]/ka[2];  //k1-p . k2+p
   xa[2] = -(k2*k2p+p)/ka[2]; //k2+p . -p
   xa[3] = (k2*k2p - p)/ka[3];   // k2-p . p
   xa[4] =  (k3*k3p - p)/ka[4];  // k3-p.p
   xa[6] = k2p;
   xa[7] = k3p;
   kargs[0] = sqrt(k2s + 4.*ps + 4.*k2*p*k2p); // k2 + 2p
   kargs[2] = sqrt(k3s + ps + 2.*k3*p*k3p); // k3+p
   kargs[3] = (k2*x - p*mu)/ka[3]; // k2-p . k1 = k2x - pu /k2mp
   kargs[4] = (k2*p*k2p + k1*k2*x - ps -p*k1*mu)/kargs[1]/ka[3]; // k2-p . p+k1
   kargs[6] = (k1*mu + k2*k2p - p)/kargs[2];  // p. -k3-p
   kargs[7] = -(k1+k2*x+p*mu)/ka[4]; // k3-p. k1
   kargs[8] = (k3*k1* kargs[9] + k3p*k3*p-p*k1*mu-ps)/ka[4]/kargs[1]; //k3-p.k1+p
   kargs[10] = (k1*x+p*k2p)/kargs[1]; //  k2.k1+p
   kargs[11] = -(k2+k1*x+p*k2p)/ka[4]; // k3-p.k2
   kargs[12] = (kargs[13]*k2*k3 -k2p*k2*p+ k3p*k3*p - ps)/ka[4]/ka[2]; // k2+p.k3-p
   kargs[15] = (k1*x - p*k2p)/ka[1]; // k1-p.k2
   kargs[16] = (k1*x + k2 - p*k2p)/kargs[2]; // k2. -k3-p
   kargs[17] = (kargs[13]*k2*k3 - p*k3*k3p)/k3/ka[3]; // k2-p.k3
   kargs[18] = (k2p*k2*p + kargs[13]*k3*k2 - ps -k3*p*k3p)/ka[3]/kargs[2]; // k2-p.p+k3
   kargs[19] = -(k1*k2*x + k1s - k2*p*k2p - k1*p*mu)/ka[1]/k3; // k1-p . k3
   kargs[20] =  (k1*p*mu + k3*k1*kargs[9] - ps - p*k3*k3p)/kargs[2]/ka[1];// k1-p.p+k3
   kargs[21] = (k2*x + p*mu)/ka[2]; // k1.k2+p
   kargs[22] = (-k1s-k1*p*mu-k1*k2*x-k2*k2p*p)/kargs[1]/k3; // k1+p.k3
   kargs[23] = (-k1*k2*x - k2s + k3p*p*k3)/k3/ka[2]; // k2+p.k3
   kargs[24] = (-k1-k2*x+p*mu)/kargs[2];//k1.k3+p
   pkp = P_L(p);
   pk1 = P_L(k1);
   pk2= P_L(k2);
   pk3 = P_L(k3);
   pka1 = P_L(ka[1]);
   pka3 = P_L(ka[3]);
   pka4 = P_L(ka[4]);

     bsptn.initnb1_dgp(scalef, ka, xa, kargs, omega0, mg1, mg2, mg3,epars);

     B222  = 8.*ps*pkp* F2b_k13 * F2b_k12 * F2b_k23 * P_L(ka[2])* pka1;
     B321a = 6.*ps*pkp*( F2b_12a * (pk1 * pk2 * F3b_2pp * F1b_1 + pk2 * pk1 * F3b_1pp * F1b_2)
 											            +F2b_13a * (pk3 * pk1 * F3b_1pp * F1b_3 + pk1 * pk3 * F3b_3pp * F1b_1)
 												          +F2b_23a * (pk3 * pk2 * F3b_2pp * F1b_3 + pk2 * pk3 * F3b_3pp * F1b_2));
     B321b = 6.*ps*pkp*(pk1 * F1b_1 *(F2b_p2mp * F3b_12mp * pka3  + F2b_p3mp * F3b_13mp * pka4)
 											          +  pk2 * F1b_2 *(F2b_k12 * F3b_21mp * pka1  +  F2b_p3mp * F3b_23mp * pka4)
 												        +  pk3 * F1b_3 *(F2b_k12 * F3b_31mp * pka1  +  F2b_p2mp * F3b_32mp * pka3));
     B411 = 12.*ps*pkp*(F4b_12pp * F1b_1 * F1b_2 * pk1 * pk2
                                + F4b_13pp * F1b_1 * F1b_3 * pk1 * pk3
                                + F4b_23pp * F1b_2 * F1b_3 * pk2 * pk3);

    return B222+B321a+B321b+B411;

      }


      static real bloopfr(const PowerSpectrum& P_L, double vars[], double k1, double k2, double x, double k3, double p, double mu, double phi){
        BSPTN bsptn;
        double B411,B222,B321a,B321b,pkp,pk1,pk2,pk3,pka1,pka4,pka3;

        double scalef = vars[4];
        double omega0= vars[0];
        double mg1= vars[1];
        double mg2= vars[2];
        double mg3= vars[3];
        double epars[3];
        epars[0]=vars[7];
        epars[1]=vars[8];
        epars[2]=1.; //vars[9];

         double k2p,k2pp,k3p,ps,ka[8],xa[8],kargs[36];
         double sx = (1.-pow2(x));
         double k1s = pow2(k1);
         double k2s = pow2(k2);
         double k3s = pow2(k3);
         kargs[5] = x;
         ka[5] = k1;
         ka[6] = k2;
         ka[7] = k3;
         kargs[9] = -(k1 + k2*x)/k3; // k1.k3
         kargs[13] = -(k2 + k1*x)/k3; // k2.k3
        ka[0] = p;
        ps = pow2(p);
        k2pp = sqrt(sx*(1.-pow2(mu)));
        ka[1] = sqrt(k1s+ps-2.*p*k1*mu); // k1-p
        xa[1] = (k1*mu-p)/ka[1]; // p. k1-p
        xa[5] = mu;
        kargs[1] = sqrt(k1s + ps + 2.*p*k1*mu);  // p + k1
        kargs[14] = -(k1*mu + p)/kargs[1]; // - p . k1+p
        k2p = mu*x - k2pp*cos(phi*M_PI+M_PI); // k2.p
        k3p = (-k1*mu-k2*k2p)/k3; // k3.p
        ka[2] = sqrt(k2s+ps+2.*p*k2*k2p); // k2+p
        ka[3] = sqrt(k2s + ps - 2.*k2*p*k2p); // k2 - p
        ka[4] = sqrt(k3s + ps - 2.*k3*p*k3p); // k3-p
        xa[0] = (k1*k2*x + k1*p*mu - p*k2*k2p - ps)/ka[1]/ka[2];  //k1-p . k2+p
        xa[2] = -(k2*k2p+p)/ka[2]; //k2+p . -p
        xa[3] = (k2*k2p - p)/ka[3];   // k2-p . p
        xa[4] =  (k3*k3p - p)/ka[4];  // k3-p.p
        xa[6] = k2p;
        xa[7] = k3p;
        kargs[0] = sqrt(k2s + 4.*ps + 4.*k2*p*k2p); // k2 + 2p
        kargs[2] = sqrt(k3s + ps + 2.*k3*p*k3p); // k3+p
        kargs[3] = (k2*x - p*mu)/ka[3]; // k2-p . k1 = k2x - pu /k2mp
        kargs[4] = (k2*p*k2p + k1*k2*x - ps -p*k1*mu)/kargs[1]/ka[3]; // k2-p . p+k1
        kargs[6] = (k1*mu + k2*k2p - p)/kargs[2];  // p. -k3-p
        kargs[7] = -(k1+k2*x+p*mu)/ka[4]; // k3-p. k1
        kargs[8] = (k3*k1* kargs[9] + k3p*k3*p-p*k1*mu-ps)/ka[4]/kargs[1]; //k3-p.k1+p
        kargs[10] = (k1*x+p*k2p)/kargs[1]; //  k2.k1+p
        kargs[11] = -(k2+k1*x+p*k2p)/ka[4]; // k3-p.k2
        kargs[12] = (kargs[13]*k2*k3 -k2p*k2*p+ k3p*k3*p - ps)/ka[4]/ka[2]; // k2+p.k3-p
        kargs[15] = (k1*x - p*k2p)/ka[1]; // k1-p.k2
        kargs[16] = (k1*x + k2 - p*k2p)/kargs[2]; // k2. -k3-p
        kargs[17] = (kargs[13]*k2*k3 - p*k3*k3p)/k3/ka[3]; // k2-p.k3
        kargs[18] = (k2p*k2*p + kargs[13]*k3*k2 - ps -k3*p*k3p)/ka[3]/kargs[2]; // k2-p.p+k3
        kargs[19] = -(k1*k2*x + k1s - k2*p*k2p - k1*p*mu)/ka[1]/k3; // k1-p . k3
        kargs[20] =  (k1*p*mu + k3*k1*kargs[9] - ps - p*k3*k3p)/kargs[2]/ka[1];// k1-p.p+k3
        kargs[21] = (k2*x + p*mu)/ka[2]; // k1.k2+p
        kargs[22] = (-k1s-k1*p*mu-k1*k2*x-k2*k2p*p)/kargs[1]/k3; // k1+p.k3
        kargs[23] = (-k1*k2*x - k2s + k3p*p*k3)/k3/ka[2]; // k2+p.k3
        kargs[24] = (-k1-k2*x+p*mu)/kargs[2];//k1.k3+p

        kargs[25] = ps;
        kargs[26] = pow2(ka[1]);
        kargs[27] = pow2(ka[2]);
        kargs[28] = pow2(ka[3]);
        kargs[29] = pow2(ka[4]);
        kargs[30] = k1s;
        kargs[31] = k2s;
        kargs[32] = k3s;
        kargs[33] = pow2(kargs[0]);
        kargs[34] = pow2(kargs[1]);
        kargs[35] = pow2(kargs[2]);

        pkp = P_L(p);
        pk1 = P_L(k1);
        pk2= P_L(k2);
        pk3 = P_L(k3);
        pka1 = P_L(ka[1]);
        pka3 = P_L(ka[3]);
        pka4 = P_L(ka[4]);

          bsptn.initnb1_fr(scalef, ka, xa, kargs, omega0, mg1, mg2, mg3,epars);

          B222  = 8.*ps*pkp*F2b_k13 * F2b_k12 * F2b_k23 * P_L(ka[2])* pka1;

          B321a = 6.*ps*pkp*( F2b_12a * (pk1 * pk2 * F3b_2pp * F1b_1 + pk2 * pk1 * F3b_1pp * F1b_2)
      											            +F2b_13a * (pk3 * pk1 * F3b_1pp * F1b_3 + pk1 * pk3 * F3b_3pp * F1b_1)
      												          +F2b_23a * (pk3 * pk2 * F3b_2pp * F1b_3 + pk2 * pk3 * F3b_3pp * F1b_2));
          B321b = 6.*ps*pkp*(pk1 * F1b_1 *(F2b_p2mp * F3b_12mp * pka3  + F2b_p3mp * F3b_13mp * pka4)
      											          +  pk2 * F1b_2 *(F2b_k12 * F3b_21mp * pka1  +  F2b_p3mp * F3b_23mp * pka4)
      												        +  pk3 * F1b_3 *(F2b_k12 * F3b_31mp * pka1  +  F2b_p2mp * F3b_32mp * pka3));
          B411 = 12.*ps*pkp*(F4b_12pp * F1b_1 * F1b_2 * pk1 * pk2
                                     + F4b_13pp * F1b_1 * F1b_3 * pk1 * pk3
                                     + F4b_23pp * F1b_2 * F1b_3 * pk2 * pk3);

         return  B321a + B321b + B222 + B411 ;

           }

// vars :
// 0 = omega0, 1 = mg1, 2=mg2, 3=mg3, 4=a, 5=which model , 6 = integration accuracy
double BSPT::Bloopn(double vars[], double k1, double k2, double x)const{
      BSPTN bsptn;
      int check = vars[5];
      double loop,tree,prefac;
      double c[3] = {QMINp,-0.99999999,-1.};
      double d[3] = {QMAXp,0.9999999,1.};
      double k3 = sqrt(pow2(k1)+pow2(k2)+2.*k1*k2*x);
      double ka[3],kargs[3],xa[3];
      ka[0]=k1;
      ka[1]=k2;
      ka[2]=k3;
      k3 = ka[2];
      kargs[0] = k3;
      kargs[1] = k2;
      kargs[2] = k1;
      xa[0] = -(k1*x+k2)/k3; // k2.k3
      xa[1] = x; //k1.k2
      xa[2] = -(k2*x+k1)/k3;// k1.k3
      switch (check) {
         case 1:
          prefac = M_PI/pow3(pow2(dnorm_spt)*2.*M_PI);
          if (k1<0.005 && k2 <0.005 && k3<0.005) {
            bsptn.initnb0_dgp(vars[4], ka, xa, kargs, vars[0], vars[1], vars[2],vars[3]);
            tree = 2./pow4(dnorm_spt)*(P_L(k1)*P_L(k2)*F1b_k1*F1b_k2*F2b_k12 + P_L(k3)*P_L(k1)*F1b_k1*F1b_k3*F2b_k13 + P_L(k2)*P_L(k3)*F1b_k2*F1b_k3*F2b_k23);
            loop = 0.;
          }
          else{
          loop = prefac*Integrate<3>(bind(bloopdgp,cref(P_L),vars,k1,k2,x,k3,std::placeholders::_1,std::placeholders::_2,std::placeholders::_3), c, d, vars[6]);
          tree = 2./pow4(dnorm_spt)*(P_L(k1)*P_L(k2)*F1b_1*F1b_2*F2b_12a + P_L(k3)*P_L(k1)*F1b_1*F1b_3*F2b_13a + P_L(k2)*P_L(k3)*F1b_2*F1b_3*F2b_23a);
          }
          return tree+loop;

        case 2:
          prefac = M_PI/pow3(pow2(dnorm_spt)*2.*M_PI);
          if (k1<0.005 && k2 <0.005 && k3<0.005) {
            bsptn.initnb0_fofr(vars[4], ka, xa, kargs, vars[0], vars[1], vars[2],vars[3]);
            loop = 0.;
            tree = 2./pow4(dnorm_spt)*(P_L(k1)*P_L(k2)*F1b_k1*F1b_k2*F2b_k12 + P_L(k3)*P_L(k1)*F1b_k1*F1b_k3*F2b_k13 + P_L(k2)*P_L(k3)*F1b_k2*F1b_k3*F2b_k23);
          }
          else{
          loop = prefac*Integrate<3>(bind(bloopfr,cref(P_L),vars,k1,k2,x,k3,std::placeholders::_1,std::placeholders::_2,std::placeholders::_3), c, d, vars[6]);
          tree = 2./pow4(dnorm_spt)*(P_L(k1)*P_L(k2)*F1b_1*F1b_2*F2b_12a + P_L(k3)*P_L(k1)*F1b_1*F1b_3*F2b_13a + P_L(k2)*P_L(k3)*F1b_2*F1b_3*F2b_23a);
          }
          return  tree+loop;
          default:
              warning("BSPT: invalid indices, a = %d \n", check);
              return 0;
        }
  }
