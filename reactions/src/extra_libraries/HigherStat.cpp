#if HAVE_CONFIG_H
# include <config.h>
#endif

#include "HigherStat.h"
#include "PowerSpectrum.h"
#include "LinearPS.h"
#include "Quadrature.h"
#include "SpecialFunctions.h"
#include "array.h"
#include "Spline.h"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_deriv.h>

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cerrno>
#include <cstdlib>
#include <stdlib.h>
#include <sys/stat.h>
#include <math.h>
#include <vector>
#include <boost/bind.hpp>
using boost::cref;

const real KMIN_e = 0.0001;
const real KMAX_e = 30.;

// Load Linear PS
HigherStat::HigherStat(const Cosmology& C, const PowerSpectrum& P_l, real epsrel_)
: C(C), P_l(P_l),kintegral(*this)
{
MonteCarloIntegral::Range domain[5] = {
    { KMIN_e, KMAX_e },
    { KMIN_e, KMAX_e },
    { KMIN_e, KMAX_e },
    { -0.999999, 0.999999 },
    { -0.999999, 0.999999 },
};
kintegral.SetDomain(domain);
epsrel = epsrel_;
}
HigherStat::~HigherStat() {
}


//Filters
inline double w_filter(double k, int a){
  switch (a) {
  case 1:
    return 3./k*SphericalBesselJ1(k); //TOPHAT
    break;
  case 2:
	  return exp(-pow2(k)/2.); //GAUSSIAN
    break;
  case 3:
	  return 1.; // NO FILTER
    break;
    default:
     warning("Filter: invalid indices");
     return 0;
    }
}


/*Gamma 1 and Gamma 2 (logarithmic derivatives of variance) */

// Create splines functions of sigma^2
Spline sig2_spline;
void HigherStat::sigsqr_init(int a) const{
  int loop_N = 1000;
  double RMAX =250;
  double RMIN = 1.;
  double Rval, sigval, vars[4];
  vector<double> Rval_tab, sigval_tab;
for(int i = 0; i< loop_N; i++){
    Rval = RMIN + i*(RMAX-RMIN)/(loop_N-1);//RMIN * exp(i*log(RMAX/RMIN)/(loop_N-1));
    sigval = sig_sqr(Rval,a,1,vars);
    Rval_tab.push_back(Rval);
    sigval_tab.push_back(sigval);
      }
  sig2_spline = LinearSpline(Rval_tab,sigval_tab);
}

//function for derivative

static double mysigsqr(double R,void * params){
  (void)(params);
  return sig2_spline(R);
}

static double gam1(double R, void * params1){
  (void)(params1);
gsl_function F;
double result, abserr;
  F.params = 0;
  F.function = &mysigsqr;
  gsl_deriv_central(&F, R, 1e-3, &result, &abserr);
return R*result/sig2_spline(R);
}

static double gam2(double R){
gsl_function F;
double result, abserr;
  F.params = 0;
  F.function = &gam1;
  gsl_deriv_central(&F, R, 1e-4, &result, &abserr);
return R*result;
}



// b selects analytic (1) or numerical (2)
// a selects window function form

/* Variance */
static double var_integranda(const PowerSpectrum& P_l, double R, int a,double k){
	return pow2(k*w_filter(k*R,a)*D_spt/dnorm_spt)*P_l(k);
}

static double gam1_test(double R, const PowerSpectrum& P_l){
double step = 0.00001;
double top = Integrate(bind(var_integranda, cref(P_l), R + step, 1, _1), KMIN_e , KMAX_e, 1e-4)/(2.*pow2(M_PI));
double bot = Integrate(bind(var_integranda, cref(P_l), R - step, 1, _1), KMIN_e , KMAX_e, 1e-4)/(2.*pow2(M_PI));
double mid = Integrate(bind(var_integranda, cref(P_l), R, 1, _1), KMIN_e , KMAX_e, 1e-4)/(2.*pow2(M_PI));
return R*(top-bot)/(2.*step)/mid;
}

static double gam2_test(double R, double step, const PowerSpectrum& P_l){
double top = gam1_test(R+step, cref(P_l));
double bot = gam1_test(R-step, cref(P_l));
double mid = gam1_test(R, cref(P_l));
return R*(top-bot)/(2.*step)/mid;
}




static double var_integrandn(const PowerSpectrum& P_l, double R, int a, double vars[], double k){
  IOW iow;
  iow.initn_lin(0, vars[4], k, vars[0], vars[1], vars[2], vars[3]);
	return pow2(k*w_filter(k*R,a)*F1_nk/dnorm_spt)*P_l(k);
}

// vars is for numerical skewness only : vars[0]=omega0, vars[1-3] = MG params and vars[4]=scale factor
// a chooses filter, b chooses analytic or numerical (1 or 2 resp)
double HigherStat::sig_sqr(double R, int a, int b, double vars[]) const{
  switch (b) {
    case 1:
        return  Integrate(bind(var_integranda, cref(P_l), R, a, _1), KMIN_e , KMAX_e, epsrel)/(2.*pow2(M_PI));
    break;
    case 2:
         return  Integrate(bind(var_integrandn, cref(P_l), R, a, vars, _1), KMIN_e , KMAX_e, epsrel)/(2.*pow2(M_PI));
    break;
  }

}

/* Skewness */
// analytic (GR AND nDGP)
static double s3_num_integranda(const PowerSpectrum& P_l, double R, int a, double k1, double k2, double x){
   double k1s = pow2(k1);
   double k2s = pow2(k2);
	return k1s*k2s*w_filter(k1*R,a)*w_filter(k2*R,a)*w_filter(sqrt(k1s+k2s+2.*k1*k2*x)*R,a)*(pow2(D_spt)*F2eds(k1,k2,x) + F_spt*(1.-x*x))*pow2(D_spt)/pow4(dnorm_spt)*P_l(k1)*P_l(k2);
}
// numerical
static double s3_num_integrandn(const PowerSpectrum& P_l, double R, int a, double vars[], double k1, double k2, double x){
  double kargs[4],kv[3],xv[3], p22,p13,d;
  double k1s = pow2(k1);
  double k2s = pow2(k2);
  IOW iow;
        kv[0] = k1;
        kv[1] = k2;
        kv[2] = kv[1];
        xv[0] = -0.99999999;
        xv[1] = x;
        xv[2] = -x;
        kargs[0] = sqrt(kv[1]*kv[1]+kv[0]*kv[0]-2.*kv[1]*kv[0]*xv[1]);
        kargs[2] = sqrt(kv[1]*kv[1]+2.*kv[1]*kv[2]*xv[0]+kv[2]*kv[2]);
        kargs[1] = sqrt(kv[2]*kv[2]+2.*kv[2]*kv[0]*xv[2]+kv[0]*kv[0]);
        kargs[3] = sqrt(kv[0]*kv[0]+2.*kv[0]*kv[1]*xv[1]+kv[1]*kv[1]);
        iow.initn2(vars[4],kv,xv,kargs,vars[0],vars[1],vars[2],vars[3]);
	return k1s*k2s*w_filter(k1*R,a)*w_filter(k2*R,a)*w_filter(sqrt(k1s+k2s+2.*k1*k2*x)*R,a)*F2A_nk[0]*F1p_nk[0]*F1_nk/pow4(dnorm_spt)*P_l(k1)*P_l(k2);
}

// vars is for numerical skewness only : vars[0]=omega0, vars[1-3] = MG params and vars[4]=scale factor
// a chooses window functions (see w_filter)
// b chooses (1) analytic (2) numerical (3) gamma1 approx for top-hat filter
double HigherStat::skewness(double R, int a, int b, double vars[]) const{
  double c[3] = {KMIN_e,KMIN_e,-0.99999999};
  double d[3] = {KMAX_e,KMAX_e,0.99999999};
  void * params;

  switch (b) {
    case 1:
	     return  3./4./pow4(M_PI) * Integrate<3>(bind(s3_num_integranda, cref(P_l), R, a, _1, _2, _3), c,d, epsrel);
    break;
    case 2:
      return  3./4./pow4(M_PI) * Integrate<3>(bind(s3_num_integrandn, cref(P_l), R, a, vars, _1, _2, _3), c,d, epsrel);
    break;
    case 3:
      return (34./7. + gam1(R,params)) * pow2(sig2_spline(R));
    break;
  }
}

/* Analytic Kurtosis */

static double s4_num_integranda1(const PowerSpectrum& P_l, double R, int a, double args[], double k[], double x12, double x23, double x13){
   double k1 = k[0];
   double k2 = k[1];
   double k3 = k[2];
   double k1s = args[0];
   double k2s = args[1];
   double k3s = args[2];
   double filt1 = args[3];
   double filt2 = args[4];
   double filt3 = args[5];
   double k12 = sqrt(k1s+k2s+2.*k1*k2*x12);
   double k23 = sqrt(k2s+k3s-2.*k2*k3*x23);
   double k13 = sqrt(k1s+k3s+2.*k1*k3*x13);
   double filt12  = w_filter(k12*R,a);
   double filt23 = w_filter(k23*R,a);
   double term1, term2;

   term1 = F2eds(k1,k2,x12)*F2eds(k2,k3,-x23)*filt23*filt12*filt1*filt3;
   term2 = 0.;
   if (a>=2) {
     double k123 = sqrt(k1s + k2s + k3s + 2.*k1*k2*x12 + 2.*k2*k3*x23 + 2.*k1*k3*x13);
     double filt123 = w_filter(k123*R,a);
     term2 = 0.5*filt1*filt2*filt3*F3edsb(k1,k2,k3,k23,k12,k13,x23,x12,x13)*filt123;
   }
	return term1 + term2;
}

static double s4_num_integranda2(const PowerSpectrum& P_l, double R, int a, double k1, double k2, double k3){
   double c[3] = {-1.,-1.,-1.};
   double d[3] = {1.,1.,1.};
   double k[3];
   k[0] = k1;
   k[1] = k2;
   k[2] = k3;
   double args[6];
   args[0]= pow2(k1);
   args[1]= pow2(k2);
   args[2]= pow2(k3);
   args[3]= w_filter(k1*R,a);
   args[4]= w_filter(k2*R,a);
   args[5]= w_filter(k3*R,a);
	return args[0]*args[1]*args[2]*P_l(k1)*P_l(k2)*P_l(k3)*Integrate<3>(bind(s4_num_integranda1, cref(P_l), R, a,args,k, _1, _2, _3), c,d, 1e-4);
}

// vars is for numerical kurtosis only : vars[0]=omega0, vars[1-3] = MG params and vars[4]=scale factor
// a chooses window functions (see w_filter)
// b chooses (1) analytic (2) numerical (3) gamma1 + gamma2 approx for top-hat filter


HigherStat::kurtosis_integral::kurtosis_integral(const HigherStat& hs_)
    : MonteCarloIntegral(5), hs(hs_)
{
    maxeval = 10000000;
}

void HigherStat::kurtosis_integral::Integrand(const double x[], double* f, double* param) const {
  double R = *param;
  double k1 = x[0];
  double k2 = x[1];
  double k3 = x[2];
  double x12 = x[3];
  double x23 = x[4];
  double k1s= pow2(k1);
  double k2s = pow2(k2);
  double k3s =  pow2(k3);
  double filt1 = w_filter(k1*R,1);
  double filt3 = w_filter(k3*R,1);
  double k12 = sqrt(k1s+k2s+2.*k1*k2*x12);
  double k23 = sqrt(k2s+k3s+2.*k2*k3*x23);
  double filt12  = w_filter(k12*R,1);
  double filt23  = w_filter(k23*R,1);
  double term1 =  F2eds(k1,k2,x12)*F2eds(k2,k3,x23)*filt23*filt12*filt1;
  f[0] = 2.*k1s*k2s*k3s*hs.P_l(k1)*hs.P_l(k2)*hs.P_l(k3)*filt3*(term1);
}


/* kbar definition (EQ. A.25 OF 9312026)*/
static double kbar_integrand(const PowerSpectrum& P_l, double R, double k){
	return pow2(pow2(k)*w_filter(k*R,1)*D_spt/dnorm_spt)*P_l(k);
}

double HigherStat::kurtosis(double R, int a, int b, double vars[]) const{
  double c[3] = {KMIN_e,KMIN_e,KMIN_e};
  double d[3] = {KMAX_e,KMAX_e,KMAX_e};
  double term4 = 0.;
  double kbar;
  double sigsqr = sig_sqr(R, 1, 1, vars);
  if (b==1 && a ==1) {
    kbar = Integrate(bind(kbar_integrand, cref(P_l), R, _1), KMIN_e , KMAX_e, epsrel)/(2.*pow2(M_PI)*sigsqr);
    term4 =  2./189.*pow3(sigsqr)*(1364. + 9.*gam1_test(R,cref(P_l)) * (60. + 7. * gam1_test(R,cref(P_l))) - 126. * kbar);
  }
  switch (b) {
    case 1:
      return  term4;//0.75 * pow6(D_spt/(M_PI*dnorm_spt)) * Integrate<3>(bind(s4_num_integranda2, cref(P_l), R, a, _1, _2, _3), c,d, epsrel); //+ term4;
    case 2:
      return (60712./1323. + 62./3.*gam1_test(R,cref(P_l)) + 7./3.*pow2(gam1_test(R,cref(P_l))) + 2./3.*gam2_test(R,epsrel,cref(P_l))) * pow3(sigsqr);
    case 3:
    return 0.75 * pow6(D_spt/(M_PI*dnorm_spt)) * Integrate<3>(bind(s4_num_integranda2, cref(P_l), R, a, _1, _2, _3), c,d, epsrel);
  }
}
