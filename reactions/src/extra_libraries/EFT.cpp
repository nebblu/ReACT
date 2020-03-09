#if HAVE_CONFIG_H
# include <config.h>
#endif

#include "EFT.h"
#include "SPT.h"
#include "NoWigglePS.h"
#include "PowerSpectrum.h"
#include "LinearPS.h"
#include "Quadrature.h"
#include "SpecialFunctions.h"
#include "array.h"
#include "Spline.h"
#include <gsl/gsl_sf_bessel.h>

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


// Load Linear PS
EFT::EFT(const Cosmology& C, const PowerSpectrum& P_l, double epsrel_)
: C(C), P_l(P_l)
{
    epsrel = epsrel;
}

double error_eft = 1e-3;

//Atsushi limit on k-angular integration
inline double ATS(double k, double r){
 double KMAX = QMAXp/k;
 double KMIN = QMINp/k;
 if(r>=0.5) {
 return 1./(2.*r);
 }
 else{
 return  Min(0.999999999999, (1.+r*r-KMIN*KMIN)/2./r);
 }
}


// No wiggle spectrum
Spline RENW;
// Spline E-HU No-Wiggle PS
Spline EHUNW;
// linear growth splined
Spline ling_spline;
real sigma_v;
double ANW;


// K magnitude integrals for ANW
const real KMIN_e = 1e-5;
const real KMAX_e = 50.;


// function to set normalization of linear growth to fiducial for data comparisons
double sig_8_eft;
double rempf_eft;

// Normalise linear and logarithmic growth

void EFT::Geft(real dfid) const {
    sig_8_eft = dfid/D_spt;
         }

void EFT::rempeft(real ffid) const{
  rempf_eft = ffid/fl_spt;
}

static double ANW_integrand(double q, double k){
return ling_spline(k)*RENW(k)*(1.-j0(k*q))*pow2(q);
}

static double nowig_integrand(const PowerSpectrum& P_l, double k, double q){
double lambda = 0.25*pow(k/0.05,0.04);
return P_l(q)/EHUNW(q)*exp(-pow2(log(k/q)/lambda)/2.)/q;
}


static real eft_exp(const PowerSpectrum& P_l, double q){
  double growth = pow2(sig_8_eft*D_spt/dnorm_spt);
  return  P_l(q)*growth/(6.*pow2(M_PI));
}
// Initialization of sigma_v^2 k^2 term
// a = 0 :analytic
// a = 1 numeric


void EFT::sigmav_init(double scalef, double omega0, double mg1, double mg2, double mg3) const{
    sigma_v = Integrate(bind(eft_exp, cref(P_l), _1), KMIN_e, KMAX_e, error_eft);//LinearSpline(kval_table,sigd_table);
  }


// Initialize the NW power spectrum and A^{nw,0l}
void EFT::NWinit(double scalef, double omega0, double mg1, double mg2, double mg3) const{
NoWigglePS now(C, 0. , EisensteinHu);
IOW iow;
vector<double> kval_table, now_table, ehu_table, ling_table;
double k, nowps, ehups, ling;
const double qmin = 10.;
const double qmax = 300.;
real c[2] = {qmin,KMIN_e};
real d[2] = {qmax,KMAX_e};
int n1 = 2000;
ling = pow2(sig_8_eft*D_spt/dnorm_spt);



for(int i = 0; i<n1; i++ ){
  k = KMIN_e*exp(i*log(KMAX_e/KMIN_e)/(n1-1.));
  ehups = now.Evaluate(k);
  kval_table.push_back(k);
  ehu_table.push_back(ehups);
    }
  EHUNW = LinearSpline(kval_table,ehu_table);


  for(int i = 0; i<n1; i++ ){
    k = KMIN_e*exp(i*log(KMAX_e/KMIN_e)/(n1-1.));
    double lambda = 0.25*pow(k/0.05,0.04);
    nowps = EHUNW(k)/sqrt(2*M_PI*pow2(lambda))*Integrate<ExpSub>(bind(nowig_integrand,cref(P_l),k,_1),KMIN_e,KMAX_e,error_eft);
    now_table.push_back(nowps);
    ling_table.push_back(ling);

      }
    RENW = LinearSpline(kval_table,now_table);
    ling_spline = LinearSpline(kval_table,ling_table);
    ANW = 1./((pow3(qmax)-pow3(qmin))*pow2(M_PI))*Integrate<2>(bind(ANW_integrand,_1,_2),c,d,error_eft,error_eft);
  }


 // Global initialization - see below for all functions
  void EFT::eft_init(double s8, double scalef, double omega0, double mg1, double mg2, double mg3)const{
    IOW iow;
    SPT spt(C,P_l,1e-3);
    iow.inite(scalef, omega0 ,mg1, mg2, mg3); // evolution factors and normalization
    Geft(D_spt*s8);
    rempeft(fl_spt);
    spt.G(D_spt*s8);
    spt.remp(fl_spt);
    NWinit(scalef,omega0,mg1,mg2,mg3); // EHU NW initialisation
    sigmav_init(scalef,omega0,mg1,mg2,mg3); // sigma_{v,lin}
  }


// No Wiggle 1-loop : analytic
static double Df22_dd(double k, double r, double x) {
    double d = 1. + r*r - 2.*r*x;
    if(d < 1e-5)
        return 0;
    else
        return RENW(r*k) * RENW(k*sqrt(d))* 2 * r*r* (pow4(D_spt)*pow2((3*r + 7*x - 10*r*x*x)/(14.*r*d)) + pow2(F_spt) * pow2((x*x-1)/d) + pow2(D_spt)*F_spt*(3*r + 7*x - 10*r*x*x) * (1-x*x) / (7.*r*pow2(d)));}


static double Dmidintdd(double k, double r) {
    double KMAX = QMAXp/k;
    double KMIN = QMINp/k;
    double YMIN = Max(-1., (1.+r*r-KMAX*KMAX)/2./r);
    double YMAX = Min(0.99999999999, (1.+r*r-KMIN*KMIN)/2./r);
return Integrate(bind(Df22_dd, k, r, _1), YMIN, YMAX , error_eft);
}


/* P_{\delta\delta}^{(22)} */
double EFT::Pres22_dd(double k) const{
 	double KMAX = QMAXp/k;
	double KMIN = QMINp/k;
    return k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt)*Integrate<ExpSub>(bind(Dmidintdd, k, _1), KMIN, KMAX, error_eft);
}

static double Df13_dd(double k, double r) {
    double s;
	 if(r < 1e-2)
		s = pow4(D_spt)/252.*(-168 + (928./5.)*pow2(r) - (4512./35.)*pow4(r) + (416./21.)*pow6(r)) +  (D_spt*Cx_spt/6.) * (-96./5. * r*r + 96./35. * pow4(r) + 32./105.*pow6(r)) + (D_spt*I_spt/6.) * (32.*r*r - 96./5. * pow4(r) + 96./35. * pow6(r)) + (D_spt*KL_spt/12.) * ( 256./5. * r*r - 768./35. * pow4(r) + 256./105. * pow6(r)) ;

	else if(fabs(r-1) < 1e-10)

		s = pow4(D_spt)/252.*(-88 + 8*(r-1)) +  D_spt*Cx_spt/6. * (-16 - 24*(r-1)) + D_spt*I_spt/6. *(16 + 8*(r-1)) +  D_spt*KL_spt/12. *(32 + 32*(r-1)) ;

	else if(r>100)
		s = pow4(D_spt)/252.*(-488./5. + (96./5.)/pow2(r) - (160./21.)/pow4(r) - (1376./1155.)/pow6(r)) + D_spt*Cx_spt/6. * (-32 + 96./(5.*r*r) - 96./(35.*pow4(r)) - 32./(105.*pow6(r)))+ D_spt*I_spt/6. *(96./5. - 96./(35.*pow2(r)) - 32./(105.*pow4(r)) - 32./(385.*pow6(r))) + (D_spt*KL_spt/12.) * (256./5. - 768./(35.*pow2(r)) + 256./(105.*pow4(r)) + 256./(1155.*pow6(r))) ;

	else
   s = pow4(D_spt)/252.*(12./pow2(r) - 158. + 100.*pow2(r) - 42.*pow4(r) + 3./pow3(r) * pow3(r*r - 1) * (7*r*r + 2) * log((1+r)/fabs(1-r))) + D_spt*Cx_spt*(1./pow2(r) - 16./6. - pow2(r) + 1./(2.*pow3(r)) * pow3(r*r - 1) * log((1+r)/fabs(1-r))) + D_spt*I_spt*(1. + 16./6. * pow2(r) - pow4(r) + 1./(2.*r) * pow3(r*r - 1) * log((1+r)/fabs(1-r))) + D_spt*KL_spt/(12. * pow3(r)) * ( -6*r + 22*pow3(r) + 22*pow5(r) - 6*pow7(r) +  3*pow4(r*r - 1) * log((1+r)/fabs(1-r))) ;

    return RENW(k*r) * s;
}

/* P_{\delta\delta}^{(13)} */
double EFT::Pres13_dd(double k) const{
	double KMAX = QMAXp/k;
	double KMIN = QMINp/k;
    return  k*k*k/(4*M_PI*M_PI)*RENW(k)/pow4(dnorm_spt)*Integrate<ExpSub>(bind(Df13_dd, k, _1), KMIN, KMAX, error_eft,error_eft);
}


static real Df22_dt(real k, real r, real x) {
    real d = 1. + r*r - 2.*r*x;
    if(d < 1e-5)
        return 0;
    else
	 return  RENW(k*r) * RENW(k*sqrt(d)) * r * r * (rempf_eft*fdgp_spt*pow4(D_spt)*(3*r + 7*x - 10*r*x*x)*(r - 7*x + 6*x*x*r) / (98.*r*r*pow2(d)) - F_spt*Fd_spt/H_spt * 2 * pow2(x*x-1) / pow2(d) +  Fd_spt*pow2(D_spt)/H_spt * (3*r + 7*x - 10*r*x*x) * (x*x-1) / (7.*r*pow2(d)) - F_spt*rempf_eft*fdgp_spt*pow2(D_spt)*(r - 7*x + 6*r*x*x)*(x*x-1)/(7.*r*pow2(d)));
//  return RENW(r*k) * RENW(k*sqrt(d))* 2 * r*r* (pow4(D_spt)*pow2((3*r + 7*x - 10*r*x*x)/(14.*r*d)) + 0*pow2(F_spt) * pow2((x*x-1)/d) + 0*pow2(D_spt)*F_spt*(3*r + 7*x - 10*r*x*x) * (1-x*x) / (7.*r*pow2(d)));}
}


static real Dmidintdt(real k, real r) {
    real KMAX = QMAXp/k;
    real KMIN = QMINp/k;
    real YMIN = Max(-1., (1.+r*r-KMAX*KMAX)/2./r);
    real YMAX = Min(0.99999999999999, (1.+r*r-KMIN*KMIN)/2./r);
return Integrate(bind(Df22_dt, k, r, _1), YMIN, YMAX ,  error_eft);
}


/* P_{\delta\theta}^{(22)} */
real EFT::Pres22_dt(real k) const {
      real KMAX = QMAXp/k;
    	real KMIN = QMINp/k;
    return k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt) * Integrate<ExpSub>(bind(Dmidintdt, k, _1), KMIN, KMAX, error_eft);
}


static real Df13_dt( real k, real r) {
real s;

if(r < 1e-2)

s = -rempf_eft*fdgp_spt*pow4(D_spt)/252.*(-168 + (416./5.)*pow2(r) - (2976./35.)*pow4(r) + (224./15.)*pow6(r)) - (D_spt*Cdx_spt + Cx_spt*rempf_eft*fdgp_spt*D_spt*H_spt)/(12.*H_spt) * (-96./5. * r*r + 96./35. * pow4(r) + 32./105.*pow6(r)) - (rempf_eft*fdgp_spt*D_spt*H_spt*I_spt + Id_spt*D_spt -Fd_spt*pow2(D_spt))/(H_spt*12.)* (32*r*r - 96./5. * pow4(r) + 96./35. * pow6(r)) - (D_spt*KLd_spt + D_spt*rempf_eft*fdgp_spt*H_spt*KL_spt)/(H_spt*24.) * ( 256./5. * r*r - 768./35. * pow4(r) + 256./105. * pow6(r)) ;

else if(fabs(r-1) < 1e-10)

s = -rempf_eft*fdgp_spt*pow4(D_spt)/252.*(-152 - 56*(r-1)) - (D_spt*Cdx_spt + Cx_spt*rempf_eft*fdgp_spt*D_spt*H_spt)/(12*H_spt) * (-16 - 24*(r-1)) - (rempf_eft*fdgp_spt*D_spt*H_spt*I_spt + Id_spt*D_spt -Fd_spt*pow2(D_spt))/(H_spt*12.)*(16 + 8*(r-1)) - (D_spt*KLd_spt + D_spt*rempf_eft*fdgp_spt*H_spt*KL_spt)/(H_spt*24.) *(32 + 32*(r-1)) ;

else if(r>100)
s = -rempf_eft*fdgp_spt*pow4(D_spt)/252.*(-200 + (2208./35.)/pow2(r) - (1312./105.)/pow4(r) - (1888./1155.)/pow6(r)) -  (D_spt*Cdx_spt + Cx_spt*rempf_eft*fdgp_spt*D_spt*H_spt)/(12.*H_spt) *( -32 + 96./(5.*r*r) - 96./(35.*pow4(r)) - 32./(105.*pow6(r))) - (rempf_eft*fdgp_spt*D_spt*H_spt*I_spt + Id_spt*D_spt -Fd_spt*pow2(D_spt))/(H_spt*12.)* (96./5. - 96./(35.*pow2(r)) - 32./(105.*pow4(r)) - 32./(385.*pow6(r))) - (D_spt*KLd_spt + D_spt*rempf_eft*fdgp_spt*H_spt*KL_spt)/(12.*H_spt) * (256./5. - 768./(35.*pow2(r)) + 256./(105.*pow4(r)) + 256./(1155.*pow6(r))) ;

else
s = -rempf_eft*fdgp_spt*pow4(D_spt)/252.*(24/pow2(r) - 202 + 56*pow2(r) - 30*pow4(r) + 3/pow3(r) * pow3(r*r - 1) * (5*r*r + 4) * log((1+r)/fabs(1-r))) - (D_spt*Cdx_spt + Cx_spt*rempf_eft*fdgp_spt*D_spt*H_spt)/(H_spt*12.*pow3(r))*(6*r - 16*pow3(r) - 6*pow5(r) + 3*pow3(r*r-1)*log((1+r)/fabs(1-r))) - (rempf_eft*fdgp_spt*D_spt*H_spt*I_spt + Id_spt*D_spt -Fd_spt*pow2(D_spt))/(H_spt*12.*r)*(6*r + 16*pow3(r) - 6*pow5(r) + 3*pow3(r*r-1)*log((1+r)/fabs(1-r))) - (D_spt*KLd_spt + D_spt*rempf_eft*fdgp_spt*H_spt*KL_spt)/(H_spt*24.*pow3(r))*(-6*r + 22*pow3(r) + 22*pow5(r) - 6*pow7(r) + 3*pow4(r*r-1)*log((1+r)/fabs(1-r)));

return RENW(k*r) * s;
}

//return  RENW(k*r) * s;
//}

/* P_{\delta\theta}^{(13)} */
real EFT::Pres13_dt(real k) const {
   	real KMAX = QMAXp/k;
	  real KMIN = QMINp/k;
    return k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt) *RENW(k) * Integrate<ExpSub>(bind(Df13_dt, k, _1), KMIN, KMAX,  error_eft,  error_eft);
}


static real Df22_tt(real k, real r, real x) {
    real d = 1. + r*r - 2.*r*x;
    if(d < 1e-5)
        return 0;
    else
        return  RENW(k*r) * RENW(k*sqrt(d)) * r * r * 2 * ( pow2(rempf_eft*fdgp_spt*pow2(D_spt)) * pow2((7*x - r - 6*r*x*x)/(14*r*d)) + pow2(Fd_spt/H_spt)*pow2((x*x-1)/d) + (rempf_eft*fdgp_spt*pow2(D_spt)*Fd_spt/H_spt)*(r-7*x+6*x*x*r)*(x*x-1)/(7*r*pow2(d)));}


static real Dmidintt(real k, real r) {
    real KMAX = QMAXp/k;
    real KMIN = QMINp/k;
    real YMIN = Max(-1., (1.+r*r-KMAX*KMAX)/2./r);
    real YMAX = Min(0.99999999999999, (1.+r*r-KMIN*KMIN)/2./r);
return Integrate(bind(Df22_tt, k, r, _1), YMIN, YMAX , error_eft);
}

/* P_{\theta\theta}^{(22)} */
real EFT::Pres22_tt(real k) const {
	real KMAX = QMAXp/k;
	real KMIN = QMINp/k;
    return k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt) * Integrate<ExpSub>(bind(Dmidintt, k, _1), KMIN, KMAX, error_eft);
}


static real Df13_tt(real k, real r) {
real s;

if(r < 1e-2)

s =  pow2(rempf_eft*fdgp_spt*pow2(D_spt))/84.*(-56 - (32./5.)*pow2(r) - (96./7.)*pow4(r) + (352./105.)*pow6(r)) +  rempf_eft*fdgp_spt*D_spt*Cdx_spt/(H_spt*6.) * (-96./5. * r*r + 96./35. * pow4(r) + 32./105.*pow6(r)) + rempf_eft*fdgp_spt*D_spt*(Id_spt-Fd_spt*D_spt)/(6.*H_spt) * (32*r*r - 96./5. * pow4(r) + 96./35. * pow6(r)) + rempf_eft*fdgp_spt*D_spt*KLd_spt/(H_spt*12.) * ( 256./5.*r*r - 768./35. * pow4(r) + 256./105. * pow6(r)) ;

else if(fabs(r-1) < 1e-10)

s =  pow2(rempf_eft*fdgp_spt*pow2(D_spt))/84.*(-72 - 40*(r-1)) +   rempf_eft*fdgp_spt*D_spt*Cdx_spt/(H_spt*6.) * (-16 - 24*(r-1)) + rempf_eft*fdgp_spt*D_spt*(Id_spt-Fd_spt*D_spt)/(6.*H_spt) *(16 + 8*(r-1)) +  rempf_eft*fdgp_spt*D_spt*KLd_spt/(H_spt*12.) *(32 + 32*(r-1)) ;

else if(r>100)
s =  pow2(rempf_eft*fdgp_spt*pow2(D_spt))/84.*(-504./5. + (1248./35.)/pow2(r) - (608./105.)/pow4(r) - (160./231.)/pow6(r)) + rempf_eft*fdgp_spt*D_spt*Cdx_spt/(H_spt*6.)*(-32 + 96./(5.*r*r) - 96./(35.*pow4(r)) - 32./(105.*pow6(r)))+ rempf_eft*fdgp_spt*D_spt*(Id_spt-Fd_spt*D_spt)/(6.*H_spt) *(96./5. - 96./(35.*pow2(r)) - 32./(105.*pow4(r)) - 32./(385.*pow6(r)) ) + rempf_eft*fdgp_spt*D_spt*KLd_spt/(H_spt*12.)  * (256./5. - 768./(35.*pow2(r)) + 256./(105.*pow4(r)) + 256./(1155.*pow6(r))) ;

else
s =  pow2(rempf_eft*fdgp_spt*pow2(D_spt))/84.*(12./pow2(r) - 82 + 4*pow2(r) - 6*pow4(r) + 3/pow3(r) * pow3(r*r - 1) * (r*r + 2) * log((1+r)/fabs(1-r))) + rempf_eft*fdgp_spt*D_spt*Cdx_spt/H_spt * 1./(6.*pow3(r))*(6*r - 16*pow3(r) - 6*pow5(r) + 3*pow3(r*r-1)* log((1+r)/fabs(1-r))) + rempf_eft*fdgp_spt*D_spt*(Id_spt-Fd_spt*D_spt)/(H_spt*6.*r)*(6*r + 16*pow3(r) - 6*pow5(r) + 3*pow3(r*r-1)*log((1+r)/fabs(1-r))) + rempf_eft*fdgp_spt*D_spt*(KLd_spt)/(H_spt*12.*pow3(r))*(-6*r + 22*pow3(r) + 22*pow5(r) - 6*pow7(r) + 3*pow4(r*r-1)* log((1+r)/fabs(1-r))) ;

return RENW(k*r) * s;}




/* P_{\theta\theta}^{(13)} */

real EFT::Pres13_tt(real k) const {
  real KMAX = QMAXp/k;
	real KMIN = QMINp/k;
    return k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt) *RENW(k) * Integrate<ExpSub>(bind(Df13_tt, k, _1), KMIN, KMAX, error_eft,error_eft);
}


/* No Wiggle 1-loop power spectra */

double EFT::PNWloop(double k, int a) const{
  switch(a) {
      case 1:
          return   pow2(sig_8_eft*D_spt/dnorm_spt)*RENW(k) + pow4(sig_8_eft)*(Pres13_dd(k) + Pres22_dd(k));
          break;
      case 2:
          return  -rempf_eft*fdgp_spt*pow2(sig_8_eft*D_spt/dnorm_spt)*RENW(k) + pow4(sig_8_eft)*(Pres13_dt(k) + Pres22_dt(k));
          break;
      case 3:
          return  pow2(sig_8_eft*rempf_eft*fdgp_spt*D_spt/dnorm_spt)*RENW(k) + pow4(sig_8_eft)*(Pres13_tt(k) + Pres22_tt(k));
          break;
      default:
          warning("SPT: invalid indices, a = %d\n", a);
          return 0;
  }
}

/* Resummed 1-loop spectra a la Fonseca */

double EFT::Presum(double k, int a) const{
  SPT spt(C,P_l,1e-3);
  switch(a){
  case 1:
   return   PNWloop(k,1) + exp(-0.5*pow2(sig_8_eft*k)*ANW)*(spt.PLOOP(k,4)-PNWloop(k,1) + 0.5*pow2(sig_8_eft*k)*ANW*(pow2(sig_8_eft*D_spt/dnorm_spt)*(P_l(k) - RENW(k))));
      break;
  case 2:
  //  return   PNWloop(k,2) + exp(-0.5*pow2(sig_8_eft*k)*ANW)*(spt.PLOOP(k,1)-PNWloop(k,2) + 0.5*pow2(sig_8_eft*k)*ANW*(pow2(sig_8_eft*D_spt/dnorm_spt)*(P_l(k) - RENW(k))));
   return   PNWloop(k,2) + exp(-0.5*pow2(sig_8_eft*k)*ANW)*(spt.PLOOP(k,5)-PNWloop(k,2) + 0.5*pow2(sig_8_eft*k)*ANW*(-rempf_eft*fdgp_spt*pow2(sig_8_eft*D_spt/dnorm_spt)*(P_l(k) - RENW(k))));
      break;
  case 3:
      return PNWloop(k,3) + exp(-0.5*pow2(sig_8_eft*k)*ANW)*(spt.PLOOP(k,6)-PNWloop(k,3) + 0.5*pow2(sig_8_eft*k)*ANW*(pow2(rempf_eft*fdgp_spt)*pow2(sig_8_eft*D_spt/dnorm_spt)*(P_l(k) - RENW(k))));
      break;
  case 4:
      return pow2(sig_8_eft*D_spt/dnorm_spt)*RENW(k);
  case 5:
      return exp(-0.5*pow2(sig_8_eft*k)*ANW)*pow2(sig_8_eft*D_spt/dnorm_spt)*(P_l(k)-RENW(k));
  default:
          warning("SPT: invalid indices, a = %d\n", a);
          return 0;
}
}



/* EFT-Resummed multipole rsd spectra a la Fonseca */

// u1<1 and e=0: 2d-spectrum with u1 = u (not correct currently - need to decompose the counter term into cs0,cs2,cs4)
// u1 > 1 and 0<e<4 : multipoles
// e decides multipole
// e=1 : mono
// e=2 : quad
// e=3 : hex
// barr[0]=b1, barr[1]=b2, barr[2]=N

real EFT::PTNSeft(real k, double barr[], real u1, int e, double ds1) const {
    SPT spt(C,P_l,1e-3);
    double sv, nw_part1,nw_part2,nw_part3, w_part1, w_part2, f1,f2,f3,f4,f5,f6,f7,f8,f9,c10,p1,p2,p3,abc1,abc2,abc3, lbnw1, lbnw2;
    double F0 = rempf_eft*fdgp_spt;
    sv = pow2(k)*sigma_v; // = (ksigma_v)^2
    if (u1<1.) {
      f1 = (1-sv*pow2(u1*F0));
      f2 = pow2(u1)*f1;
      f3 = pow4(u1)*f1;
      f4 = f1 * exp(-0.5*pow2(k)*pow2(barr[0])*ANW*(1+F0*(2+F0)*pow2(u1)));
      f5 = pow2(u1)*f4;
      f6 = pow4(u1)*f4;
      f7 = f4*pow2(k)*pow2(barr[0])*ANW*(1+F0*(2+F0)*pow2(u1));
      f8 = pow2(u1)*f7;
      f9 = pow4(u1)*f7;
      c10 = exp(-0.5*pow2(k)*pow2(barr[0])*ANW*(1+F0*(2+F0)*pow2(u1)));
      abc1 = ABCnow(k, barr[0], u1, pow2(barr[0])*ANW, 1, 1, e); // ABC w no prefactor w no wiggle spectrum
      abc2 = ABCnow(k, barr[0], u1, pow2(barr[0])*ANW, 2, 1, e); // ABC w exp(-0.5pow2(barr[0])*ANW_s) w no wiggle spectrum
      abc3 = ABCnow(k, barr[0], u1, pow2(barr[0])*ANW, 2, 2, e); // ABC w exp(-0.5pow2(barr[0])*ANW_s) w P_l spectrum
    }
    else{
    f1 = factL(k,sv,F0,pow2(barr[0])*ANW,0,e,1);
    f2 = factL(k,sv,F0,pow2(barr[0])*ANW,1,e,1);
    f3 = factL(k,sv,F0,pow2(barr[0])*ANW,2,e,1);
    f4 = factL(k,sv,F0,pow2(barr[0])*ANW,0,e,2);
    f5 = factL(k,sv,F0,pow2(barr[0])*ANW,1,e,2);
    f6 = factL(k,sv,F0,pow2(barr[0])*ANW,2,e,2);
    f7 = factL(k,sv,F0,pow2(barr[0])*ANW,0,e,3);
    f8 = factL(k,sv,F0,pow2(barr[0])*ANW,1,e,3);
    f9 = factL(k,sv,F0,pow2(barr[0])*ANW,2,e,3);
    c10 = 1.;
    abc1 = ABCnow(k, barr[0], 0., pow2(barr[0])*ANW, 1, 1, e); // ABC w no prefactor w no wiggle spectrum
    abc2 = ABCnow(k, barr[0], 0., pow2(barr[0])*ANW, 2, 1, e); // ABC w exp(-0.5pow2(barr[0])*ANW_s) w no wiggle spectrum
    abc3 = ABCnow(k, barr[0], 0., pow2(barr[0])*ANW, 2, 2, e); // ABC w exp(-0.5pow2(barr[0])*ANW_s) w P_l spectrum
    }

    p1 = PNWloop(k,1);
    p2 = PNWloop(k,2);
    p3 = PNWloop(k,3);

    lbnw1 = Lag_bias_nw(1,k,barr);
    lbnw2 = Lag_bias_nw(2,k,barr);

    nw_part1 =  f1*(pow2(barr[0])*p1 + lbnw1 + barr[2])  - 2.*f2*(barr[0]*p2 + lbnw2) + f3*p3 + abc1; // no wiggle part 1
    w_part1  =  f4*(pow2(barr[0])*(spt.PLOOP(k,4) - p1) + spt.Lag_bias(1,k,barr) - lbnw1 + barr[2]) - 2*f5*(barr[0]*(spt.PLOOP(k,5) - p2) + spt.Lag_bias(2,k,barr) - lbnw2) + f6*(spt.PLOOP(k,6)-p3) + c10*(abc3-abc2); // wiggle part
    w_part2  =  pow2(sig_8_eft*D_spt/dnorm_spt)*(pow2(barr[0])*f7*(P_l(k)-RENW(k)) +  2.*f8*(F0*barr[0])*(P_l(k)-RENW(k)) + pow2(F0)*f9*(P_l(k)-RENW(k)));

    return nw_part1 + w_part1 + w_part2 - 2.*ds1*pow2(k)*P_l(k)*pow2(sig_8_eft*D_spt/dnorm_spt);
//  pow2(bl)*factL(k,sv,0,a,2)*Presum(k,1) - 2*factL(k,sv,1,a,2)*bl*Presum(k,2) + factL(k,sv,2,a,2)* Presum(k,3) +  ABCnow(k, bl, 0., 2, a) - 2.*ds1*pow2(k)*RENW(k)*pow2(D_spt/dnorm_spt);
}



/* TNS power spectrum in SPT (sigma v calculated within SPT) */

real EFT::PTNSspt(real k, double barr[], real u1, int e) const {
    SPT spt(C,P_l,1e-3);
    double sv,f4,f5,f6,abc;
    double F0 = rempf_eft*fdgp_spt;
    sv = pow2(k)*sigma_v;//(k); // = (ksigma_v)^2
    if (u1<1) {
      f4 = (1.-sv*pow2(u1*F0));
      f5 = pow2(u1)*f4;
      f6 = pow4(u1)*f4;
      abc = ABCnow(k,  barr[0], u1, 0., 1, 2, e);
    }
    else{
    f4 = factL(k,sv,F0,0.,0,e,1);
    f5 = factL(k,sv,F0,0.,1,e,1);
    f6 = factL(k,sv,F0,0.,2,e,1);
    abc = ABCnow(k, barr[0], 0., 0., 1, 2, e);
    }
    return f4*(pow2(barr[0])*spt.PLOOP(k,4) + spt.Lag_bias(1,k,barr) + barr[2]) - 2.*f5*(barr[0]*spt.PLOOP(k,5) + spt.Lag_bias(2,k,barr)) + f6*spt.PLOOP(k,6) + abc;
}


/* Matsubara version of the SPT rsd power spectrum (should be the same as the TNS version ) */
// u=0, e=1,2,3 gives multipoles
//u neq 0 , e=0 gives 2d spectrum
real EFT::PTNSspt_mat(real k, real bl, real u1, int e) const {
    double f4,nlm;
    double F0 = rempf_eft*fdgp_spt;
    if (u1<1) {
      f4 = pow2(1+pow2(u1)*F0);
      nlm = NL_mat(k, bl, u1, 2, e);
    }
    else{
    f4 = factL(k,0.,F0,0.,0,e,4);
    nlm = NL_mat(k, bl, 0., 2, e);
    }
    return pow2(sig_8_eft*D_spt/dnorm_spt)*f4*P_l(k) + nlm;
}




/* TAYLOR EXPANSION BIAS MODEL terms: NO Wiggle version (used for resummation)*/

static real lbias_selec_nw(int a,  real bias[], real k, real r, real x) {
    real d = 1. + r*r - 2.*r*x;
    if(d < 1e-5)
        return 0;
    else{
      double sker = 2./3.-(1.-pow2((x-r)/sqrt(d)));
      double dker = -2./7.*(1.-x*x);
      double plkr = RENW(r*k);
      double bs2 = -4./7.*(bias[0]-1.);
      double b3n = 32./315.*(bias[0]-1.);
      double G2ker = G2eds(k*r,k*sqrt(d),(x-r)/sqrt(d));
      double F2ker = F2eds(k*r,k*sqrt(d),(x-r)/sqrt(d));
      switch (a) {
      case 1:
      return  plkr*RENW(k*sqrt(d))*pow2(r)*
              (2.*bias[0]*bias[1]*F2ker // b2
              +2.*bias[0]*bs2*F2ker*sker //bs2
              +pow2(bias[1])/2. //b22 a
              +bias[1]*bs2*sker //b2s2 a
              +pow2(bs2)*pow2(sker)/2.) //bs22 a
              - pow2(r*plkr)*(pow2(bias[1])/2. //b22 b
                                +2.*bias[1]*bs2/3. // b2s2 b
                                +pow2(bs2)*2./9.) // bs22 b
              +RENW(k)*b3n*bias[0]*plkr*2.*105./16.*pow2(r)*(dker*sker+8./63.) // b3 terms
                ;  // shot noise
     case 2 :
     return   -rempf_eft*fl_spt*plkr*RENW(k*sqrt(d))*pow2(r)*(bias[1]*G2ker //b2t
                                                 +bs2*G2ker*sker) // bs2t
              -rempf_eft*fl_spt*RENW(k)*b3n*plkr*105./16.*pow2(r)*(dker*sker+8./63.);
            }
}
}

/* Integrate them */
static real midpb2_nw(int a, real bias[], real k, real r) {
    real KMAX = QMAXp/k;
    real KMIN = QMINp/k;
    real YMIN = Max(-1., (1.+r*r-KMAX*KMAX)/2./r);
    real YMAX = Min(0.9999999999999, (1.+r*r-KMIN*KMIN)/2./r);
return Integrate(bind(lbias_selec_nw, a, bias, k, r, _1), YMIN, YMAX ,  error_eft);
}

//P_bs2
real EFT::Lag_bias_nw(int a, real k, real bias[]) const {
 	real KMAX = QMAXp/k;
	real KMIN = QMINp/k;
    return k*k*k/(4*M_PI*M_PI) * pow4(sig_8_eft*Dl_spt/dnorm_spt) * Integrate<ExpSub>(bind(midpb2_nw, a, bias, k, _1), KMIN, KMAX, error_eft);
}



/* RSD CORRECTION TERMS: Non vanishing A, B and C terms and DGP C terms */

// Notes for analytic ABC terms : DEFINITION OF THETA= - DEL V / aHf (and is accounted for in definition of Cross Bispectrum)
// See http://arxiv.org/pdf/1006.0699v1.pdf for derivation

// NOTE ON PARAMETERS: x =k1.k/(k^2*r)

// c chooses wiggle or non wiggle linear spectrum

static real ABC(int c, int a, const PowerSpectrum& P_l, real k, real r, real x) {
  double pk1,pk2,pk3;
  real d = 1. + r*r - 2.*r*x;
  if(d < 1e-5)
        return 0;
  else{
  if (c==1) {
    pk1 = RENW(k);
    pk2 = RENW(k*sqrt(d));
    pk3 = RENW(k*r);
  }
  else{
    pk1 = P_l(k);
    pk2 = P_l(k*sqrt(d));
    pk3 = P_l(k*r);
  }
		switch(a) {
		// A terms
				//A11
			case 1:
				return (-r*r*r/7.)*(x+6*x*x*x+r*r*x*(-3+10*x*x)+r*(-3+x*x-12*x*x*x*x))*pk1 * pk2 / pow2(d);
        break;
      	//A12
			case 2:
				return  (r*r*r*r/14.)*(x*x-1)*(-1+7*r*x-6*x*x) * pk1 * pk2 / pow2(d);
        break;
      	//A22
			case 3:
				return  (r*r*r/14.)*(r*r*x*(13-41*x*x)-4*(x+6*x*x*x)+r*(5+9*x*x+42*x*x*x*x))*pk1*pk2/ pow2(d);
        break;
      	//A33
			case 4:
				return (r*r*r/14.)*(1-7*r*x+6*x*x)*(-2*x+r*(-1+3*x*x))*pk1*pk2/ pow2(d);
        break;
      	//A11t
			case 5:
				return  (1./7.)*(x+r-2*r*x*x)*(3*r+7*x-10*r*x*x)*pk3*pk2/pow2(d);
        break;
      	//A12t
			case 6:
				return (r/14.)*(x*x-1)*(3*r+7*x-10*r*x*x)*pk3*pk2/pow2(d);
        break;
      	//A22t
			case 7:
				return (1./14.)*(28*x*x+r*x*(25-81*x*x)+r*r*(1-27*x*x+54*x*x*x*x))*pk3*pk2/pow2(d);
        break;
      	//A23t
			case 8:
				return (r/14.)*(-x*x+1)*(r-7*x+6*r*x*x)*pk3*pk2/pow2(d);
        break;
		//A33t
			case 9:
				return (1./14.)*(-2*x-r+3*r*x*x)*(r-7*x+6*r*x*x)*pk3*pk2/pow2(d);
        break;
		//a11
			case 10:
				return (-7*x*x+pow3(r)*x*(-3+10*x*x)+3*r*(x+6*pow3(x))+r*r*(6-19*x*x-8*pow4(x)))*pk1*pk3/(7.*d);
        break;
		//a12
			case 11:
				return r*(-1+x*x)*(6*r-7*(1+r*r)*x + 8*r*x*x)*pk1*pk3/(14.*d);
        break;
			//a22
			case 12:
				return (-28*x*x+r*r*r*x*(-13+41*x*x)+r*x*(11+73*x*x)-2*r*r*(-9+31*x*x+20*x*x*x*x))*pk1*pk3/(14.*d);
        break;
		//a33
			case 13:
				return (7*x + r*(-6+7*r*x-8*x*x))*(-2*x+r*(-1+3*x*x))*pk1*pk3/(14.*d);
        break;

		//B terms
				//B111
			case 14:
				return (r*r/2.)*(x*x-1)*pk3*pk2/d;
        break;
	//B112
			case 15:
				return (3*r*r/8.)*pow2(x*x-1)*pk3*pk2/d;
        break;
		//B121
			case 16:
				return (3.*r*r*r*r/8.)*pow2(x*x-1)*pk3*pk2/pow2(d);
        break;
			//B122
			case 17:
				return (5*r*r*r*r/16.)*pow3(x*x-1)*pk3*pk2/pow2(d);
        break;
			//B211
			case 18:
				return (r/2.)*(r+2*x-3*r*x*x)*pk3*pk2/d;
        break;
		//B212
			case 19:
				return (-3*r/4.)*(x*x-1)*(-r-2*x+5*r*x*x)*pk3*pk2/d;
        break;
			//B221
			case 20:
				return (3*r*r/4.)*(x*x-1)*(-2+r*r+6*r*x-5*r*r*x*x)*pk3*pk2/pow2(d);
        break;
			//B222
			case 21:
				return (-3*r*r/16.)*pow2(x*x-1)*(6-30*r*x-5*r*r+35*r*r*x*x)*pk3*pk2/pow2(d);
        break;
		//B312
			case 22:
				return (r/8.)*(4*x*(3-5*x*x)+r*(3-30*x*x+35*x*x*x*x))*pk3*pk2/d;
        break;
		//B321
			case 23:
				return (r/8.)*(-8*x+r*(-12+36*x*x+12*r*x*(3-5*x*x)+r*r*(3-30*x*x+35*x*x*x*x)))*pk3*pk2/pow2(d);
        break;
		//B322
			case 24:
				return (3*r/16.)*(x*x-1)*(-8*x+r*(-12+60*x*x+20*r*x*(3-7*x*x)+5*r*r*(1-14*x*x+21*x*x*x*x)))*pk3*pk2/pow2(d);
        break;
		//B422 // HAVE INCLUDED THE 'TYPO' IN 3RD TERM : 6rr -> 6rrx
			case 25:
				return (r/16.)*(8*x*(-3+5*x*x)-6*r*(3-30*x*x+35*x*x*x*x)+6*r*r*x*(15-70*x*x+63*x*x*x*x)+r*r*r*(5-21*x*x*(5-15*x*x+11*x*x*x*x)))*pk3*pk2/pow2(d);
		// nDGP extra terms
    break;
			//C11
			case 26:
				return 2.*((r*r*r)*(x*x-1)*(Dd_spt*F_spt*r*(r*x-1)-D_spt*Fd_spt*x*d)*pk1*pk2)/ (Dd_spt*pow2(d));
        break;
			//C12
			case 27:
				return -D_spt*Fd_spt*((r*r*r*r)*pow2(x*x-1)*pk1*pk2)/(Dd_spt*pow2(d));
        break;
			//C22
			case 28:
				return ((r*r*r)*(x*x-1)*(2*Dd_spt*F_spt*r*(r*x-1)-D_spt*Fd_spt*(r+4*x+2*r*r*x-7*r*x*x))*pk1*pk2)/ (Dd_spt*pow2(d));
        break;
			//C23
			case 29:
				return -D_spt*Fd_spt*((r*r*r*r)*pow2(x*x-1)*pk1*pk2)/ (Dd_spt*pow2(d));
        break;
				//C33
			case 30:
				return D_spt*Fd_spt*((r*r*r)*(x*x-1)*(-2*x+r*(-1+3*x*x))*pk1*pk2)/ (Dd_spt*pow2(d));
        break;
				//C11t
			case 31:
				return 2*F_spt*(r*(x*x-1)*(-x+r*(-1+2*x*x))*pk3*pk2)/ pow2(d);
        break;
				//C12t
			case 32:
				return -F_spt*(r*r*pow2(x*x-1)*pk3*pk2)/ pow2(d);
        break;
				//C22t
			case 33:
				return r*(x*x-1)*(2*D_spt*Fd_spt*(-x+r*(-1+2*x*x)) + Dd_spt*F_spt*(-2*x+r*(-1+3*x*x)))*pk3*pk2/ (Dd_spt*pow2(d));
        break;
				//C23t
			case 34:
				return -D_spt*Fd_spt*(r*r*pow2(x*x-1)*pk3*pk2)/ (Dd_spt*pow2(d));
        break;
				//C33t
			case 35:
				return D_spt*Fd_spt*(r*(x*x-1)*(-2*x+r*(-1+3*x*x))*pk3*pk2)/ (Dd_spt*pow2(d));
        break;
				//c11
			case 36:
				return  -2*r*(1-x*x)*(D_spt*Fd_spt*r*(r*x-1)-Dd_spt*F_spt*x*d)*pk3*pk1/(Dd_spt*d);
        break;
				//c12
			case 37:
				return -D_spt*Fd_spt*r*r*pow2(-1+x*x)*pk3*pk1/(Dd_spt*d);
        break;
				//c22
			case 38:
				return r*(-1+x*x)*(-2*Dd_spt*F_spt*x*d+D_spt*Fd_spt*(-2*x+2*r*r*x+3*r*(-1+x*x)))*pk3*pk1/(Dd_spt*d);
        break;
				//c33
			case 39:
				return D_spt*Fd_spt*r*(-1+x*x)*(-2*x+r*(-1+3*x*x))*pk3*pk1/(Dd_spt*d);
        break;
/* C terms for full SPT calculation */
        case 40: // u^2f^2
          return pow4(r)*(1-pow2(x))*pk3*pk2/(2*pow2(d));
      //    return (1-pow2(x))*pk3*pk2/2.;
          break;
        case 41: // u^2f^3
          return  pow4(r)*3.*(1-2*pow2(x)+pow4(x))*pk3*pk2/(4*pow2(d));
        //  return 3*(pow2(r)-2*pow2(r*x)+pow2(r*pow2(x)))*pk3*pk2/(4*d);
          break;
        case 42: // u^2f^4
        //  return 5*(pow4(r)-3*pow2(pow2(r)*x)+3*pow4(x*r)-pow4(r)*pow6(x))*pk3*pk2/(16*pow2(d));
          return pow4(r)*5.*(1-3*pow2(x)+3*pow4(x)-pow6(x))*pk3*pk2/(16*pow2(d));
          break;
        case 43: // u^4f^2
    //    return (-1+3*pow2(x))*pk3*pk2/2.;
          return pow2(r)*(16.-8*pow2(r)-32*r*x+24*pow2(x*r))*pk3*pk2/(16*pow2(d));
          break;
        case 44: // u^4f^3
    //    return (2-3*pow2(r)-12*r*x-2*pow2(x)+18*pow2(r*x)+12*r*pow3(x)-15*pow2(r*pow2(x)))*pk3*pk2/(2*d);
          return  pow2(r)*(16-24*pow2(r)-96*r*x-16*pow2(x)+144*pow2(r*x)+96*r*pow3(x)-120*pow2(r)*pow4(x))*pk3*pk2/(16*pow2(d));
          break;
        case 45: // u^4f^4
      //  return (-120*pow3(r)*x+240*pow3(r*x)-120*pow3(r)*pow5(x)+36*pow2(r)*(1-2*pow2(x)+pow4(x))+pow4(r)*(-15+135*pow2(x)-225*pow4(x)+105*pow6(x)))*pk3*pk2/(16*pow2(d));
          return  pow2(r)*(6-15*pow2(r)-60*r*x-12*pow2(x)+135*pow2(r*x)+120*r*pow3(x)+6*pow4(x)-225*pow2(r)*pow4(x)-60*r*pow5(x)+105*pow2(r)*pow6(x))*pk3*pk2/(16*pow2(d));
          break;
        case 46: // u^6f^3
    //    return (-4+3*pow2(r)+24*r*x+12*pow2(x)-30*pow2(r*x)-40*r*pow3(x)+35*pow2(r*pow2(x)))*pk3*pk2/(4*d);
          return pow2(r)*(-16+12*pow2(r)+96*r*x+48*pow2(x)-120*pow2(x*r)-160*r*pow3(x)+140*pow2(r)*pow4(x))*pk3*pk2/(16*pow2(d));
          break;
        case 47: // u^6f^4
  //      return (8-72*pow2(r)+15*pow4(r)-96*r*x+240*pow3(r)*x-8*pow2(x)+432*pow2(r*x)-225*pow2(pow2(r)*x)+96*r*pow3(x)-800*pow3(r*x)
//                -360*pow2(r*pow2(x))+525*pow4(r*x)+560*pow3(r)*pow5(x)-315*pow4(r)*pow6(x))*pk3*pk2/(16*pow2(d));
          return pow2(r)*(-12+15*pow2(r)+120*r*x+72*pow2(x)-225*pow2(x*r)-400*r*pow3(x)-60*pow4(x)+525*pow2(r)*pow4(x)+280*r*pow5(x)-315*pow2(r)*pow6(x))*pk3*pk2/(16*pow2(d));
          break;
        case 48: // u^8f^4
  //      return (-8+36*pow2(r)-5*pow4(r)+96*r*x-120*pow3(r)*x+24*pow2(x)-360*pow2(r*x)+105*pow2(pow2(r)*x)-160*r*pow3(x)+560*pow3(r*x)
    //            +420*pow2(r*pow2(x))-315*pow4(r*x)-504*pow3(r)*pow5(x)+231*pow4(r)*pow6(x))*pk3*pk2/(16*pow2(d));
         return  pow2(r)*(6-5*pow2(r)-60*r*x-60*pow2(x)+105*pow2(x*r)+280*r*pow3(x)+70*pow4(x)-315*pow2(r)*pow4(x)-252*r*pow5(x)+231*pow2(r)*pow6(x))*pk3*pk2/(16*pow2(d));
         break;
      default:
			warning("SPT: invalid indices, a = %d\n", a);
				return 0;
}
}
}


/* function to select multipole or u-dependent Redshift PS */
// it gives either u^(2n) or factL- see SpecialFunctions.cpp

// a = 1 :LCDM
// a = 2 :nDGP
// a = 3 :Numerical (Arbitrary MG)

// e = 0:  RSD PS
// e = 1 : Monopole
// e = 2 : Quadrupole
// e = 3 : Hexdecapole

static real u0(real k, real u, real F0, real anw, int n, int e,  int a){
	if (e >= 1) {
		return factL(k, u, F0, anw, n, e, a); }
		else{
			return pow(u,2.*n);
    }
}



static real midintabc(double u0x[], real bl, real anw, int a, int c, int e, const PowerSpectrum& P_l,  real k, real u, real r) {
  real KMAX = QMAXp/k;
  real KMIN = QMINp/k;
	real YMIN = Max(-1., (1.+r*r-KMAX*KMAX)/2./r);
  real YMAX1 = Min(0.99999999999, (1.+r*r-KMIN*KMIN)/2./r);
  real YMAX = ATS(k,r);
	real F0,D1,X1;
			F0 = rempf_eft*fdgp_spt;
			D1 = sig_8_eft*D_spt;
			X1 = 1.;
      real u01 = u0x[0];//u0(k,u,F0,anw,1,e,a);
      real u02 = u0x[1];//u0(k,u,F0,anw,2,e,a);
      real u03 = u0x[2];//u0(k,u,F0,anw,3,e,a);
      real u04 = u0x[3];//u0(k,u,F0,anw,4,e,a);
	//A term
		return     pow4(D1/dnorm_spt)*2.*(F0*u01*pow2(bl)*(Integrate(bind(ABC, c,  1,  cref(P_l), k, r, _1),  YMIN, YMAX ,  error_eft)
										   + Integrate(bind(ABC, c,  5,  cref(P_l), k,  r, _1),  YMIN, YMAX ,  error_eft)
										   + Integrate(bind(ABC, c,  10,  cref(P_l), k, r, _1),  YMIN, YMAX ,  error_eft))

							+ F0*F0*u01*bl*(Integrate(bind(ABC, c,  2,  cref(P_l), k,  r, _1),  YMIN, YMAX , error_eft)
										   + Integrate(bind(ABC, c,  6,  cref(P_l), k,  r, _1), YMIN, YMAX ,  error_eft)
										   + Integrate(bind(ABC, c,  11,  cref(P_l), k,  r, _1), YMIN, YMAX ,  error_eft))

							+ F0*F0*u02*bl*(Integrate(bind(ABC, c,  3,  cref(P_l), k,  r, _1),   YMIN, YMAX , error_eft)
										   + Integrate(bind(ABC, c,  7,  cref(P_l), k,  r, _1),   YMIN, YMAX ,  error_eft)
										   + Integrate(bind(ABC, c,  12,  cref(P_l), k,  r, _1),  YMIN, YMAX ,  error_eft))

					     + F0*F0*F0*u02*(Integrate(bind(ABC, c,  2,  cref(P_l), k,  r, _1),  YMIN, YMAX ,  error_eft)
										   + Integrate(bind(ABC, c,  8,  cref(P_l), k,  r, _1),   YMIN, YMAX ,  error_eft)
										   + Integrate(bind(ABC, c,  11,  cref(P_l), k, r, _1), YMIN, YMAX ,  error_eft))

						   + F0*F0*F0*u03*(Integrate(bind(ABC, c,  4,  cref(P_l), k, r, _1),   YMIN, YMAX , error_eft)
										   + Integrate(bind(ABC, c,  9,  cref(P_l), k,  r, _1),  YMIN, YMAX ,  error_eft)
										   + Integrate(bind(ABC, c,  13,  cref(P_l), k,  r, _1),YMIN, YMAX ,  error_eft)))

	//B term
		  +  pow4(D1/dnorm_spt)*2*(F0*F0*pow2(bl)*u01*Integrate(bind(ABC, c,  14,   cref(P_l), k, r, _1), YMIN, YMAX , error_eft)

						  -F0*F0*F0*bl*u01*(Integrate(bind(ABC, c,  15,  cref(P_l), k, r, _1), YMIN, YMAX , error_eft)
											    +Integrate(bind(ABC, c,  16,   cref(P_l), k,  r, _1),YMIN, YMAX , error_eft))

						+F0*F0*F0*F0*u01*Integrate(bind(ABC, c,  17,   cref(P_l), k, r, _1), YMIN, YMAX , error_eft)

							  +F0*F0*pow2(bl)*u02*Integrate(bind(ABC, c,  18,   cref(P_l), k,  r, _1),YMIN, YMAX ,  error_eft)

						 -F0*F0*F0*bl*u02*(Integrate(bind(ABC, c,  19,  cref(P_l), k,  r, _1), YMIN, YMAX , error_eft)
										      + Integrate(bind(ABC, c,  20,   cref(P_l), k,  r, _1), YMIN, YMAX ,  error_eft))

					  +F0*F0*F0*F0*u02*Integrate(bind(ABC, c,  21,   cref(P_l), k,  r, _1),YMIN, YMAX , error_eft)

					    - F0*F0*F0*bl*u03*(Integrate(bind(ABC, c,  22,   cref(P_l), k,  r, _1), YMIN, YMAX ,  error_eft)
											  + Integrate(bind(ABC, c,  23,   cref(P_l), k,  r, _1), YMIN, YMAX , error_eft))

					  + F0*F0*F0*F0*u03*Integrate(bind(ABC, c,  24,   cref(P_l), k,  r, _1),YMIN, YMAX ,  error_eft)
							+F0*F0*F0*F0*u04*Integrate(bind(ABC, c,  25,   cref(P_l), k,  r, _1), YMIN, YMAX , error_eft))

// C term- full SPT calculation @ tree level

      + pow4(D1/dnorm_spt)*pow2(bl)*(pow2(F0)*u01*Integrate(bind(ABC, c,  40,  cref(P_l), k, r, _1), YMIN, YMAX1 , error_eft)
                                  +pow3(F0)*u01*Integrate(bind(ABC, c,  41,  cref(P_l), k, r, _1), YMIN, YMAX1 , error_eft)
                                  +pow4(F0)*u01*Integrate(bind(ABC, c,  42,  cref(P_l), k, r, _1), YMIN, YMAX1 , error_eft)
                                  +pow2(F0)*u02*Integrate(bind(ABC, c,  43,  cref(P_l), k, r, _1), YMIN, YMAX1 , error_eft)
                                  +pow3(F0)*u02*Integrate(bind(ABC, c,  44,  cref(P_l), k, r, _1), YMIN, YMAX1 , error_eft)
                                  +pow4(F0)*u02*Integrate(bind(ABC, c,  45,  cref(P_l), k, r, _1), YMIN, YMAX1 , error_eft)
                                  +pow3(F0)*u03*Integrate(bind(ABC, c,  46,  cref(P_l), k, r, _1), YMIN, YMAX1 , error_eft)
                                  +pow4(F0)*u03*Integrate(bind(ABC, c,  47,  cref(P_l), k, r, _1), YMIN, YMAX1 , error_eft)
                                  +pow4(F0)*u04*Integrate(bind(ABC, c,  48,  cref(P_l), k, r, _1), YMIN, YMAX1 , error_eft))
//nDGP analytic

		  +       X1*D1*D1*2/pow4(dnorm_spt)*(F0*pow2(bl)*u01*(Integrate(bind(ABC, c,  26,  cref(P_l), k, r, _1),   YMIN, YMAX ,  error_eft)
										+ Integrate(bind(ABC, c,  31,  cref(P_l), k, r, _1),   YMIN, YMAX ,  error_eft)
										+ Integrate(bind(ABC, c,  36,  cref(P_l), k, r, _1), YMIN, YMAX ,  error_eft)) +

							F0*F0*bl*u01*(Integrate(bind(ABC, c,  27,  cref(P_l), k, r, _1),   YMIN, YMAX ,  error_eft)
										+ Integrate(bind(ABC, c,  32,  cref(P_l), k, r, _1),   YMIN, YMAX ,  error_eft)
										+ Integrate(bind(ABC, c,  37,  cref(P_l), k, r, _1),  YMIN, YMAX ,  error_eft)) +

						   F0*F0*bl*u02*(Integrate(bind(ABC, c,  28,  cref(P_l), k, r, _1),   YMIN, YMAX ,  error_eft)
										+ Integrate(bind(ABC, c,  33,  cref(P_l), k, r, _1),   YMIN, YMAX ,  error_eft)
										+ Integrate(bind(ABC, c,  38,  cref(P_l), k,  r, _1), YMIN, YMAX ,  error_eft)) +

						F0*F0*F0*u02*(Integrate(bind(ABC, c,  29,  cref(P_l), k, r, _1),   YMIN, YMAX ,  error_eft)
										+ Integrate(bind(ABC, c,  34,  cref(P_l), k, r, _1),   YMIN, YMAX ,  error_eft)
										+ Integrate(bind(ABC, c,  37,  cref(P_l), k,  r, _1),  YMIN, YMAX ,  error_eft)) +

					   F0*F0*F0*u03*(Integrate(bind(ABC, c,  30,  cref(P_l), k, r, _1),   YMIN, YMAX ,  error_eft)
									   + Integrate(bind(ABC, c,  35,  cref(P_l), k, r, _1),   YMIN, YMAX ,  error_eft)
									   + Integrate(bind(ABC, c,  39,  cref(P_l), k,  r, _1),  YMIN, YMAX ,  error_eft)));

}

/*Analytical A and B and C*/
// e = 0:  RSD PS
// e = 1 : Monopole
// e = 2 : Quadrupole
// e = 3 : Hexdecapole

// c = 0:  P_l
// c = 1 : RENW

//a gives multipole factor

//U = u for 2d PS
//U = sigma_v for multipoles

// For A or B or C seperated : just remove as necessary.

real EFT::ABCnow(real k, real bl, real U, real anw, int a, int c, int e) const{
      int n3 = 300;
      real KMAX = QMAXp/k;
      real KMIN = QMINp/k;
      double y[n3];
      double integrand[n3];
      double u0x[4];
      double F0 = rempf_eft*fdgp_spt;
      u0x[0] = u0(k,U,F0,anw,1,e,a);
      u0x[1] = u0(k,U,F0,anw,2,e,a);
      u0x[2] = u0(k,U,F0,anw,3,e,a);
      u0x[3] = u0(k,U,F0,anw,4,e,a);
      for (int i = 0; i<n3; i++){
      y[i] = KMIN * exp(i*log(KMAX/KMIN)/(n3-1.));
      integrand[i] = midintabc(u0x,bl, anw, a, c, e, cref(P_l), k, U, y[i]);
      }
    double res = 0.;
      for( int i = 1; i < n3; ++ i ){
    res += 0.5 * (y[i] - y[i-1])*(integrand[i] + integrand[i-1]);
        }
      return  k*k*k/(4*M_PI*M_PI) * res;
      }




static real Mterms(int c, int a, const PowerSpectrum& P_l, real k, real r, real x) {
        double pk1,pk2,pk3;
        real d = 1. + r*r - 2.*r*x;
        if(d < 1e-5)
              return 0;
        else{
        if (c==1) {
          pk1 = RENW(k);
          pk2 = RENW(k*sqrt(d))/pow2(d);
          pk3 = RENW(k*r);
        }
        else{
          pk1 = P_l(k);
          pk2 = P_l(k*sqrt(d))/pow2(d);
          pk3 = P_l(k*r);
        }
      		switch(a) {
      		// A terms
      			case 1:
      				return 1./98.*pow2(3*r+7*x-10*r*pow2(x))*pk2*pk3;
              break;
            case 2:
              return 1./28.*(1-pow2(x))*(7-6*pow2(r)-42*r*x+48*pow2(x*r))*pk2*pk3;
              break;
            case 3:
              return 1./196.*(-49+637*pow2(x)+42*r*x*(17-45*pow2(x))+6*pow2(r)*(19-157*pow2(x)+236*pow4(x)))*pk2*pk3;
              break;
            case 4:
              return 3./16.*pow2(r)*pow2(1-pow2(x))*pk2*pk3;
              break;
            case 5:
              return 1./14.*(-7+35*pow2(x)+54*r*x-110*r*pow3(x)+6*pow2(r)-66*pow2(r*x)+88*pow2(r*pow2(x)))*pk2*pk3;
              break;
            case 6:
             return 1./8.*(1-pow2(x))*(2-3*pow2(r)-12*r*x+15*pow2(r*x))*pk2*pk3;
              break;
            case 7:
              return 1./16.*(-4+12*pow2(x)+24*r*x-40*r*pow3(x)+3*pow2(r)-30*pow2(r*x)+35*pow2(r*pow2(x)))*pk2*pk3;
              break;
      		//B terms
      			case 8:
      				return 1./252.*(12/pow2(r)-158+100*pow2(r)-42*pow4(r)+3/pow3(r)*pow3(pow2(r)-1)*(7*pow2(r)+2)*log((1+r)/fabs(1-r)))*pk3*pk1;
              break;
            case 9:
              return 1./168.*(18/pow2(r)-178-66*pow2(r)+18*pow4(r)-9/pow3(r)*pow4(pow2(r)-1)*log((1+r)/fabs(1-r)))*pk3*pk1;
              break;
            case 10:
              return 1./168.*(18/pow2(r)-218+126*pow2(r)-54*pow4(r)+9/pow3(r)*pow3(pow2(r)-1)*(3*pow2(r)+1)*log((1+r)/fabs(1-r)))*pk3*pk1;
              break;
            case 11:
              return -2./3.*pk3*pk1;
              break;

            default:
      			warning("EFT: invalid indices, a = %d\n", a);
      				return 0;
      }
      }
    }



 static real midintmat(double u0x[], real bl, int c, int e, const PowerSpectrum& P_l,  real k, real u, real r) {
      real KMAX = QMAXp/k;
      real KMIN = QMINp/k;
    	real YMIN = Max(-1., (1.+r*r-KMAX*KMAX)/2./r);
      real YMAX = Min(0.99999999999, (1.+r*r-KMIN*KMIN)/2./r);
    	real F0,D1;
      real u00 = u0x[0];//u0(k,u,F0,0.,0,e,1);
      real u01 = u0x[1];//u0(k,u,F0,0.,1,e,1);
      real u02 = u0x[2];//u0(k,u,F0,0.,2,e,1);
      real u03 = u0x[3];//u0(k,u,F0,0.,3,e,1);
      real u04 = u0x[4];//u0(k,u,F0,0.,4,e,1);

      real u10 = u0x[5];//u0(k,u,F0,0.,0,e,5);
      real u11 = u0x[6];//u0(k,u,F0,0.,1,e,5);
      real u12 = u0x[7];//u0(k,u,F0,0.,2,e,5);
    			F0 = rempf_eft*fdgp_spt;
    			D1 = pow4(sig_8_eft*D_spt/dnorm_spt);
    	//A term
    		return     D1*(u00*Integrate(bind(Mterms, c,  1,  cref(P_l), k, r, _1),  YMIN, YMAX ,  error_eft)
    					      + 4.*F0*u01*Integrate(bind(Mterms, c,  1,  cref(P_l), k,  r, _1),  YMIN, YMAX ,  error_eft)
    						+ pow2(F0)*u01*Integrate(bind(Mterms, c,  2,  cref(P_l), k, r, _1),  YMIN, YMAX ,  error_eft)
                + pow2(F0)*u02*Integrate(bind(Mterms, c,  3,  cref(P_l), k, r, _1),  YMIN, YMAX ,  error_eft)
              + 2.*pow3(F0)*u02*Integrate(bind(Mterms, c,  2,  cref(P_l), k, r, _1),  YMIN, YMAX ,  error_eft)
                + pow4(F0)*u02*Integrate(bind(Mterms, c,  4,  cref(P_l), k, r, _1),  YMIN, YMAX ,  error_eft)
                + pow3(F0)*u03*Integrate(bind(Mterms, c,  5,  cref(P_l), k, r, _1),  YMIN, YMAX ,  error_eft)
                + pow4(F0)*u03*Integrate(bind(Mterms, c,  6,  cref(P_l), k, r, _1),  YMIN, YMAX ,  error_eft)
                + pow4(F0)*u04*Integrate(bind(Mterms, c,  7,  cref(P_l), k, r, _1),  YMIN, YMAX ,  error_eft)
    	//B terms
    		                              + u10*Mterms(c, 8, cref(P_l), k, r, -10.)
                                      + 3.*F0*u11*Mterms(c, 8, cref(P_l), k, r, -10.)
                                      + pow2(F0)*u11*Mterms(c, 9, cref(P_l), k, r, -10.)
                                      + pow2(F0)*u12*Mterms(c, 10, cref(P_l), k, r, -10.)
                                      + pow3(F0)*u12*Mterms(c, 11, cref(P_l), k, r, -10.));
    }


real EFT::NL_mat(real k, real bl, real U, int c, int e) const{
          int n3 = 300;
          real KMAX = QMAXp/k;
          real KMIN = QMINp/k;
          double y[n3];
          double integrand[n3];
          double u0x[8];
          double F0 = rempf_eft*fdgp_spt;
          u0x[0]=u0(k,U,F0,0.,0,e,1);
          u0x[1]=u0(k,U,F0,0.,1,e,1);
          u0x[2]=u0(k,U,F0,0.,2,e,1);
          u0x[3]=u0(k,U,F0,0.,3,e,1);
          u0x[4]=u0(k,U,F0,0.,4,e,1);
          u0x[5]=u0(k,U,F0,0.,0,e,5);
          u0x[6]=u0(k,U,F0,0.,1,e,5);
          u0x[7]=u0(k,U,F0,0.,2,e,5);
          for (int i = 0; i<n3; i++){
          y[i] = KMIN * exp(i*log(KMAX/KMIN)/(n3-1.));
          integrand[i] = midintmat(u0x,bl, c, e, cref(P_l), k, U, y[i]);
          }
        double res = 0.;
          for( int i = 1; i < n3; ++ i ){
        res += 0.5 * (y[i] - y[i-1])*(integrand[i] + integrand[i-1]);
            }
          return  k*k*k/(4*M_PI*M_PI) * res;
          }





/* GC-EFT WORKING ZONE */

/* TAYLOR EXPANSION BIAS MODEL terms: NO Wiggle version (used for resummation)*/
static real lbias_selec_decomp(const PowerSpectrum& P_l, int a, int c,  double k, double r, double x) {
  real d = 1 + r*r - 2*r*x;
  if( d < 1e-5)
      return 0.;
  else{
    double d1 = sqrt(d);
    double sker = 2./3.-(1.-pow2((x-r)/d1));
    double dker = -2./7.*(1.-x*x);
    double pld,plk,plkr;
    double G2ker = G2eds(k*r,k*d1,(x-r)/d1);
    double F2ker = F2eds(k*r,k*d1,(x-r)/d1);

    switch (c) {
      case 1:
      plk = RENW(k);
      pld = RENW(k*d1);
      plkr = RENW(k*r);
      case 2:
      plk = P_l(k);
      pld = P_l(k*d1);
      plkr = P_l(k*r);
    }

switch (a) {
  // b1 b2
  case 1:
      return  plkr*pld*pow2(r)*(2.*F2ker +(-4./7.)*sker)- pow2(r*plkr)*2.*(-4./7.)/3.;
      break;
  // b1^2
  case 2:
      return plkr*pld*pow2(r)*(-8./7.*sker*F2ker + 16./49.*pow2(sker)/2.) - pow2(r*plkr)*(16./49.*2./9.) +plk*32./315.*plkr*2.*105./16.*pow2(r)*(dker*sker+8./63.);
      break;
  // b1
  case 3:
      return plkr*pld*pow2(r)*(2.*4./7.*sker*F2ker - 2.*16/49.*pow2(sker)/2.) + pow2(r*plkr)*(32./49.*2./9.) -plk*32./315.*plkr*2.*105./16.*pow2(r)*(dker*sker+8./63.); // b3 terms
      break;
  // b2^2
  case 4:
      return plkr*pld*pow2(r)/2.- pow2(r*plkr)/2.;
      break;
 // b2
  case 5:
      return plkr*pld*pow2(r)*sker*4./7.  - pow2(r*plkr)*2.*4./7./3.;
      break;
// none
  case 6:
      return plkr*pld*pow2(r)*16./49.*pow2(sker)/2. - pow2(r*plkr)*(16./49.*2./9.);
      break;
// 2nd grouping P_dt
// both multiplied by f
// b2
   case 7:
      return -plkr*pld*pow2(r)*G2ker;
      break;
// b1
  case 8:
      return -plkr*pld*pow2(r)*((-4./7.)*G2ker*sker) - plk*32./315.*plkr*105./16.*pow2(r)*(dker*sker+8./63.);
      break;
// none
  case 9:
      return -plkr*pld*pow2(r)*G2ker*sker*4./7. -plk*32./315.*(-1.)*plkr*105./16.*pow2(r)*(dker*sker+8./63.);
      break;

}
}
}

/* Integrate them */
static real midpb2_decomp(const PowerSpectrum& P_l, int a, int c,  real k, real r) {
  real KMAX = QMAXp/k;
  real KMIN = QMINp/k;
  real YMIN = Max(-1., (1.+r*r-KMAX*KMAX)/2./r);
  real YMAX = Min(0.9999999999999, (1.+r*r-KMIN*KMIN)/2./r);
return Integrate(bind(lbias_selec_decomp, cref(P_l), a, c, k, r, _1), YMIN, YMAX ,  error_eft);
}


//a chooses the component and c chooses P_L or RENW linear spectrum
real EFT::Lag_bias_decomp(real k, int a, int c) const {
  int n3 = 300;
  real KMAX = QMAXp/k;
  real KMIN = QMINp/k;
 double y[n3];
 double integrand[n3];
 for (int i = 0; i<n3; i++){
 y[i] = KMIN * exp(i*log(KMAX/KMIN)/(n3-1.));
 integrand[i] = midpb2_decomp(cref(P_l), a, c, k, y[i]);
 }
double res = 0.;
 for( int i = 1; i < n3; ++ i ){
res += 0.5 * (y[i] - y[i-1])*(integrand[i] + integrand[i-1]);
   }
 return  k*k*k/(4.*M_PI*M_PI) * res;
}



/* Decomposition of TNS 1-LOOP RSD spectrum into u^i and bias terms */
// c chooses normal spectrum (1) or no wiggle spectrum (2)
// a chooses the term
static real select_abc_decomp(real k, const PowerSpectrum& P_l,  int a, int c, real r){
  real KMAX = QMAXp/k;
  real KMIN = QMINp/k;
  real YMIN = Max(-1., (1.+r*r-KMAX*KMAX)/2./r);
  real YMAX1 = Min(0.99999999999, (1.+r*r-KMIN*KMIN)/2./r);
  real YMAX = ATS(k,r);

switch (a) {

/* ABC terms */
//b u^2 f^2
    case 5:
      return 2.*(Integrate(bind(ABC, c,  2,  cref(P_l), k,  r, _1),  YMIN, YMAX , error_eft)
            + Integrate(bind(ABC, c,  6,  cref(P_l), k,  r, _1), YMIN, YMAX ,  error_eft)
            + Integrate(bind(ABC, c,  11,  cref(P_l), k,  r, _1), YMIN, YMAX ,  error_eft));
// b u^2 f^3
    case 6:
      return  -2.*(Integrate(bind(ABC, c,  15,  cref(P_l), k, r, _1), YMIN, YMAX , error_eft)
                  +Integrate(bind(ABC, c,  16,   cref(P_l), k,  r, _1),YMIN, YMAX , error_eft));


// b u^4 f^2
    case 7:
      return  2.*(Integrate(bind(ABC, c,  3,  cref(P_l), k,  r, _1),   YMIN, YMAX , error_eft)
                + Integrate(bind(ABC, c,  7,  cref(P_l), k,  r, _1),   YMIN, YMAX ,  error_eft)
                + Integrate(bind(ABC, c,  12,  cref(P_l), k,  r, _1),  YMIN, YMAX ,  error_eft));

// b u^4 f^3
    case 8:
      return -2.*(Integrate(bind(ABC, c,  19,  cref(P_l), k,  r, _1), YMIN, YMAX , error_eft)
                + Integrate(bind(ABC, c,  20,   cref(P_l), k,  r, _1), YMIN, YMAX ,  error_eft));

// b u^6 f^3
     case 9:
      return   -2.*(Integrate(bind(ABC, c,  22,   cref(P_l), k,  r, _1), YMIN, YMAX ,  error_eft)
                    + Integrate(bind(ABC, c,  23,   cref(P_l), k,  r, _1), YMIN, YMAX , error_eft));


// b^2 u^2 f
    case 10:
      return   2.*(Integrate(bind(ABC, c,  1,  cref(P_l), k, r, _1),  YMIN, YMAX ,  error_eft)
                  + Integrate(bind(ABC, c,  5,  cref(P_l), k,  r, _1),  YMIN, YMAX ,  error_eft)
                  + Integrate(bind(ABC, c,  10,  cref(P_l), k, r, _1),  YMIN, YMAX ,  error_eft));


// b^2 u^2 f^2
    case 11:
      return   2.*(Integrate(bind(ABC, c,  14,   cref(P_l), k, r, _1), YMIN, YMAX , error_eft)
                  +Integrate(bind(ABC, c,  40,  cref(P_l), k, r, _1), YMIN, YMAX1 , error_eft)/2.);


// b^2 u^2 f^3
    case 12:
      return    Integrate(bind(ABC, c,  41,  cref(P_l), k, r, _1), YMIN, YMAX1 , error_eft);

// b^2 u^2 f^4
    case 13:
      return    Integrate(bind(ABC, c,  42,  cref(P_l), k, r, _1), YMIN, YMAX1 , error_eft);

// b^2 u^4 f^2
    case 14:
      return   Integrate(bind(ABC, c,  43,  cref(P_l), k, r, _1), YMIN, YMAX1 , error_eft)
              +2.*Integrate(bind(ABC, c,  18,   cref(P_l), k,  r, _1),YMIN, YMAX ,  error_eft);

// b^2 u^4 f^3
    case 15:
      return  Integrate(bind(ABC, c,  44,  cref(P_l), k, r, _1), YMIN, YMAX1 , error_eft);

// b^2 u^4 f^4
    case 16:
      return  Integrate(bind(ABC, c,  45,  cref(P_l), k, r, _1), YMIN, YMAX1 , error_eft);

// b^2 u^6 f^3
    case 17:
      return Integrate(bind(ABC, c,  46,  cref(P_l), k, r, _1), YMIN, YMAX1 , error_eft);

// b^2 u^6 f^4
    case 18:
      return Integrate(bind(ABC, c,  47,  cref(P_l), k, r, _1), YMIN, YMAX1 , error_eft);

// b^2 u^8 f^4
    case 19:
      return Integrate(bind(ABC, c,  48,  cref(P_l), k, r, _1), YMIN, YMAX1 , error_eft);


// b^0 u^2 f^4
    case 20:
      return 2.*Integrate(bind(ABC, c,  17,   cref(P_l), k, r, _1), YMIN, YMAX , error_eft);


// b^0 u^4 f^3
    case 21:
      return 2.*(Integrate(bind(ABC, c,  2,  cref(P_l), k,  r, _1),  YMIN, YMAX ,  error_eft)
                   + Integrate(bind(ABC, c,  8,  cref(P_l), k,  r, _1),   YMIN, YMAX ,  error_eft)
                   + Integrate(bind(ABC, c,  11,  cref(P_l), k, r, _1), YMIN, YMAX ,  error_eft));

// b^0 u^4 f^4
    case 22:
      return 2.*Integrate(bind(ABC, c,  21,   cref(P_l), k,  r, _1),YMIN, YMAX , error_eft);

// b^0 u^6 f^3
    case 23:
      return 2.*(Integrate(bind(ABC, c,  4,  cref(P_l), k, r, _1),   YMIN, YMAX , error_eft)
                    + Integrate(bind(ABC, c,  9,  cref(P_l), k,  r, _1),  YMIN, YMAX ,  error_eft)
                    + Integrate(bind(ABC, c,  13,  cref(P_l), k,  r, _1),YMIN, YMAX ,  error_eft));

// b^0 u^6 f^4
    case 24:
      return 2.*Integrate(bind(ABC, c,  24,   cref(P_l), k,  r, _1),YMIN, YMAX ,  error_eft);

// b^0 u^8 f^4
    case 25:
      return 2.*(Integrate(bind(ABC, c,  25,   cref(P_l), k,  r, _1), YMIN, YMAX , error_eft));


}
}


real EFT::TNS_decomp(real k, int a, int c) const{
  SPT spt(C,P_l,1e-3);

if (a<5 || a>25) {
switch (a) {
      case 0:
        return   P_l(k);

      case 1:
        return   (spt.P13D_dd(k)+spt.P22D_dd(k))/pow4(D_spt/dnorm_spt);

      case 2:
        return   (spt.P13D_dt(k) + spt.P22D_dt(k))/pow4(D_spt/dnorm_spt)/(-fl_spt);

      case 3:
            return   (spt.P13D_tt(k) + spt.P22D_tt(k))/pow4(D_spt/dnorm_spt)/pow2(fl_spt);

      case 4:
        return   sigma_v/pow2(D_spt/dnorm_spt);


      case 26:
        return   RENW(k);

      case 27:
        return   (Pres13_dd(k) + Pres22_dd(k))/pow4(sig_8_eft*D_spt/dnorm_spt);

      case 28:
        return   (Pres13_dt(k) + Pres22_dt(k))/pow4(sig_8_eft*D_spt/dnorm_spt)/(-fl_spt);

      case 29:
        return   (Pres13_tt(k) + Pres22_tt(k))/pow4(sig_8_eft*D_spt/dnorm_spt)/pow2(fl_spt);

      case 30:
        return   ANW/pow2(sig_8_eft*D_spt/dnorm_spt);

}}
else{
       int n3 = 300;
       real KMAX = QMAXp/k;
       real KMIN = QMINp/k;
      double y[n3];
      double integrand[n3];
      for (int i = 0; i<n3; i++){
      y[i] = KMIN * exp(i*log(KMAX/KMIN)/(n3-1.));
      integrand[i] = select_abc_decomp(k, cref(P_l), a, c, y[i]);
      }
    double res = 0.;
      for( int i = 1; i < n3; ++ i ){
    res += 0.5 * (y[i] - y[i-1])*(integrand[i] + integrand[i-1]);
        }
      return  k*k*k/(4.*M_PI*M_PI) * res;
      }

}

/* EFT and TNS model spline construction */
// ps terms
Spline mlpdd, mlpdt, mlptt, mlpddr, mlpdtr, mlpttr;
// abc terms
Spline abc1,abc2,abc3,abc4,abc5,abc6,abc7,abc8,abc9,abc10,abc11,abc12,abc13,abc14,abc15,abc16,abc17,abc18,abc19,abc20,abc21;
// bias terms
Spline bias1,bias2,bias3,bias4,bias5,bias6,bias7,bias8,bias9;

// c decides for normal (2) or resummed (1)
real EFT::Model_init(int wig) const{
  vector<double> kval_table, mlpddt, mlpdtt, mlpttt, mlpddrt, mlpdtrt, mlpttrt;
  vector<double> abc1t,abc2t,abc3t,abc4t,abc5t,abc6t,abc7t,abc8t,abc9t,abc10t,abc11t,abc12t,abc13t,abc14t,abc15t,
                 abc16t,abc17t,abc18t,abc19t,abc20t,abc21t;
  vector<double> bias1t,bias2t,bias3t,bias4t,bias5t,bias6t,bias7t,bias8t,bias9t;
  int n1 = 150;
  double kmin = 1e-4;
  double kmax = 0.5;
  for(int i = 0; i<n1; i++ ){
    double k = kmin*exp(i*log(kmax/kmin)/(n1-1.));
    kval_table.push_back(k);
    mlpddt.push_back(TNS_decomp(k,1,wig));
    mlpdtt.push_back(TNS_decomp(k,2,wig));
    mlpttt.push_back(TNS_decomp(k,3,wig));
    mlpddrt.push_back(TNS_decomp(k,27,wig));
    mlpdtrt.push_back(TNS_decomp(k,28,wig));
    mlpttrt.push_back(TNS_decomp(k,29,wig));

    abc1t.push_back(TNS_decomp(k,5,wig));
    abc2t.push_back(TNS_decomp(k,6,wig));
    abc3t.push_back(TNS_decomp(k,7,wig));
    abc4t.push_back(TNS_decomp(k,8,wig));
    abc5t.push_back(TNS_decomp(k,9,wig));
    abc6t.push_back(TNS_decomp(k,10,wig));
    abc7t.push_back(TNS_decomp(k,11,wig));
    abc8t.push_back(TNS_decomp(k,12,wig));
    abc9t.push_back(TNS_decomp(k,13,wig));
    abc10t.push_back(TNS_decomp(k,14,wig));
    abc11t.push_back(TNS_decomp(k,15,wig));
    abc12t.push_back(TNS_decomp(k,16,wig));
    abc13t.push_back(TNS_decomp(k,17,wig));
    abc14t.push_back(TNS_decomp(k,18,wig));
    abc15t.push_back(TNS_decomp(k,19,wig));
    abc16t.push_back(TNS_decomp(k,20,wig));
    abc17t.push_back(TNS_decomp(k,21,wig));
    abc18t.push_back(TNS_decomp(k,22,wig));
    abc19t.push_back(TNS_decomp(k,23,wig));
    abc20t.push_back(TNS_decomp(k,24,wig));
    abc21t.push_back(TNS_decomp(k,25,wig));

    bias1t.push_back(Lag_bias_decomp(k,1,wig));
    bias2t.push_back(Lag_bias_decomp(k,2,wig));
    bias3t.push_back(Lag_bias_decomp(k,3,wig));
    bias4t.push_back(Lag_bias_decomp(k,4,wig));
    bias5t.push_back(Lag_bias_decomp(k,5,wig));
    bias6t.push_back(Lag_bias_decomp(k,6,wig));
    bias7t.push_back(Lag_bias_decomp(k,7,wig));
    bias8t.push_back(Lag_bias_decomp(k,8,wig));
    bias9t.push_back(Lag_bias_decomp(k,9,wig));
      }

    mlpdd = LinearSpline(kval_table,mlpddt);
    mlpdt = LinearSpline(kval_table,mlpdtt);
    mlptt = LinearSpline(kval_table,mlpttt);
    mlpddr = LinearSpline(kval_table,mlpddrt);
    mlpdtr = LinearSpline(kval_table,mlpdtrt);
    mlpttr = LinearSpline(kval_table,mlpttrt);

    abc1 = LinearSpline(kval_table,abc1t);
    abc2 = LinearSpline(kval_table,abc2t);
    abc3 = LinearSpline(kval_table,abc3t);
    abc4 = LinearSpline(kval_table,abc4t);
    abc5 = LinearSpline(kval_table,abc5t);
    abc6 = LinearSpline(kval_table,abc6t);
    abc7 = LinearSpline(kval_table,abc7t);
    abc8 = LinearSpline(kval_table,abc8t);
    abc9 = LinearSpline(kval_table,abc9t);
    abc10 = LinearSpline(kval_table,abc10t);
    abc11 = LinearSpline(kval_table,abc11t);
    abc12 = LinearSpline(kval_table,abc12t);
    abc13 = LinearSpline(kval_table,abc13t);
    abc14 = LinearSpline(kval_table,abc14t);
    abc15 = LinearSpline(kval_table,abc15t);
    abc16 = LinearSpline(kval_table,abc16t);
    abc17 = LinearSpline(kval_table,abc17t);
    abc18 = LinearSpline(kval_table,abc18t);
    abc19 = LinearSpline(kval_table,abc19t);
    abc20 = LinearSpline(kval_table,abc20t);
    abc21 = LinearSpline(kval_table,abc21t);

    bias1 = LinearSpline(kval_table,bias1t);
    bias2 = LinearSpline(kval_table,bias2t);
    bias3 = LinearSpline(kval_table,bias3t);
    bias4 = LinearSpline(kval_table,bias4t);
    bias5 = LinearSpline(kval_table,bias5t);
    bias6 = LinearSpline(kval_table,bias6t);
    bias7 = LinearSpline(kval_table,bias7t);
    bias8 = LinearSpline(kval_table,bias8t);
    bias9 = LinearSpline(kval_table,bias9t);

}

// eft model with bias a la Roy and Donald
real EFT::pkmu_eft(double k, int a, int c, double params[]) const{
double f,Dl,b,b2,n0,cs0,cs2,cs4;
  f=params[0];
  Dl=params[1];
  b=params[2];
  b2=params[3];
  n0=params[4];
  cs0=params[5];
  cs2=params[6];
  cs4=params[7];
  double d2 = pow2(Dl);
  double d4 = pow2(d2);
  double f2 = pow2(f);
  double f4 = pow2(f2);

  double b02 = pow2(b);
  double b22 = pow2(b2);

  double k2 = pow2(k);
  // double u2 = pow2(u);
  // double u4 = pow2(u2);
  // double u6 = pow2(u2)*u2;
  // double u8 = pow2(u4);

  double sv = TNS_decomp(k,4,1);
  double anw = TNS_decomp(k,30,1);

  double pddresum, pdtresum, pttresum;


  pddresum = d2*RENW(k)+d4*mlpddr(k)+ exp(-0.5*k2*d2*anw)*(d2*P_l(k)+d4*mlpdd(k)-(d2*RENW(k)+d4*mlpddr(k))+ 0.5*k2*d2*anw*(d2*(P_l(k)-RENW(k))));
  pdtresum = d2*f*RENW(k)+d4*f*mlpdtr(k)+ exp(-0.5*k2*d2*anw)*(d2*f*P_l(k)+d4*f*mlpdt(k)-(d2*f*RENW(k)+d4*f*mlpdtr(k))+ 0.5*k2*d2*anw*(d2*f*(P_l(k)-RENW(k))));
  pttresum = d2*f2*RENW(k)+d4*f2*mlpttr(k)+ exp(-0.5*k2*d2*anw)*(d2*f2*P_l(k)+d4*f2*mlptt(k)-(d2*f2*RENW(k)+d4*f2*mlpttr(k))+ 0.5*k2*d2*anw*(d2*f2*(P_l(k)-RENW(k))));

  double u0i = factL(k, d2*sv*k2, f, 1., 0 , c,  1);
  double u2i = factL(k, d2*sv*k2, f, 1., 1 , c,  1);//pow2(u);
  double u4i = factL(k, d2*sv*k2, f, 1., 2 , c,  1);//pow2(u2);

  double u2 = factL(k,  0., f, 1., 1 , c,  1);//pow2(u);
  double u4 = factL(k,  0., f, 1., 2 , c,  1);//pow2(u2);
  double u6 = factL(k,  0., f, 1., 3 , c,  1); //pow2(u2)*u2;
  double u8 = factL(k,  0., f, 1., 4 , c,  1); // pow2(u4);


  // counter term type
  // double cterm;
  //   switch (a) {
  //     case 1:
  //     cterm = (cs0*u0 + cs2*u2 + cs4*u4 + u6*(f*f2*cs0 - f2*cs2 + f*cs4));
  //     break;
  //     case 2:
  //     switch (c) {
  //       case 1:
  //       cterm =cs0;
  //       break;
  //       case 2:
  //       cterm =cs2;
  //       break;
  //       case 3:
  //       cterm =cs4;
  //       break;
  //         }
  //     break;
  //     default:
  //         warning("EFT: invalid indices, a = %d\n", a);
  //     return 0;
  //   }


  return     u0i*(b02*pddresum  +  d4*(b*b2*bias1(k) + b02*bias2(k) + b*bias3(k) + b22*bias4(k)
                                                     + b2*bias5(k) + bias6(k)) + n0)
                       + 2.*u2i*(b*pdtresum  -  d4*f*(b2*bias7(k) + b*bias8(k) + bias9(k)))
                       + u4i*pttresum

        + d4*(b*u2*f2*abc1(k) + b*u2*f*f2*abc2(k) + b*u4*f2*abc3(k) + b*u4*f*f2*abc4(k) + b*u6*f*f2*abc5(k)
            + b02*u2*f*abc6(k) + b02*u2*f2*abc7(k) + b02*u2*f*f2*abc8(k) + b02*u2*f4*abc9(k) + b02*u4*f2*abc10(k)
            + b02*u4*f*f2*abc11(k) + b02*u4*f4*abc12(k)+ b02*u6*f*f2*abc13(k) + b02*u6*f4*abc14(k) +b02*u8*f4*abc15(k)
            + u2*f4*abc16(k) + u4*f*f2*abc17(k) + u4*f4*abc18(k) + u6*f*f2*abc19(k) + u6*f4*abc20(k) + u8*f4*abc21(k))

        -2.*k2*d2*b02*RENW(k)*(cs0 + cs2*u2 + cs4*u4 + u6*(f*f2*cs0 - f2*cs2 + f*cs4));

}


// tns model with bias a la Roy and Donald
// a chooses Lorentzian (2) or Gaussian (1) form of FoG factor
// c chooses multipole (1 mono, 2 quad, 3 hex)
real EFT::pkmu_tns(double k, int a, int c,  double params[]) const{
double f,Dl,b,b2,n0,sigv;
  f=params[0];
  Dl=params[1];
  b=params[2];
  b2=params[3];
  n0=params[4];
  sigv=params[5];

  double d2 = pow2(Dl);
  double d4 = pow2(d2);
  double f2 = pow2(f);
  double f4 = pow2(f2);

  double b02 = pow2(b);
  double b22 = pow2(b2);

  double k2 = pow2(k);
  double u0 = factL(k, sigv, 1., 1., 0 , c, a +5 );
  double u2 = factL(k, sigv, 1., 1., 1 , c, a +5 );//pow2(u);
  double u4 = factL(k, sigv, 1., 1., 2 , c, a +5);//pow2(u2);
  double u6 = factL(k, sigv, 1., 1., 3 , c, a +5 ); //pow2(u2)*u2;
  double u8 = factL(k, sigv, 1., 1., 4 , c, a +5 ); // pow2(u4);

  // double FoG;
  // switch (a) {
  //   case 1:
  //   FoG = 1./(1. + (k2*u2*pow2(sigv))/2.);
  //   case 2:
  //   FoG = exp(-k2*u2*pow2(sigv));
  // }
  return  pow2(b)*u0*(d2*P_l(k)+d4*mlpdd(k))  +  d4*u0*(b*b2*bias1(k) + b02*bias2(k) + b*bias3(k) + b22*bias4(k)
                                                     + b2*bias5(k) + bias6(k)) + u0*n0
                       + 2.*u2*b*(d2*f*P_l(k)+d4*f*mlpdt(k)  -  d4*f*(b2*bias7(k) + b*bias8(k) + bias9(k)))
                       + u4*(d2*f2*P_l(k)+d4*f2*mlptt(k))
        + d4*(b*u2*f2*abc1(k) + b*u2*f*f2*abc2(k) + b*u4*f2*abc3(k) + b*u4*f*f2*abc4(k) + b*u6*f*f2*abc5(k)
            + b02*u2*f*abc6(k) + b02*u2*f2*abc7(k) + b02*u2*f*f2*abc8(k) + b02*u2*f4*abc9(k) + b02*u4*f2*abc10(k)
            + b02*u4*f*f2*abc11(k) + b02*u4*f4*abc12(k)+ b02*u6*f*f2*abc13(k) + b02*u6*f4*abc14(k) +b02*u8*f4*abc15(k)
            + u2*f4*abc16(k) + u4*f*f2*abc17(k) + u4*f4*abc18(k) + u6*f*f2*abc19(k) + u6*f4*abc20(k) + u8*f4*abc21(k));
          }


// Kaiser-pheno model
// real EFT::pkmu_kp(double k, int c,  double params[]) const{
// double f,Dl,b,c1,c2,c3;
//   f=params[0];
//   Dl=params[1];
//   b=params[2];
//   c1=params[3];
//   c2=params[4];
//   c3=params[5];
//
//   double d2 = pow2(Dl);
//   double f2 = pow2(f);
//   double k2 = pow2(k);
//
//   double u0 = factL(k,  0., 1., 1., 0 , c, 6 );
//   double u2 = factL(k,  0., 1., 1., 1 , c, 6 );//pow2(u);
//   double u4 = factL(k,  0., 1., 1., 2 , c, 6);//pow2(u2);
//   double u6 = factL(k,  0., 1., 1., 3 , c, 6 ); //pow2(u2)*u2;
//   double u8 = factL(k,  0., 1., 1., 4 , c, 6 ); // pow2(u4);
//
//
//   return  b2*u0*d2*P_l(k)  + 2.*u2*b*d2*f*P_l(k) +  u4*d2*f2*P_l(k) + c1*k2*d2*P_l(k) + c2*k2*d2*u2*P_l(k) + c3*k3 ;
//
//   -2.*k2*d2*b02*RENW(k)*(cs0 + cs2*u2 + cs4*u4 + u6*(f*f2*cs0 - f2*cs2 + f*cs4));
//
//           }
