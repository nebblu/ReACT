#if HAVE_CONFIG_H
# include <config.h>
#endif
#include "CorrelationFunction.h"
#include "PowerSpectrum.h"
#include "Quadrature.h"
#include "SpecialFunctions.h"
#include "array.h"
#include "Spline.h"
#include "SPT.h"
#include "RegPT.h"

#include <cuba.h>
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

//#include<omp.h>

// LOAD LINEAR POWER SPECTRUM ///
CorrelationFunction::CorrelationFunction(const Cosmology& C, const PowerSpectrum& P, real kmin_, real kmax_)
    :C(C), P(P), kmin(kmin_), kmax(kmax_)
{ }


// K magnitude integrals for GSM components
const real KMIN_e = 1e-4;
const real KMAX_e = 10.;

// K magnitude integrals for PS Transforms
const real KMIN_ps = 1e-4;
const real KMAX_ps = 30.;

// Errors on magnitude integrals of gsm terms (1 for analytic, 2 numeric - should be the same)
const double error1 = 1e-3;

//Error on A13 k integral
const double error3 = 1e-3;

// Error on magnitude integral of PS loop
const double error4 = 1e-3;

// Array holding F2[k,p], G1[p] and G1[k] to use in A13
//Note the 1st index needs to be matched to Aterms_init number of loops
double A13_F2A[100][n2][n1];
double A13_G1P[100][n2];
double A13_G1K[100];


/*Selection functions for numerical functions */
inline int num_selec_mag(double kmin, double kmax, double y){
		int mag_int;
	//	 mag_int = (n2*1.-1)*(y-QMINp/0.5)*0.005/QMAXp + 1./n2; // linear
	//	 mag_int = sqrt((y-QMINp/kmax)*kmin/QMAXp)*(n2-1); //quadratic
			mag_int = (int)round((n2-1)*log((y*kmax)/KMIN_e)/log(KMAX_e*kmax/(kmin*KMIN_e))); //exponential
			return mag_int;
	}


// Old stuff

void ComputeXiLM(int l, int m, const PowerSpectrum& P,
                 int Nr, const double r[], double xi[],
                 int Nk, double kmin, double kmax)
 {
    assert(Nk > 0 && (Nk % 2) == 0);
    const double dk = (kmax - kmin)/Nk;

    /* Choose appropriate spherical Bessel function */
    double (*sj)(double x);
    if(l == 0)      sj = SphericalBesselJ0;
    else if(l == 1) sj = SphericalBesselJ1;
    else if(l == 2) sj = SphericalBesselJ2;
    else if(l == 3) sj = SphericalBesselJ3;
    else if(l == 4) sj = SphericalBesselJ4;
    else {
        fprintf(stderr, "ComputeXiLM: l = %d not supported\n", l);
        return;
    }

    array k = array::linspace(kmin, kmax, Nk+1);
    array mult(Nk+1);

    #pragma omp parallel for
    for(int j = 0; j <= Nk; j++) {
        /* Multiplicative factor for Simpson's rule: either 1, 2, or 4 */
        mult[j] = 2 + 2*(j % 2) - (j == 0) - (j == Nk);
        /* All other purely k-dependent factors */
        mult[j] *= P(k[j]) * pow(k[j], m) * (dk/3) / (2*M_PI*M_PI);
    }

    /* Integrate $P(k) k^m j_l(kr) dk$ over the interval $[kmin,kmax]$ using Simpson's rule */
    #pragma omp parallel for
    for(int i = 0; i < Nr; i++) {
        xi[i] = 0;
        for(int j = 0; j <= Nk; j++)
            xi[i] += mult[j] * sj(k[j]*r[i]);
    }
}
static real f(const PowerSpectrum& P, real r, real k) {
    return k*sin(k*r)*P(k);
}
real CorrelationFunction::Evaluate(real r) const {
    return 1./(2*M_PI*M_PI*r) * Integrate(bind(f, cref(P), r, _1), kmin, kmax);
}
array CorrelationFunction::EvaluateMany(const array& r) const {
    int Nr = (int) r.size();
    array xi(Nr);
    #pragma omp parallel for
    for(int i = 0; i < Nr; i++)
        xi[i] = Evaluate(r[i]);

    return xi;
}


/* Gaussian Streaming Model - See Reid and White paper: arXiv:1105.4165v1 */
// a=1 analytic, a=2 numeric
// b = 0 (scale indep)
// b = 1 (scale dep)
// smo gives the number of samples to use in loop splining (should be >=200)
// scalef is the scale factor
// omega0 is the matter density
// mg1,2,3 are modified gravity params. (mg2 and mg3 used as sigma_v and gal bias! - need to change)

// Global initialization - see below for all functions
void CorrelationFunction::gsm_init(int a, int b, double s8, double scalef, double omega0, double smo, double mg1, double mg2, double mg3)const{
  IOW iow;
  RegPT rpt(C, P, 1e-3);
  iow.inite(scalef, omega0 ,mg1, mg2, mg3); // evolution factors and normalization
  rpt.Greg(s8*D_spt);
  loop_init(a,b,scalef,omega0,smo,mg1,mg2,mg3); // spline RegPT loop terms
  Aterms_init(a,b,scalef,omega0,mg1,mg2,mg3); // spline a terms
  sig_init(a,b); // spline velocity dispersion (perp and parallel to LOS)
}

/* CALCULATE THE FOURIER TRANSFORM OF 1-LOOP PS  and TNS multipoles*/

/*Precompute RegPT 1-loop PS for v12 and sigma_12 and TNS values */
// b = 0 (analytic or scale indep)
// b = 1 (scale dep)
// smo gives the number of samples

//RegPT splines
Spline ddloopy, dtloopy, ttloopy;
Spline ddloopy_num, dtloopy_num, ttloopy_num;

// Included the TNS monopole and quadrupole (and hex?) for FT to Xi multipoles
Spline tns0_num, tns2_num, tns4_num;
// F_1(k) and G1(k) used in v12_L and xi_real (linear)
Spline f1k_num,g1k_num;
// for multipoles p2 is the vel dispersion and p3 is the linear galaxy bias
// TODO on case 5: Can include the numerical part with a_term_init -> no need to initialize numerical kernels 80+62 times!
void CorrelationFunction::loop_init(int a, int b, double scalef, double omega0, double smo, double p1, double p2, double p3)const{
  RegPT rpt(C, P, 1e-3);
  IOW iow;
  int loop_N = (int)smo; // Chosen as min for convergence of xi(r) - similarly for numerical kernels with n2 = 80
  double loopy1,loopy2,loopy3,k,sd1,sd2,sd3,sd4;
  double ktns[loop_N],tns0[loop_N],tns2[loop_N],tns4[loop_N];
  vector<double> kval_table,loop_table1,loop_table2,loop_table3;
  vector<double> tns0_table,tns2_table,tns4_table;
  vector<double> f1nk_table,g1nk_table;
 switch (a) {
    case 1:
        iow.inite(scalef, omega0, p1,  p2, p3);
        rpt.sigmad_init(); // initialize sigma_d damping term
        for (int i = 0; i<loop_N ; i++) {
        k = KMIN_ps*exp(i*log(KMAX_ps/KMIN_ps)/(loop_N-1.));
        if (k>1.) {
          loopy1 = 0.;
          loopy2 = 0.;
          loopy3 = 0.;
         }
        else{
              loopy1 =  rpt.PLOOPr(1,  k);;
              loopy2 =  rpt.PLOOPr(2,  k);;
              loopy3 =  rpt.PLOOPr(3,  k);;
            }
	            kval_table.push_back(k);
              loop_table1.push_back(loopy1);
              loop_table2.push_back(loopy2);
              loop_table3.push_back(loopy3);
            }
    ddloopy = LinearSpline(kval_table,loop_table1);
    dtloopy = LinearSpline(kval_table,loop_table2);
    ttloopy = LinearSpline(kval_table,loop_table3);
       break;
  /*Numerical PS initialization  */
  /*We also initialize the linear growth */
    case 2:
      sd1 = KMIN_ps;
      sd2 = 2.;
      iow.inite(scalef, omega0, p1,  p2, p3);
      rpt.sigmad_init(); // initialize sigma_d damping term
      if (b==0 && a==2) {
        iow.initn(scalef, QMINp/sd2 ,QMAXp/sd1, 1., omega0, p1 , p2, p3);
        }
      for (int i = 0; i<loop_N ; i++) {
       k = KMIN_ps*exp(i*log(KMAX_ps/KMIN_ps)/(loop_N-1.));
       if (k>1.) {
         loopy1 = 0.;
         loopy2 = 0.;
         loopy3 = 0.;
        }
       else{
       if (b==1 && a==2) {
       iow.initn(scalef, QMINp/k ,QMAXp/k, k, omega0, p1 , p2, p3); //Initialise kernels for k
       sd1 = k;
       sd2 = k;
         }
          loopy1 = rpt.PLOOPnr(sd1, sd2, 1,  k); //delta delta
          loopy2 = rpt.PLOOPnr(sd1, sd2, 2,  k); //delta -theta
          loopy3 = rpt.PLOOPnr(sd1, sd2, 3,  k); // theta-theta
	         sd3 = F1_nk;
	         sd4 = G1_nk;
	       }
         kval_table.push_back(k);
         loop_table1.push_back(loopy1);
         loop_table2.push_back(loopy2);
         loop_table3.push_back(loopy3);
         f1nk_table.push_back(sd3);
         g1nk_table.push_back(sd4);
        }
        ddloopy_num = LinearSpline(kval_table,loop_table1);
        dtloopy_num = LinearSpline(kval_table,loop_table2);
        ttloopy_num = LinearSpline(kval_table,loop_table3);
        f1k_num = LinearSpline(kval_table,f1nk_table);
        g1k_num = LinearSpline(kval_table,g1nk_table);
        break;
  /*Analytic TNS P_i initialization  */
  // RegPT initialization of TNS multipoles (analytic)
    case 3:
  //  #pragma omp parallel for schedule(dynamic)
    iow.inite(scalef, omega0, p1,  p2, p3);
    rpt.sigmad_init(); // initialize sigma_d damping term
      for (int i = 0; i<loop_N ; i++) {
       ktns[i] = KMIN_ps*exp(i*log(KMAX_ps/KMIN_ps)/(loop_N-1.));
       if (ktns[i]>1.) {
         tns0[i] = 0.;
         tns2[i] = 0.;
         tns4[i] = 0.;
        }
       else{
       tns0[i] = rpt.PTNSMnDGPr(ktns[i], p3, p2, 1);
       tns2[i] = rpt.PTNSMnDGPr(ktns[i], p3, p2, 2);
       tns4[i] = rpt.PTNSMnDGPr(ktns[i], p3, p2, 3);
       }
	}
      for (int i = 0; i<loop_N ; i++) {
       kval_table.push_back(ktns[i]);
       tns0_table.push_back(tns0[i]);
       tns2_table.push_back(tns2[i]);
       tns4_table.push_back(tns4[i]);
        }

      tns0_num = LinearSpline(kval_table,tns0_table);
      tns2_num = LinearSpline(kval_table,tns2_table);
      tns4_num = LinearSpline(kval_table,tns4_table);
      break;
  // RegPT initialization of TNS multipoles (numeric)
      case 4:
            sd1 = KMIN_ps;
            sd2 = 2.;
            iow.inite(scalef, omega0, p1,  p2, p3);
            rpt.sigmad_init(); // initialize sigma_d damping term
            if (b==0 && a==4) {
              iow.initn(scalef, QMINp/sd2 ,QMAXp/sd1, 1., omega0, p1 , p2, p3);
              }
         //  #pragma omp parallel for schedule(dynamic)
	          for (int i = 0; i<loop_N ; i++) {
             ktns[i] = KMIN_ps*exp(i*log(KMAX_ps/KMIN_ps)/(loop_N-1.));
             if (ktns[i]>=1.) {
               tns0[i] = 0.;
               tns2[i] = 0.;
               tns4[i] = 0.;
              }
             else{
             if (b==1 && a==4) {
             iow.initn(scalef, QMINp/ktns[i] ,QMAXp/ktns[i], ktns[i], omega0, p1 , p2, p3); //Initialise kernels for k
             sd1 = ktns[i];
             sd2 = ktns[i];
               }
                 tns0[i] = rpt.PTNSMmgr( ktns[i], sd1, sd2, 1, p3, p2);
                 tns2[i] = rpt.PTNSMmgr( ktns[i], sd1, sd2, 2, p3, p2);
                 tns4[i] = rpt.PTNSMmgr( ktns[i], sd1, sd2, 3, p3, p2);
                }
		            }
         for (int i = 0; i<loop_N ; i++) {
		             kval_table.push_back(ktns[i]);
                 tns0_table.push_back(tns0[i]);
                 tns2_table.push_back(tns2[i]);
                 tns4_table.push_back(tns4[i]);
               }
                tns0_num = LinearSpline(kval_table,tns0_table);
                tns2_num = LinearSpline(kval_table,tns2_table);
                tns4_num = LinearSpline(kval_table,tns4_table);
            break;
          // 1st order kernels initialization : b = 0 for scale indep (initialize in gsm2.cpp), b=1 scale dep.
          // Use if we want to output xi_real linear by itself
          case 5:
          if (b==0) {
            iow.initn(scalef, KMIN_ps/KMAX_ps ,KMAX_ps/KMIN_ps, 1., omega0, p1 , p2, p3);
            }
          for (int i = 0; i<loop_N ; i++) {
           k = KMIN_ps*exp(i*log(KMAX_ps/KMIN_ps)/(loop_N-1.));
           if (b==1) {
           iow.initn(scalef, KMIN_ps , KMAX_ps, k, omega0, p1 , p2, p3); //Initialise kernels for k
            }
           sd1 = F1_nk;
           sd2 = G1_nk;
           kval_table.push_back(k);
           f1nk_table.push_back(sd1);
           g1nk_table.push_back(sd2);
          }
          f1k_num = LinearSpline(kval_table,f1nk_table);
          g1k_num = LinearSpline(kval_table,g1nk_table);
          break;
}
}
//Correlation function : Linear, 1-loop and FT of TNS multipoles
//a =1 : Linear(GR)
//a =2 : Regpt 1-loop(GR)
//a =3 : linear(MG)
//a =4 : Regpt 1-loop(MG)
//a =5 : TNS_0
//a =6 : TNS_2
//a =7 : TNS_4
static real xi_integrand(const PowerSpectrum& P, real r, int a, real k) {
	switch (a) {

//real space
		case 1:
			return  SphericalBesselJ0(k*r)*pow2(k)*pow2(D_spt/dnorm_spt)*P(k);
      break;
  	case 2:
      return SphericalBesselJ0(k*r)*pow2(k)*ddloopy(k);
      break;
    case 3:
      return  SphericalBesselJ0(k*r)*pow2(k)*pow2(f1k_num(k)/dnorm_spt)*P(k);
      break;
    case 4:
      return  SphericalBesselJ0(k*r)*pow2(k)*ddloopy_num(k);
      break;
// multipoles
    case 5:
      return  SphericalBesselJ0(k*r)*pow2(k)*tns0_num(k);
      break;
    case 6:
      return -SphericalBesselJ2(k*r)*pow2(k)*tns2_num(k);
      break;
    case 7:
      return SphericalBesselJ4(k*r)*pow2(k)*tns4_num(k);
      break;
      default:
    throw "Invalid index: Choose 1-7";
}
}

real CorrelationFunction::xi_real(double r, int a) const {
   return 1./(2.*pow2(M_PI))*Integrate(bind(xi_integrand, cref(P), r, a, _1), KMIN_ps, KMAX_ps , error4);
}

/*GAUSSIAN STREAMING MODEL COMPONENTS AS IN BETH/REID 2011 - ANALYTIC -GR */

// Note in general a = 1 is GR (analytic) and a = 2 is general MG (numerical)
// Linear infall velocity and velocity dispersion (Eqs. 7 and 11)

//Linear infall velocity :
static double v12_L_integrand(int a, double r, const PowerSpectrum& P, double k){
  switch (a) {
    case 1:
	     return -fl_spt*pow2(Dl_spt/dnorm_spt)*P(k)*SphericalBesselJ1(k*r)*k;
       break;
    case 2:
       return g1k_num(k)*f1k_num(k)*pow2(1./dnorm_spt)*P(k)*SphericalBesselJ1(k*r)*k ;
       break;
}
}
real CorrelationFunction::v12_L(int a, double b, double r) const{
	return b/pow2(M_PI)*Integrate(bind(v12_L_integrand, a, r, cref(P), _1), KMIN_e , KMAX_e, error1);
}


// APPENDIX TERMS:
//Infall velocity terms
//A5
static double a5_kernel(const PowerSpectrum& P, double k, double y){
  double ysqr = pow2(y);
  double ysqrinv;
  double atanhdiff;
  double myresult;
  if(y < 0.5) {
    myresult = -1./3.+4./7.*ysqr-12./35.*pow2(ysqr)+1./420.*(27.*pow3(ysqr)-12.*pow4(ysqr)+9.*pow5(ysqr));
  if(y > 0.001) { //to avoid problems with 1/y
    atanhdiff = atanh(2.*y/(1.+ysqr)) - 2.*y*(1.+ysqr/3.+pow2(ysqr)/5.);
    myresult += 9./168.*pow3(-1.+ysqr)*(atanhdiff)/y;
    }
    }
  else {
    if(y > 2.0) {
      ysqrinv = 1./ysqr;
      myresult = 1./105. - 12./245.*ysqrinv -17./980.*pow2(ysqrinv) + 6./245.*pow3(ysqrinv) - 3./196.*pow4(ysqrinv);
      atanhdiff = atanh(2.*y/(1.+ysqr)) - 2./y*(1.+ysqrinv/3.+pow2(ysqrinv)/5.+pow3(ysqrinv)/7.);
      if(y<30.) { //after 30, it becomes a 1.5e-6 correction; this is causing us a convergence problem anyway.
        myresult += 9./168.*pow3(-1.+ysqr)*(atanhdiff)/y;
        }
      }
    else {
      //just evaluate the function.
      myresult = 2.*y*(9.+52.*ysqr-9.*pow2(ysqr)) - 56.*y*(1.+ysqr);
      if(fabs(y-1) > 0.001) { //stay away from singularity at y=1; evaluates to 0 anyway for this term.
        myresult += 9.*pow3(-1.+ysqr)*atanh(2.*y/(1.+ysqr));
        }
      myresult = myresult/168./y;
      }
    }
    return myresult*P(k)*P(k*y);
  }
//A6
static double a6_integrand(const PowerSpectrum& P, double k, double y, double x){
  double ysqr,xsqr,d;
  ysqr = pow2(y);
  xsqr = pow2(x);
  d = (1.+ysqr-2.*y*x);
  if (d<1e-5){
    return 0.;
  }
  if(x == 1.) {  //at x=1, diverges!!
  return (0.5/(1.-y))*P(k*sqrt(d))*P(k*y);
    }
  return ((3.*y*x-10.*y*x*xsqr+7.*xsqr)/(14.*d))*P(k*sqrt(d))*P(k*y);
}
// Perform angular integration
static double a6_kernel(const PowerSpectrum& P,  double k, double y){
	return Integrate(bind(a6_integrand ,cref(P), k ,y, _1), XMIN , XMAX, error1);
}
//Total kernel and include evolution
static double a56_kernel(const PowerSpectrum& P,  double k, double y){
 return -fl_spt*pow4(k*Dl_spt/dnorm_spt)*(a5_kernel(cref(P),k,y) + a6_kernel(cref(P),k,y));
}
//Will include bias in v_in expression
static double a56_integrand(const PowerSpectrum& P,  double k){
 double YMAX=KMAX_e/k;
 double YMIN=KMIN_e/k;
	return Integrate(bind(a56_kernel,cref(P),k, _1), YMIN , YMAX , error1)/(2.*pow4(M_PI));
}

// Numerical a56
static double a56_num_kernel(int y1, const PowerSpectrum& P, double k, double y){
  double myresult=0.;
  double temp_a56;
  double check;
  double ker[6];
for( int i = 0; i < n1; i++ )
      {
      double d = 1.+ y*y - 2.*y*x128[i];
      check =1.;
      if(d < 1e-5){
      check =0.;
      }
      else {
        // 1st order kernels F1/G1(k-p)
         ker[0] = F1kmp_nk[i*n2 + y1];
        // 1st order kernels F1/G1(p)
         ker[1] = F1p_nk[i*n2 + y1];
         ker[2] = G1p_nk[i*n2 + y1];
        //symmetrized 2nd order kernels for ps F2/G2(p,k-p)
         ker[3] = F2_nk[i*n2 + y1];
        //symmetrized 2nd order kernels F2/G2(-p,k)
         ker[4] = F2B_nk[i*n2 + y1];
         ker[5] = G2B_nk[i*n2 + y1];

      temp_a56 =y*(x128[i]*ker[0]*ker[2]*ker[3]*P(k*sqrt(d))*check*P(k*y)+(F1_nk*ker[1]*ker[5]*y*(1.-y*x128[i])/d +ker[2]*F1_nk*ker[4]*x128[i])*P(k)*P(k*y));
      }
      myresult += w128[i] * temp_a56;
      }
      return myresult;
}

static double a56_integrand_num(const PowerSpectrum& P, double k1, double k2, double k) {
    	real KMAX = KMAX_e/k;
    	real KMIN = KMIN_e/k;
      int y1 = num_selec_mag(k1, k2, KMIN);
      int y2 = num_selec_mag(k1, k2, KMAX);
      double y[n2];
      double integrand[n2];
      for (int i = y1; i<=y2; i++){
      y[i] = KMIN_e/k2 * exp(i*log(KMAX_e*k2/(KMIN_e*k1))/(n2*1.-1.));
      integrand[i] = a56_num_kernel(i, cref(P), k, y[i]);
}
double res = 0.;
  for( int i = y1+1; i <= y2; ++ i ){
res += 0.5 * (y[i] - y[i-1])*(integrand[i] + integrand[i-1]);
}
return  1./(2.*pow4(dnorm_spt*M_PI))*pow4(k)*res;
}

//Velocity dispersion terms
//A10 (contributes a constant)
static double a10_kernel(const PowerSpectrum& P, double k, double y){
  double myresult;
  double ysqr = pow2(y);
  myresult = 2.*y*(3.+8.*y*y-3.*pow4(y));
          if(fabs(y-1) > 0.001) { //stay away from singularity at y=1; evaluates to 0 anyway for this term.
            myresult += 3.*pow3(-1.+ysqr)*atanh(2.*y/(1.+ysqr));
            }
  return  P(k*y)*P(k)*pow3(k)*myresult/168./y;
}
static double a10_integrand(const PowerSpectrum& P, double k){
  double YMIN = KMIN_e/k;
  double YMAX = KMAX_e/k;
  return  pow2(fl_spt)*pow4(Dl_spt/dnorm_spt)*Integrate(bind(a10_kernel,cref(P),k,_1),YMIN,YMAX,error1);
}
//2 factor from Eq.30 already accounted for
//ONLY NEED TO COMPUTE THIS ONCE - ex. INCLUDE IN GSM2.CPP IN EXAMPLES FILE BEFORE LOOP BEGINGS
static double a10(const PowerSpectrum& P){
  return 1./pow4(M_PI)*Integrate(bind(a10_integrand,cref(P),_1),KMIN_e,KMAX_e,error1);
}


static double a10_num_kernel(int y1 , int k_int, const PowerSpectrum& P, double k, double y){
  double myresult=0.;
  double temp_a10;
  double ker[4];
for( int i = 0; i < n1; i++ )
      {
      // Load array for F2(p,k),G1(p) and G1(k) for A13 computation later
        ker[0] = G1p_nk[i*n2 + y1];
        ker[2] = F2A_nk[i*n2 + y1];
        A13_F2A[k_int][y1][i] = ker[2];
        A13_G1P[k_int][y1] = ker[0];
        A13_G1K[k_int] = G1_nk;
      double d = 1+ y*y + 2*y*x128[i];
      if(d < 1e-5){
      temp_a10 = 0.;
      }
      else {
        // 1st order kernels F1(p)
         ker[1] = F1p_nk[i*n2 + y1];
        //symmetrized 2nd order kernels F2/G2(p,k)
         ker[3] = G2A_nk[i*n2 + y1];
      temp_a10 = G1_nk*(2.*ker[3]*ker[1]*(1.+y*x128[i])/d-x128[i]*ker[0]*ker[2]/y);
      }
      myresult += w128[i] * temp_a10;
      }
      return P(k*y)*P(k)*pow3(k)*pow2(y)*myresult;
}

static double a10_integrand_num(const PowerSpectrum& P, int k_int, double k1, double k2, double k){
  real KMAX = KMAX_e/k;
  real KMIN = KMIN_e/k;
  int y1 = num_selec_mag(k1, k2, KMIN);
  int y2 = num_selec_mag(k1, k2, KMAX);
  double y[n2];
  double integrand[n2];
  for (int i = y1; i<=y2; i++){
    y[i] = KMIN_e/k2 * exp(i*log(KMAX_e*k2/(KMIN_e*k1))/(n2*1.-1.));
    integrand[i] = a10_num_kernel(i, k_int, cref(P), k, y[i]);
  }
  double res = 0.;
    for( int i = y1+1; i <= y2; ++ i ){
  res += 0.5 * (y[i] - y[i-1])*(integrand[i] + integrand[i-1]);
  }
  return  1./(6.*pow4(dnorm_spt*M_PI))*res;
}


// A12,A13,A14,A15 all need a 2b factor (Eq.30) - we will account for 2 factor here and add in bias later.

// A12 - Note typo in paper with H^2f^2. Also we still need prefactor of 1/2pi^4 (2 will cancel with factor mentioned above)
static double a12_kernel(int a, const PowerSpectrum& P, double k, double y){
  double prefac = pow2(fl_spt)*pow4(Dl_spt/dnorm_spt)*P(k*y)*P(k)*pow3(k);
  double ysqr = pow2(y);
  double ysqrinv;
  double atanhdiff;
  double myresult;
  switch (a) {
    case 1:
    if(y < 0.5) {
      myresult = 8./35.*ysqr-24./245.*pow2(ysqr)+8./735.*(pow3(ysqr)+pow4(ysqr)/11.+3./143.*pow5(ysqr));
      return prefac*myresult;
      }
    else {
      if(y > 2.0) {
        ysqrinv = 1./ysqr;
        myresult = 8./35.-24./245.*ysqrinv+8./735.*(pow2(ysqrinv)+pow3(ysqrinv)/11.+3./143.*pow4(ysqrinv)+1./143.*pow5(ysqrinv));
        return prefac*myresult;
        }
      else {
        //just evaluate the function.
        myresult = 2.*y*(-9.+33.*ysqr+33.*pow2(ysqr)-9.*pow3(ysqr));
        if(fabs(y-1) > 0.001) { //stay away from singularity at y=1; evaluates to 0 anyway for this term.
          myresult += 9.*pow4(-1.+ysqr)*atanh(2.*y/(1.+ysqr));
          }
        myresult = myresult/(672.*ysqr*y);
        return prefac*myresult;
        }
      }
      break;
  case 2:
    if(y < 0.5) {
      myresult = -1./3.+12./35.*ysqr-12./49.*pow2(ysqr)+4./105.*pow3(ysqr)+12./2695.*pow4(ysqr)+4./3185.*pow5(ysqr);  //agrees with true result up to 2e-7
      return prefac*myresult;
      }
    else {
      if(y > 2.0) {
        ysqrinv = 1./ysqr;
        myresult = -23./105. + 12./245.*ysqrinv -4./245.*pow2(ysqrinv) - 4./1617.*pow3(ysqrinv) - 4./5005.*pow4(ysqrinv)-12./35035.*pow5(ysqrinv);  //agrees at 5e-8 level
        return prefac*myresult;
        }
      else {
        //just evaluate the function.
        myresult = 2.*y*(9.-109.*ysqr+63.*pow2(ysqr)-27.*pow3(ysqr));
        if(fabs(y-1) > 0.001) { //stay away from singularity at y=1; evaluates to 0 anyway for this term.  leaving this out is ~1e-5 difference at 0.001
          myresult += 9.*pow3(-1.+ysqr)*(1.+3.*ysqr)*atanh(2.*y/(1.+ysqr));
          }
        myresult = myresult/(672.*ysqr*y);
        return prefac*myresult;
        break;
        }
      }
    }
}

// factor of 2 from Eq.31 already accounted for in all the following a terms
static double a12_integrand(int a ,const PowerSpectrum& P, double k){
  double YMIN = KMIN_e/k;
  double YMAX = KMAX_e/k;
    return  1./pow4(M_PI)*Integrate(bind(a12_kernel, a, cref(P), k, _1), YMIN, YMAX, error1);
  }

// Numerical
  static double a12_num_kernel(int a, int y1 ,const PowerSpectrum& P, double k, double y){
    double myresult=0.;
    double temp_a12,term_selec;
    double ker[2];
  for( int i = 0; i < n1; i++ )
        {
        double d = 1+ y*y + 2*y*x128[i];
        if(d < 1e-5){
        temp_a12 = 0.;
        }
        else {
          // 1st order kernels G1(p)
           ker[0] = G1p_nk[i*n2 + y1];
          //symmetrized 2nd order kernels G2(p,k)
           ker[1] = G2A_nk[i*n2 + y1];
        switch (a) {
          case 1:
          term_selec = y*(1-pow2(x128[i]));
          break;
          case 2:
          term_selec = -(3.*pow2(x128[i])*y-y+2.*x128[i]);
          break;
        }
        temp_a12 = F1_nk*ker[0]*ker[1]*term_selec/d;
        }
        myresult += w128[i] * temp_a12;
        }
        return P(k*y)*P(k)*pow3(k)*y*myresult;
  }


static double a12_integrand_num(int a, const PowerSpectrum& P, double k1, double k2, double k){
    real KMAX = KMAX_e/k;
    real KMIN = KMIN_e/k;
    int y1 = num_selec_mag(k1, k2, KMIN);
    int y2 = num_selec_mag(k1, k2, KMAX);
    double y[n2];
    double integrand[n2];
    for (int i = y1; i<=y2; i++){
      y[i] = KMIN_e/k2 * exp(i*log(KMAX_e*k2/(KMIN_e*k1))/(n2*1.-1.));
      integrand[i] = a12_num_kernel(a,i, cref(P), k, y[i]);
    }
    double res = 0.;
      for( int i = y1+1; i <= y2; ++ i ){
    res += 0.5 * (y[i] - y[i-1])*(integrand[i] + integrand[i-1]);
    }
    return  1./(2.*pow4(dnorm_spt*M_PI))*res;
  }


// A13 only splined over r - treated later


//A14
static double a14_kernel(const PowerSpectrum& P, double k, double y, double x){
    double d,ysqr,xsqr;
    ysqr = pow2(y);
    xsqr = pow2(x);
    d = (1.+ysqr-2.*y*x);
    double prefac = pow2(fl_spt)*pow4(Dl_spt/dnorm_spt)* P(k*y)* P(k*sqrt(d))*pow3(k);

    if (d<1e-5){
      return 0.;
    }
    if(y == 1. && x == 1.) {
      return(0.);
      }
    if(x == 1.) {  //at x=1, diverges!!
      return prefac*(7./(y-1.));
      }
      return prefac*((y*x + 6.*y*x*xsqr - 7.*xsqr)/(14.*d));
  }
static double a14_integrand(const PowerSpectrum& P, double k){
    real c[2] = {KMIN_e/k, XMIN};
    real d[2] = {KMAX_e/k, XMAX};
    return 1./pow4(M_PI)*Integrate<2>(bind(a14_kernel,cref(P),k,_1,_2),c,d,error1,error1);
    }

// Numerical

static double a14_num_kernel(int y1, const PowerSpectrum& P, double k, double y){
    double myresult=0.;
    double temp_a14;
    double ker[3];
  for( int i = 0; i < n1; i++ )
        {
        double d = 1+ y*y - 2*y*x128[i];
        if(d < 1e-5){
        temp_a14 = 0.;
        }
        else {
          // 1st order kernels F1(k-p)
           ker[0] = F1kmp_nk[i*n2 + y1];
           // 1st order kernels G1(p)
           ker[1] =  G1p_nk[i*n2 + y1];
          //symmetrized 2nd order kernels for ps G2(p,k-p)
           ker[2] = G2_nk[i*n2 + y1];
        temp_a14 =  P(k*sqrt(d))*x128[i]*ker[1]*ker[0]*ker[2];
        }
        myresult += w128[i] * temp_a14;
        }
        return -P(k*y)*pow3(k)*y*myresult;
}


static double a14_integrand_num(const PowerSpectrum& P, double k1, double k2, double k){
  real KMAX = KMAX_e/k;
  real KMIN = KMIN_e/k;
  int y1 = num_selec_mag(k1, k2, KMIN);
  int y2 = num_selec_mag(k1, k2, KMAX);
  double y[n2];
  double integrand[n2];
  for (int i = y1; i<=y2; i++){
    y[i] = KMIN_e/k2 * exp(i*log(KMAX_e*k2/(KMIN_e*k1))/(n2*1.-1.));
    integrand[i] = a14_num_kernel(i, cref(P), k, y[i]);
  }
  double res = 0.;
    for( int i = y1+1; i <= y2; ++ i ){
  res += 0.5 * (y[i] - y[i-1])*(integrand[i] + integrand[i-1]);
  }
  return res/pow4(dnorm_spt*M_PI);
}


//A15
static double a15_kernel(const PowerSpectrum& P, double k, double y){
    double prefac = -pow2(fl_spt)*pow4(Dl_spt/dnorm_spt)* P(k*y)*P(k)*pow3(k);
    double ysqr = pow2(y);
    double ysqrinv;
    double atanhdiff;
    double myresult;
    if(y < 0.5) {
      myresult = -1./3. + 4./7.*ysqr-12./35.*pow2(ysqr)+12./245.*pow3(ysqr)+4./735.*pow4(ysqr)+4./2695.*pow5(ysqr);
      return prefac*myresult;
      }
    else {
      if(y > 2.0) {
        ysqrinv = 1./ysqr;
        myresult = 1./105.-12./245.*ysqrinv-4./735.*pow2(ysqrinv)-4./2695.*pow3(ysqrinv)-4./7007.*pow4(ysqrinv)-4./15015.*pow5(ysqrinv);
        return prefac*myresult;
        }
      else {
        //just evaluate the function.
        myresult = -2.*y*(19.-24.*ysqr+9.*pow2(ysqr));
        if(fabs(y-1) > 0.001) { //stay away from singularity at y=1; evaluates to 0 anyway for this term.
          myresult += 9.*pow3(-1.+ysqr)*atanh(2.*y/(1.+ysqr));
          }
        myresult = myresult/(168.*y);
        return prefac*myresult;
        }
    }
  }
static double a15_integrand(const PowerSpectrum& P, double k){
      double YMIN = KMIN_e/k;
      double YMAX = KMAX_e/k;
      return  1./pow4(M_PI)*Integrate(bind(a15_kernel, cref(P), k, _1), YMIN, YMAX, error1);
         }


static double a15_num_kernel(int y1, const PowerSpectrum& P, double k, double y){
      double myresult=0.;
      double temp_a15;
      double ker[4];
      for( int i = 0; i < n1; i++ )
         {
          double d = 1+ y*y + 2*y*x128[i];
          if(d < 1e-5){
          temp_a15 = 0.;
           }
          else {
          // 1st order kernels G1(p)/F1(p)
          ker[0] = G1p_nk[i*n2 + y1];
          ker[1] = F1p_nk[i*n2 + y1];
          //symmetrized 2nd order kernels for ps G2/F2(p,k)
          ker[2] = F2A_nk[i*n2 + y1];
          ker[3] = G2A_nk[i*n2 + y1];

          temp_a15 = G1_nk*(ker[3]*ker[1]*y*(1.+y*x128[i])/d-x128[i]*ker[0]*ker[2]);
        }
        myresult += w128[i] * temp_a15;
      }
      return -P(k*y)*P(k)*pow3(k)*y*myresult;
    }

    static double a15_integrand_num(const PowerSpectrum& P, double k1, double k2, double k){
      real KMAX = KMAX_e/k;
      real KMIN = KMIN_e/k;
      int y1 = num_selec_mag(k1, k2, KMIN);
      int y2 = num_selec_mag(k1, k2, KMAX);
      double y[n2];
      double integrand[n2];
      for (int i = y1; i<=y2; i++){
        y[i] = KMIN_e/k2 * exp(i*log(KMAX_e*k2/(KMIN_e*k1))/(n2*1.-1.));
        integrand[i] = a15_num_kernel(i, cref(P), k, y[i]);
      }
      double res = 0.;
        for( int i = y1+1; i <= y2; ++ i ){
      res += 0.5 * (y[i] - y[i-1])*(integrand[i] + integrand[i-1]);
      }
      return res/pow4(dnorm_spt*M_PI);
    }


// Splining of appendix terms as functions of k
Spline a56,a12a,a12b,a14,a15; // analytic
Spline a56n,a10n,a12an,a12bn,a14n,a15n; // Numerical

void CorrelationFunction::Aterms_init(int a, int b, double scalef, double omega0, double p1, double p2, double p3)const{
  // scale dependent kernels take <7 mins to initialize for loop_N=100
  // Note a56n is < 7% match with analytical but <1% individually (a5 and a6) because a56 is a difference of terms
  const int loop_N = 100; // Chosen as min deviation from converged value of Integrate[ai(k),KMIN_e,KMAX_e] to within 1% - similarly for numerical kernels with n2 = 80
  double k,loopy1,loopy2,loopy3,loopy4,loopy5,loopy6, k1,k2;
  double lim1 = KMIN_e;
  double lim2 = KMAX_e;
  IOW iow;
  vector<double> kval_table1;
  vector<double> a56_table,a12a_table,a12b_table,a14_table,a15_table;
  vector<double> a56_tablen,a10_tablen,a12a_tablen,a12b_tablen,a14_tablen,a15_tablen;
  // Analytic
  switch (a) {
    case 1:
    for (int i = 0; i<loop_N ; i++) {
     k1 = lim1*exp(i*log(lim2/lim1)/(loop_N*1.-1.));

     loopy1 = a56_integrand(cref(P),k1);
     loopy2 = a12_integrand(1 ,cref(P),k1);
     loopy3 = a12_integrand(2,cref(P),k1);
     loopy4 = a14_integrand(cref(P),k1);
     loopy5 = a15_integrand(cref(P),k1);

      kval_table1.push_back(k1);
     a56_table.push_back(loopy1);
     a12a_table.push_back(loopy2);
     a12b_table.push_back(loopy3);
     a14_table.push_back(loopy4);
     a15_table.push_back(loopy5);
    }

    a56 = LinearSpline(kval_table1,a56_table);
    a12a = LinearSpline(kval_table1,a12a_table);
    a12b = LinearSpline(kval_table1,a12b_table);
    a14 = LinearSpline(kval_table1,a14_table);
    a15 = LinearSpline(kval_table1,a15_table);

    break;
// Numerical terms
    case 2:
    if (b==0 && a==2) {
      iow.initn(scalef, KMIN_e/KMAX_e ,KMAX_e/KMIN_e, 1., omega0, p1 , p2, p3);
      }
    for (int i = 0; i<loop_N ; i++) {
    k = lim1*exp(i*log(lim2/lim1)/(loop_N-1.));
    k1 = lim1;
    k2 = lim2;
    if (b==1 && a==2) {
    iow.initn(scalef, KMIN_e/k ,KMAX_e/k, k, omega0, p1 , p2, p3); //Initialise kernels for k
    k1 = k;
    k2 = k;
      }
    loopy1 = a56_integrand_num(cref(P),k1,k2,k);
    loopy2 = a10_integrand_num(cref(P),i,k1,k2,k);
    loopy3 = a12_integrand_num(1,cref(P),k1,k2,k);
    loopy4 = a12_integrand_num(2,cref(P),k1,k2,k);
    loopy5 = a14_integrand_num(cref(P),k1,k2,k);
    loopy6 = a15_integrand_num(cref(P),k1,k2,k);

    kval_table1.push_back(k);
    a56_tablen.push_back(loopy1);
    a10_tablen.push_back(loopy2);
    a12a_tablen.push_back(loopy3);
    a12b_tablen.push_back(loopy4);
    a14_tablen.push_back(loopy5);
    a15_tablen.push_back(loopy6);
    }
    a56n = LinearSpline(kval_table1, a56_tablen);
    a10n = LinearSpline(kval_table1, a10_tablen);
    a12an = LinearSpline(kval_table1, a12a_tablen);
    a12bn = LinearSpline(kval_table1, a12b_tablen);
    a14n = LinearSpline(kval_table1, a14_tablen);
    a15n = LinearSpline(kval_table1, a15_tablen);
    break;
}
}
//Numerical A10 integral
static double a10_splined_kernel(double k){
  return a10n(k);
}
static double a10_num(){
  return Integrate(bind(a10_splined_kernel,_1), KMIN_e, KMAX_e, error1);
}


//A13 terms
static double a13_terms(int a, double r, const PowerSpectrum& P, double k){
 double prefac = 1./4.*fl_spt*pow2(Dl_spt/dnorm_spt)* P(k);
  switch (a) {
    case 0 :
    return prefac*I11(k*r)*k;
    break;
    case 1 :
    return prefac*I12(k*r)*k*k;
    break;
    case 2 :
    return prefac*I12(k*r);
    break;
    case 3 :
    return prefac*I30(k*r)*k*k;
    break;
    case 4 :
    return prefac*I30(k*r);
    break;
    case 5 :
    return prefac*I13(k*r)*k;
    break;
    case 6 :
    return prefac*I31(k*r)*k;
    break;
  }
}
// Integrate them over k_i
static double a13i(int a, double r, const PowerSpectrum& P) {
    return 1./pow2(M_PI)*Integrate(bind(a13_terms, a, r, cref(P), _1), KMIN_e, KMAX_e, error3);
}
// Those with u_lr^2 factor
static double a13u(double r, const PowerSpectrum& P) {
  double myarray[6];
  for (int i=0; i<6; i++){
  myarray[i]=a13i(i,r,cref(P));
  }
  return (20./7.*pow2(myarray[0]) - 4.*myarray[1]*myarray[2]
       + myarray[3]*myarray[4] + 8./7.*pow2(myarray[5]));
}
// Those without u_lr^2 factor
static double a13nou(double r, const PowerSpectrum& P) {
  double myarray[3];
  myarray[0]=a13i(3,r,cref(P));
  myarray[1]=a13i(4,r,cref(P));
  myarray[2]=a13i(6,r,cref(P));

  return -myarray[0]*myarray[1] + 4./7.*pow2(myarray[2]);
}

// Spline F2(k,p,x) over x for MC integral
Spline f2_xp;

inline void f2_xp_init(int k_index,int y_index){
		vector<double> x_val,f2_val;
		double f2;
		for(int i = 0;i<256;i++){
			f2 = A13_F2A[k_index][y_index][i];
			x_val.push_back(x128[i]);
			f2_val.push_back(f2);
		}
		f2_xp = LinearSpline(x_val,f2_val);
}

//Numerical A13 - Monte Carlo integrated over angles (using CUBA library)
// Copter already has Cuba algorithm integrator - use that instead of below ?
static int a13_angular_integrand(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) {

	double r = *((double*)userdata);
	double k = *((double*)userdata+1);
	double y = *((double*)userdata+2);
	double k_int = *((double*)userdata+3);
	double y_int = *((double*)userdata+4);

	  #define u2j xx[0]
	  #define u3j xx[1]
	  #define o2j xx[2]
	  #define o3j xx[3]
	  #define f ff[0]
    #define f1 ff[1]

	  // Redefinition of variables to unit hypercube for Cuba's routine and defining u_i = cos(theta_i)
    real u2 = 2.*u2j-1.;
    real u3 = 2.*u3j-1.;
    real o2 = o2j*2.*M_PI;
    real o3 = o3j*2.*M_PI;

      double so2 = sin(o2);
      double so3 = sin(o3);
      double su2 = fabs(sqrt(1.-pow2(u2)));
      double su3 = fabs(sqrt(1.-pow2(u3)));

    //  k1.k2
      double x = u2*u3+su2*su3*so2*so3;
    // make sure it's within bounds
      if (x<=x128[0] || x>=x128[255]) {
          f=0.;
    			f1=0.;
    		}
    		else{
     //real part of complex exponentials
	double  a1 = cos(r*(k*u3+k*y*u2));
    // k_i.l terms
        double  c = u2*u3-su2*so2*su3*so3; // with mu
        double  c1 = su2*so2*su3*so3; // without mu

      int z = searchnearest(x128,x);
      int k_index = (int)round(k_int);
      int y_index = (int)round(y_int);
     double b1 = A13_F2A[k_index][y_index][z];

    // 16 Pi^2 factor is because integration is done on unit hypercube int_0^1.
    f  = 16.*pow2(M_PI)*a1*c*b1;//f2_xp(x);
    f1 = 16.*pow2(M_PI)*a1*c1*b1;//f2_xp(x);
  }
   return 0;
	}


// r , y index, with mu_lr (1) and without (2), k, y
static int a13_angular(double r, double k, int k_int, double y, int y_int, double myarray[]){
  // if (fabs(lin_check[0])<5e-3 && fabs(lin_check[1])<5e-3){
  //   myarray[0] =  lin_check[0]/fabs(lin_check[0])*7e5/pow4(r);
  //   myarray[1] =  lin_check[1]/fabs(lin_check[1])*7e5/pow4(r);
  //   return 0;
  // }

  if (r>135) {
    myarray[0] =  -7e5/pow4(r);
    myarray[1] =  1.5e5/pow4(r);
      return 0;
 }

  #define NDIM 4
  #define NCOMP 2
  #define USERDATA NULL
  #define NVEC 1
  #define VERBOSE 0
  #define LAST 4
  #define SEED 0
  #define MINEVAL 0

    #define NSTART 1000
    #define NINCREASE 500
    #define NBATCH 1000
    #define GRIDNO 9
    #define STATEFILE NULL
    #define SPIN NULL

    #define KEY 0

  	#define KEY1 47
  	#define KEY2 1
  	#define KEY3 1
  	#define MAXPASS 3
  	#define BORDER 0.
  	#define MAXCHISQ 10.
  	#define MINDEVIATION .25
  	#define NGIVEN 0
  	#define LDXGIVEN NDIM
  	#define NEXTRA 0

  	int comp, nregions, neval, fail;
  	cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];
    double MAXEVAL, EPSREL, EPSABS;

  	void *userdata;
    double n[5];
    n[0] = r;
    n[1] = k;
    n[2] = y;
    n[3] = k_int;
    n[4] = y_int;
  	userdata = &n;
    //
    // if (r<35.) {
    // 	MAXEVAL = 80000;
    // 	EPSREL = 1e-3;
    // 	EPSABS = 1e-6;
    //
  	//    Divonne(NDIM, NCOMP, a13_angular_integrand, userdata, NVEC,
    //       	EPSREL, EPSABS, VERBOSE, SEED,
    //         MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
    //         BORDER, MAXCHISQ, MINDEVIATION,
    //         NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
    //         STATEFILE, SPIN,
    //         &nregions, &neval, &fail, integral, error, prob);
    //       }
  //  else{
      if (r<40) {
     MAXEVAL = 200000;
     EPSREL = 1e-4;
         }
     else if (r>=40 && r<80){
       MAXEVAL = 350000;
       EPSREL = 1e-5;
     }
     else{
       MAXEVAL = 500000;
       EPSREL = 1e-6;
     }
     EPSABS = 1e-12;
           Vegas(NDIM, NCOMP, a13_angular_integrand, userdata, NVEC,
                EPSREL, EPSABS,VERBOSE, SEED,
                MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                GRIDNO, STATEFILE, SPIN,
                &neval, &fail, integral, error, prob);
  //    }
        myarray[0] = (double)integral[0];
        myarray[1] = (double)integral[1];
    return 0;
}
//a = 1 for scale dependance, a = 0 for scale independance
static void a13_integrand_num(int a, double r, double array13[], const PowerSpectrum& P){
  const int nk = 100;
  int y1,y2,param1,param2;
  double KMIN,KMAX,resu,resnou,G1c,k1,k2;
  double k[nk+1];
  double integrand_u[nk+1],integrand_nou[nk+1];
     for(int j = 0; j<=nk; j++){
     k[j] =  KMIN_e * exp(j*log(KMAX_e/KMIN_e)/(nk-1.)); // taking half of A10_num k values
     KMAX = KMAX_e/k[j];
     KMIN = KMIN_e/k[j];
     if (a==1) {
       k1 = k[j];
       k2 = k[j];
     }
     else{
       k1 = KMIN_e;
       k2 = KMAX_e;
     }
     y1 = num_selec_mag(k1, k2, KMIN);
     y2 = num_selec_mag(k1, k2, KMAX);
     double integrandy_u[y2-y1+1];
     double integrandy_nou[y2-y1+1];
     double y[y2-y1+1];
    for (int i = y1; i<=y2; i++){
      y[i-y1] = KMIN_e/k2 * exp(i*log(KMAX_e*k2/(KMIN_e*k1))/(n2*1.-1.));
      // 1st order kernels G1(p)/F1(p)
      G1c = A13_G1P[j][i];
      double myarray[2];
      a13_angular(r, k[j],j, y[i-y1],i ,myarray);
      integrandy_u[i-y1] = G1c*y[i-y1]*P(k[j]*y[i-y1])*myarray[0];
      integrandy_nou[i-y1] = G1c*y[i-y1]*P(k[j]*y[i-y1])*myarray[1];
    }
    resu = 0.;
    resnou = 0.;
    // Integrate over y using trap rule
    for( int i = y1+1; i <= y2; ++ i ){
        resu += 0.5 * (y[i-y1] - y[i-y1-1])*(integrandy_u[i-y1] + integrandy_u[i-y1-1]);
        resnou += 0.5 * (y[i-y1] - y[i-y1-1])*(integrandy_nou[i-y1] + integrandy_nou[i-y1-1]);
    }
    integrand_u[j] =   A13_G1K[j]*P(k[j])*pow3(k[j])*resu;
    integrand_nou[j] = A13_G1K[j]*P(k[j])*pow3(k[j])*resnou;
    }
     resu = 0.;
     resnou = 0.;
     // Integrate over k using trap rule
     for( int i = 1; i <nk; ++ i ){
     resu += 0.5 * (k[i] - k[i-1])*(integrand_u[i] + integrand_u[i-1]);
     resnou += 0.5 * (k[i] - k[i-1])*(integrand_nou[i] + integrand_nou[i-1]);
     }
     double prefac = -1./pow4(dnorm_spt)/pow6(M_PI)/16.;

     array13[0]=resu*prefac;
     array13[1]=resnou*prefac;
}


//Mean infall velocity
static double v12_integrand(int a, double b,double r, double k){
  switch (a) {
    case 1:
	   return  SphericalBesselJ1(k*r)*(pow2(b)*a56(k) + b*k*dtloopy(k)/pow2(M_PI));
     break;
    case 2:
     return SphericalBesselJ1(k*r)*(pow2(b)*a56n(k) + b*k*dtloopy_num(k)/pow2(M_PI));
     break;

}
}

// Non linear in fall velocity : a = 1 gives analytic, a=2 gives numerical

real CorrelationFunction::v12(int a, double b, double r) const{
  double divisor;
  if (a==1) {
    divisor = (1+pow2(b)*xi_real(r,1));
  }
  else{
    divisor = (1+pow2(b)*xi_real(r,3));
  }
	return Integrate(bind(v12_integrand,a,b, r, _1), KMIN_e , KMAX_e , error1)/divisor; // maybe change to non linear xi?
}


//Velocity dispersion terms
//Note u = u_l^2
//Split into u_l^2 scaled part and other part


// A9,A12(Bottom),A14,A15
// mybessel1 is just j0(x)-3j1(x)/x
// u_lr^2 part
static double sig_integrandu(int a, double r, const PowerSpectrum& P, double k ){
  switch (a) {
    case 1:
      return  -mybessel1(k*r)*(ttloopy(k)/pow2(M_PI)+ a12b(k)-a14(k) - a15(k));
      break;
    case 2:
      return -mybessel1(k*r)*(ttloopy_num(k)/pow2(M_PI)+ a12bn(k)-a14n(k) - a15n(k));
      break;
    }
}

// mybessel2 is just j1(x)/x
//without u_lr^2 part
static double sig_integrandnou(int a, double r, const PowerSpectrum& P, double k ){
  switch (a) {
    case 1:
      return -mybessel2(k*r)*(ttloopy(k)/pow2(M_PI)+ a12b(k)- a14(k) - a15(k))
              + pow2(fl_spt*Dl_spt/dnorm_spt)*P(k)/3./pow2(M_PI) // extra term from line 29
              + SphericalBesselJ0(k*r)*a12a(k);
    case 2:
      return -mybessel2(k*r)*(ttloopy_num(k)/pow2(M_PI)+ a12bn(k)- a14n(k) - a15n(k))
              + pow2(g1k_num(k)/dnorm_spt)*P(k)/3./pow2(M_PI) // extra term from line 29
              +SphericalBesselJ0(k*r)*a12an(k);
            }
}


//Linear velocity dispersion  : GR (analytical) (can't combine with the above because of extra bias factor)
static double sigl_integrandu(int a, double r, const PowerSpectrum& P, double k){
  switch (a) {
    case 1:
        return  -pow2(fl_spt*Dl_spt/(dnorm_spt*M_PI))*P(k)*mybessel1(k*r);
    break;
    case 2:
        return  -pow2(g1k_num(k)/(dnorm_spt*M_PI))*P(k)*mybessel1(k*r);
    break;
  }
}
static double sigl_integrandnou(int a, double r, const PowerSpectrum& P, double k){
  switch (a) {
    case 1:
      return  pow2(fl_spt*Dl_spt/(dnorm_spt*M_PI))*P(k)*(1./3.-mybessel2(k*r));
    break;
    case 2:
      return  pow2(g1k_num(k)/(dnorm_spt*M_PI))*P(k)*(1./3.-mybessel2(k*r));
    break;
  }
}


//Integrate them ....
static double sig_u(int a, double r, const PowerSpectrum& P){
return Integrate(bind(sig_integrandu, a, r, cref(P), _1), KMIN_e, KMAX_e, error1);
}
static double sig_nou(int a, double r, const PowerSpectrum& P){
return Integrate(bind(sig_integrandnou, a, r, cref(P), _1), KMIN_e, KMAX_e, error1);
}

static double sigl_u(int a, double r, const PowerSpectrum& P){
return Integrate(bind(sigl_integrandu, a, r, cref(P), _1), KMIN_e, KMAX_e, error1);
}
static double sigl_nou(int a, double r, const PowerSpectrum& P){
return Integrate(bind(sigl_integrandnou, a, r, cref(P), _1), KMIN_e, KMAX_e, error1);
}


static void choose_a13(int a, int b, double r, double array13[], const PowerSpectrum& P){
  switch (a) {
    case 1:
        array13[0] = a13u(r, cref(P));
        array13[1] = a13nou(r, cref(P));
      break;
    case 2:
      a13_integrand_num(b, r, array13, cref(P));
      break;
  }
}


double a10_const;

// Splining of sig_12^2 terms as functions of r
Spline sigmu, signomu, siglmu, siglnomu, xi_real_spline_l,xi_real_spline_nl,v12_spline,v12l_spline;
//a =1 analytic, a =2 numeric
// b= 0 for scale indep, b =1 for scale dep (a=2 case only)
// Takes 11.5 mins to initialize with Loop terms (kernels are initialized 62 + 80 times - can reduce this to 62)
void CorrelationFunction::sig_init(int a, int b)const{
  const int loop_N = 128; // matched to 1% of converged value of Integrate[sig,2,250] for analytic and numerical of GR
  const double lim1 = 2.;
  const double lim2 = 250.;
  double sig0[loop_N],sig1[loop_N],sig2[loop_N],sig3[loop_N],sig4[loop_N],sig5[loop_N],sig6[loop_N],sig7[loop_N],r[loop_N];
  vector<double> rval_table;
  vector<double> sig_table0,sig_table1,sig_table2,sig_table3,sig_table4,sig_table5,sig_table6,sig_table7;
  //#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i<loop_N ; i++) {
     double array13[2];
     double lado_test;
     r[i] = lim1 + i*(lim2-lim1)/(1.*loop_N-1.);
  //   choose_a13(a,b,r[i],array13,cref(P));
     // mean infall velocity
     sig0[i] = v12(a,1.,r[i]);
     sig7[i] = v12_L(a,1.,r[i]);
     //  perpendicular component of sig_12
     sig1[i] = sig_nou(a,r[i],cref(P))+ array13[1];

     // parallel component - perpendicular compoent
     sig2[i] = sig_u(a,r[i],cref(P)) + array13[0] - 1./2.*pow2(sig7[i]);

     //Linear vel dispersion terms
     sig3[i] = sigl_u(a,r[i],cref(P));
     sig4[i] = sigl_nou(a,r[i],cref(P));
     //Correlation functions (linear and non linear)
     int xi_selec1 = 2*a-1;
     int xi_selec2 = 2*a;
     sig5[i] = xi_real(r[i],xi_selec1);
     sig6[i] = xi_real(r[i],xi_selec2);
     // Parallel component
    //   //perp and parallel sig_12 components can be negative in smallest rbins.
    //   //Lado solution is to set them to a small positive value:
      if(sig1[i] <= 0.) {
        sig1[i] = 0.1;
        }
      lado_test = sig2[i] + sig1[i];
      if(lado_test <= 0.) {
         sig2[i] = 0.1-sig1[i];
         }
    }
    for (int i = 0; i<loop_N ; i++) {
     sig_table0.push_back(sig0[i]);
     sig_table1.push_back(sig1[i]);
     sig_table2.push_back(sig2[i]);
     sig_table3.push_back(sig3[i]);
     sig_table4.push_back(sig4[i]);
     sig_table5.push_back(sig5[i]);
     sig_table6.push_back(sig6[i]);
     sig_table7.push_back(sig7[i]);
     rval_table.push_back(r[i]);
    }
    v12_spline = LinearSpline(rval_table,sig_table0);
    signomu = LinearSpline(rval_table, sig_table1);
    sigmu = LinearSpline(rval_table, sig_table2);
    siglmu = LinearSpline(rval_table, sig_table3);
    siglnomu = LinearSpline(rval_table, sig_table4);
    xi_real_spline_l = LinearSpline(rval_table,sig_table5);
    xi_real_spline_nl = LinearSpline(rval_table,sig_table6);
    v12l_spline = LinearSpline(rval_table,sig_table7);

    if (a==1){
    a10_const = a10(cref(P));
      }
    else
    a10_const = a10_num();
  }


// Velocity dispersion.
// Other code does not include linear vel dispersion nor division by (1+xi) - have excluded them
real CorrelationFunction::sig12(double b, double r, double u) const{
  	return   b*(u*sigmu(r) + signomu(r)+ a10_const);
 	}

// Linear velocity dispersion
real CorrelationFunction::sig12_L(double r, double u) const{
	return u*siglmu(r)+siglnomu(r);
}

  // GAUSSIAN STREAMING MODEL CORRELATION FUNCTION //

  // r_sigma goes from s_min to s_max
  // r_pi goes from -s_max to s_max
  // y goes from -s_max-50 to s_max +50  (50 = y_spanning)

real CorrelationFunction::gsm_integrand_inner(int a, double b, double r_sigma, double r_pi, double y) const{
      double r = sqrt( pow( r_sigma, 2 ) + pow( y, 2 ) );
      double mu = y / r;
      if( r < 2 ){
  	   return 0.;
     }
      double musqr = pow2(mu);
      double sigma2 = sig12(b,r,musqr);
      double exponent_gsm = r_pi - y - mu * v12(a,b,r);
      if( sigma2 < 0. ){
        throw "Error: Problem with sign of s_12^2";
      }
      double exp_idx = -0.5 / sigma2 * pow2(exponent_gsm);
      double res = ( 1. + pow2(b)*xi_real_spline_nl(r))*exp( exp_idx )/sqrt( 2 * M_PI * sigma2 );
      return res;
  }

 // Integration stuff for y
  const double y_spanning =50.;
  const double dy_gsm = 0.5;
  const int lim_gsm = 2*y_spanning/dy_gsm;
  double integrand_gsm[lim_gsm];

  // r_sigma goes from s_min to s_max
  // r_pi goes from -s_max to s_max
  // y goes from -s_max-50 to s_max +50
void CorrelationFunction::gsm_xi_init(int a, double b, double s, double mu_s)const{
      double r_sigma = s * sqrt( 1 - pow( mu_s, 2 ) );
      double r_pi    = s * mu_s;
      double y[lim_gsm];
      y[0] = r_pi-y_spanning;
  for (int i = 1; i<lim_gsm; i++){
  y[i] = y[i-1]+dy_gsm;
  }
  //#pragma omp parallel for schedule(dynamic)
  for ( int i = 0; i<lim_gsm ; i++ ) {
  integrand_gsm[i] = gsm_integrand_inner(a,b, r_sigma, r_pi, y[i]);
  }
  }
//Calculate the gsm correlation function
real CorrelationFunction::gsm_xi(int a, double b, double s, double  mu_s )const{
      double r_sigma = s * sqrt( 1 - pow( mu_s, 2 ) );
      double r_pi    = s * mu_s;
      integral intg;
      double y = r_pi - y_spanning;
      gsm_xi_init(a, b,s,mu_s);
      for ( int i = 0; i<lim_gsm; i++ ) {
        intg.read( y, integrand_gsm[i]);
          y += dy_gsm;
        }
        return intg.result(  ) - 1.;
      }

// GSM multipoles : a = 1 gives analytic, a=2 gives numerical

real CorrelationFunction::gsm_multi(int a, int order, double b, double s ) const{
  double myresult = 0;
    for( int i = 0; i < 128; i++ )
    {
   double mu_s = integral::x[i];
   double temp_xi = gsm_xi(a,b, s, mu_s )* legendre( order, mu_s );
   myresult += integral::w[i] * temp_xi;
    }
    return myresult* ( 2. * order + 1. )*0.5 ;
  }



  // GAUSSIAN STREAMING MODEL CORRELATION FUNCTION //

  // r_sigma goes from s_min to s_max
  // r_pi goes from -s_max to s_max
  // y goes from -s_max-50 to s_max +50  (50 = y_spanning)

real CorrelationFunction::lsm_integrand_inner(double b, double r_sigma, double r_pi, double y) const{
      double r = sqrt( pow( r_sigma, 2 ) + pow( y, 2 ) );
      double mu = y / r;
      if( r < 2 ){
  	   return 0.;
     }
      double musqr = pow2(mu);
      double sigma2 = sig12_L(r,musqr);
      double v12_L = v12l_spline(r);
      double yrpi = r_pi - y;
      if( sigma2 < 0. ){
        throw "Error: Problem with sign of s_12^2";
      }
      double exp_idx = -0.5 / sigma2 * pow2(yrpi-mu*v12_L);
      double res = (1. + pow2(b)*xi_real_spline_l(r) +
       (mu*yrpi*v12_L/sigma2 - 1./4.*musqr*pow2(v12_L)/sigma2*(1.-pow2(yrpi)/sigma2)))
                    *exp( exp_idx )/sqrt( 2 * M_PI * sigma2 );
      return res;
  }

 // Integration stuff for y
  double integrand_lsm[lim_gsm];

  // r_sigma goes from s_min to s_max
  // r_pi goes from -s_max to s_max
  // y goes from -s_max-50 to s_max +50
void CorrelationFunction::lsm_xi_init(double b, double s, double mu_s)const{
      double r_sigma = s * sqrt( 1 - pow( mu_s, 2 ) );
      double r_pi    = s * mu_s;
      double y[lim_gsm];
      y[0] = r_pi-y_spanning;
  for (int i = 1; i<lim_gsm; i++){
  y[i] = y[i-1]+dy_gsm;
  }
  //#pragma omp parallel for schedule(dynamic)
  for ( int i = 0; i<lim_gsm ; i++ ) {
  integrand_lsm[i] = lsm_integrand_inner(b, r_sigma, r_pi, y[i]);
  }
  }
//Calculate the gsm correlation function
real CorrelationFunction::lsm_xi(double b, double s, double  mu_s )const{
      double r_sigma = s * sqrt( 1 - pow( mu_s, 2 ) );
      double r_pi    = s * mu_s;
      integral intg;
      double y = r_pi - y_spanning;
      lsm_xi_init(b,s,mu_s);
      for ( int i = 0; i<lim_gsm; i++ ) {
        intg.read( y, integrand_lsm[i]);
          y += dy_gsm;
        }
        return intg.result(  ) - 1.;
      }

//LSM multipoles

real CorrelationFunction::lsm_multi(int order, double b, double s ) const{
  double myresult = 0;
    for( int i = 0; i < 128; i++ )
    {
   double mu_s = integral::x[i];
   double temp_xi = lsm_xi(b, s, mu_s )* legendre( order, mu_s );
   myresult += integral::w[i] * temp_xi;
    }
    return myresult* ( 2. * order + 1. )*0.5 ;
  }
