#if HAVE_CONFIG_H
# include <config.h>
#endif

#include "RegPT.h"
#include "SPT.h"
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
RegPT::RegPT(const Cosmology& C, const PowerSpectrum& P_l, real epsrel_)
: C(C), P_l(P_l)
{
    epsrel = epsrel_;
}

//Error in exponential index integral
const double error1 = 1e-3;
//Atsushi limit on k-angular integration
inline double ATS(real k, real r){
 real KMAX = QMAXp/k;
 real KMIN = QMINp/k;
 if(r>=0.5) {
 return 1./(2.*r);
 }
 else{
 return  Min(XMAX, (1.+r*r-KMIN*KMIN)/2./r);
 }
}

inline int num_selec_mag(double kmin, double kmax, double y){
    	int mag_int;
    //	 mag_int = (n2*1.-1)*(y-QMINp/0.5)*0.005/QMAXp + 1./n2; // linear
    //	 mag_int = sqrt((y-QMINp/kmax)*kmin/QMAXp)*(n2-1); //quadratic
    mag_int = (int)round((n2-1)*log((y*kmax)/QMINp)/log(QMAXp*kmax/(kmin*QMINp))); //exponential
    		return mag_int;
}

// spline linear growth for sigma_d^2 integral
Spline F1_reg;
// Regularized PT exponential
static real regpt_exp(const PowerSpectrum& P_l, double q){
  return  P_l(q)/(6.*pow2(M_PI));
}
// Initialization of sigma_d damping term
Spline damp_term;
void RegPT::sigmad_init() const{
vector<double> kval_table, sigd_table;
double k,sigd;
int n3 = 200;
  for(int i = 0; i<n3; i++ ){
    k = QMINp*exp(i*log(2.*QMAXp/(QMINp))/(n3-1.));
    sigd = pow2(k)*Integrate<ExpSub>(bind(regpt_exp, cref(P_l), _1), QMINp, k/2., error1);
    kval_table.push_back(k);
    sigd_table.push_back(sigd);
      }
    damp_term = LinearSpline(kval_table,sigd_table);
  }


  // used to treat f as free parameter in GR template
double rempfr;
void RegPT::rempreg(real f) const{
rempfr = f/fl_spt;
}

// normalization of linear growth
double sig_8r;
// function to set normalization of linear growth to fiducial for data comparisons
void RegPT::Greg(real dfid) const {
  sig_8r = dfid/D_spt;
       }

/* P_{ab}(k) 1-loop in regpt: analytical LCDM and nDGP */
// 1:  P_dd, 2: P_dt , 3: P_tt ;
real RegPT::PLOOPr(int a, real k) const {
  SPT spt(C,P_l,1e-3);
  double D0,D1,D2,D3;
  if(pow2(D_spt/dnorm_spt)*damp_term(k) >= 80. ){
    return 0;
  }
  else{
  switch (a) {
    case 1:
      D0= pow2(sig_8r)*pow2(D_spt/dnorm_spt);
      D1 = pow4(sig_8r)*spt.P13D(k,1);
      D2 = pow2(D1/P_l(k)/2.)/D0;
      D3 = pow4(sig_8r);
      break;
    case 2:
      D0= -rempfr*fdgp_spt*pow2(sig_8r)*pow2(D_spt/dnorm_spt);
      D1= rempfr*pow4(sig_8r)*spt.P13D(k,2);
      D2 =  pow6(sig_8r)*spt.P13D(k,1)*spt.P13D(k,3)*pow2(dnorm_spt/(D_spt*P_l(k)*2.))/(-rempfr*fdgp_spt);
      D3 =  pow4(sig_8r)*rempfr;
      break;
    case 3:
      D0=  pow2(sig_8r)*pow2(rempfr*fdgp_spt*D_spt/dnorm_spt);
      D1 = rempfr*rempfr* pow4(sig_8r)*spt.P13D(k,3);
      D2 = pow2(D1/P_l(k)/2.)/D0;
      D3 =  pow4(sig_8r)*pow2(rempfr);
      break;
  }
  return exp(-pow2(D_spt/dnorm_spt)*damp_term(k))*(D0*P_l(k)+D1+ D3*spt.P22D(k,a)+pow2(D_spt/dnorm_spt)*damp_term(k)/2.*(2.*D0*P_l(k)+ D1) + pow2(pow2(D_spt/dnorm_spt)*damp_term(k))/4.*D0*P_l(k) + P_l(k)*D2);
}
}

/* P_{ab}(k) 1-loop regpt :  numerical for arbitrary model of gravity - kernel dependent */

// kmin and kmax are only used if the gravity model is scale independant and kernels only need be initialized once
//kmin and kmax then refer to the limits of the desired output range of scales


real RegPT::PLOOPnr(real kmin, real kmax, int a, real k) const {
  SPT spt(C,P_l,1e-3);
  double D0;
  double D1;
  double D2;
  if(pow2(F1_nk/dnorm_spt)*damp_term(k) >= 80. ){
    return 0;
  }
  else{
  switch (a) {
    case 1:
      D0= pow2(F1_nk/dnorm_spt);
      D1 = spt.P13n(kmin,kmax,1,k);
      D2 = pow2(D1/P_l(k)/2.)/D0;
      break;
    case 2:
      D0= F1_nk*G1_nk/pow2(dnorm_spt);
      D1 =  spt.P13n(kmin,kmax,2,k);
      D2 = spt.P13n(kmin,kmax,1,k)*spt.P13n(kmin,kmax,3,k)*pow2(1./(P_l(k)*2.))/D0;
      break;
    case 3:
      D0=  pow2(G1_nk/dnorm_spt);
      D1 =  spt.P13n(kmin,kmax,3,k);
      D2 = pow2(D1/P_l(k)/2.)/D0;
      break;
  }
  return exp(-pow2(F1_nk/dnorm_spt)*damp_term(k))*(D0*P_l(k) + D1 + spt.P22n(kmin,kmax,a,k) + pow2(F1_nk/dnorm_spt)*damp_term(k)/2.*(2.*D0*P_l(k)+ D1) + pow2(pow2(F1_nk/dnorm_spt)*damp_term(k))/4.*D0*P_l(k) + P_l(k)*D2);
  }
}


/* TNS  Multipoles  */
//bl is linear galaxy bias
// a = 1 : Monopole
// a = 2 : Quadrupole
// a = 3 : Hexdecapole
real RegPT::PTNSMnDGPr(real k, real bl, real sigma_v, int a) const {
  if(pow2(D_spt/dnorm_spt)*damp_term(k) >= 80. ){
    return 0;
  }
  else{
	return  (pow2(bl)*factL(k,sigma_v,fl_spt,1.,0,a,6)*PLOOPr(1, k) - 2*factL(k,sigma_v,fl_spt,1.,1,a,6)*bl*PLOOPr(2, k) + factL(k,sigma_v,fl_spt,1.,2,a,6)* PLOOPr(3, k)+  pow4(sig_8r)*ABr(k, bl, sigma_v,a));
}
}


/*  Arbitrary model TNS  Multipoles with DFoG term  */
//bl is linear galaxy bias
// a = 1 : Monopole
// a = 2 : Quadrupole
// a = 3 : Hexdecapole
real RegPT::PTNSMmgr(real k, double kmin, double kmax, int a,real bl, real sigma_v1 ) const {
  if(pow2(F1_nk/dnorm_spt)*damp_term(k) >= 80. ){
    return 0;
  }
  else{
	return  (pow2(bl)*factL(k,sigma_v1,fl_spt,1.,0,a,6)*PLOOPnr(kmin,kmax,1, k)  - 2*factL(k,sigma_v1,fl_spt,1.,1,a,6)*bl*PLOOPnr(kmin,kmax,2, k)  + factL(k,sigma_v1,fl_spt,1.,2,a,6)* PLOOPnr(kmin,kmax,3, k)+  ABnr(kmin,kmax, bl, k, sigma_v1, a));
}
}




//REDSHIFT SPACE MODEL-TNS TERMS with damping

/* Non vanishing A, B and C terms */

// Notes for analytic ABC terms : DEFINITION OF THETA= - DEL V / aHf (and is accounted for in definition of Cross Bispectrum)
// See http://arxiv.org/pdf/1006.0699v1.pdf for derivation

// NOTE ON PARAMETERS: x =k1.k/(k^2*r)


static real ABC(int a, const PowerSpectrum& P_L, real k, real r, real x) {
    real d = 1 + r*r - 2*r*x;
    double dampexpa = exp(-pow2(D_spt/dnorm_spt)*damp_term(k*r)/2.-pow2(D_spt/dnorm_spt)*damp_term(k*sqrt(d))/2.-pow2(D_spt/dnorm_spt)*damp_term(k)/2.);
    double dampexpb = exp(-pow2(D_spt/dnorm_spt)*damp_term(k*sqrt(d))-pow2(D_spt/dnorm_spt)*damp_term(k*r));
    if(d < 1e-5)
        return 0;
    else
		switch(a) {
		// A terms
				//A11
			case 1:
				return dampexpa*(-r*r*r/7.)*(x+6*x*x*x+r*r*x*(-3+10*x*x)+r*(-3+x*x-12*x*x*x*x))*P_L(k) * P_L(k*sqrt(d))/pow2(d);
      	//A12
			case 2:
				return  dampexpa*(r*r*r*r/14.)*(x*x-1)*(-1+7*r*x-6*x*x) * P_L(k) * P_L(k*sqrt(d))/pow2(d);
      	//A22
			case 3:
				return  dampexpa*(r*r*r/14.)*(r*r*x*(13-41*x*x)-4*(x+6*x*x*x)+r*(5+9*x*x+42*x*x*x*x))*P_L(k)*P_L(k*sqrt(d))/pow2(d);
      	//A33
			case 4:
				return dampexpa*(r*r*r/14.)*(1-7*r*x+6*x*x)*(-2*x+r*(-1+3*x*x))*P_L(k)*P_L(k*sqrt(d))/pow2(d);
      	//A11t
			case 5:
				return  dampexpa*(1./7.)*(x+r-2*r*x*x)*(3*r+7*x-10*r*x*x)*P_L(k*r)*P_L(k*sqrt(d))/pow2(d);
      	//A12t
			case 6:
				return dampexpa*(r/14.)*(x*x-1)*(3*r+7*x-10*r*x*x)*P_L(k*r)*P_L(k*sqrt(d))/pow2(d);
      	//A22t
			case 7:
				return dampexpa*(1./14.)*(28*x*x+r*x*(25-81*x*x)+r*r*(1-27*x*x+54*x*x*x*x))*P_L(k*r)*P_L(k*sqrt(d))/pow2(d);
      	//A23t
			case 8:
				return dampexpa*(r/14.)*(-x*x+1)*(r-7*x+6*r*x*x)*P_L(k*r)*P_L(k*sqrt(d))/pow2(d);
		//A33t
			case 9:
				return dampexpa*(1./14.)*(-2*x-r+3*r*x*x)*(r-7*x+6*r*x*x)*P_L(k*r)*P_L(k*sqrt(d))/pow2(d);
		//a11
			case 10:
				return dampexpa*(-7*x*x+pow3(r)*x*(-3+10*x*x)+3*r*(x+6*pow3(x))+r*r*(6-19*x*x-8*pow4(x)))*P_L(k)*P_L(k*r)/(7.*d);
		//a12
			case 11:
				return dampexpa*r*(-1+x*x)*(6*r-7*(1+r*r)*x + 8*r*x*x)*P_L(k)*P_L(k*r)/(14.*d);
			//a22
			case 12:
				return dampexpa*(-28*x*x+r*r*r*x*(-13+41*x*x)+r*x*(11+73*x*x)-2*r*r*(-9+31*x*x+20*x*x*x*x))*P_L(k)*P_L(k*r)/(14.*d);
		//a33
			case 13:
				return dampexpa*(7*x + r*(-6+7*r*x-8*x*x))*(-2*x+r*(-1+3*x*x))*P_L(k)*P_L(k*r)/(14.*d);
		//B terms
				//B111
			case 14:
				return dampexpb*(r*r/2.)*(x*x-1)*P_L(k*r)*P_L(k*sqrt(d))/d;
	   //B112
			case 15:
				return dampexpb*(3*r*r/8.)*pow2(x*x-1)*P_L(k*r)*P_L(k*sqrt(d))/d;
		//B121
			case 16:
				return dampexpb*(3.*r*r*r*r/8.)*pow2(x*x-1)*P_L(k*r)*P_L(k*sqrt(d))/pow2(d);
			//B122
			case 17:
				return dampexpb*(5*r*r*r*r/16.)*pow3(x*x-1)*P_L(k*r)*P_L(k*sqrt(d))/pow2(d);
			//B211
			case 18:
				return dampexpb*(r/2.)*(r+2*x-3*r*x*x)*P_L(k*r)*P_L(k*sqrt(d))/d;
		//B212
			case 19:
				return dampexpb*(-3*r/4.)*(x*x-1)*(-r-2*x+5*r*x*x)*P_L(k*r)*P_L(k*sqrt(d))/d;
			//B221
			case 20:
				return dampexpb*(3*r*r/4.)*(x*x-1)*(-2+r*r+6*r*x-5*r*r*x*x)*P_L(k*r)*P_L(k*sqrt(d))/pow2(d);
			//B222
			case 21:
				return dampexpb*(-3*r*r/16.)*pow2(x*x-1)*(6-30*r*x-5*r*r+35*r*r*x*x)*P_L(k*r)*P_L(k*sqrt(d))/pow2(d);
		//B312
			case 22:
				return dampexpb*(r/8.)*(4*x*(3-5*x*x)+r*(3-30*x*x+35*x*x*x*x))*P_L(k*r)*P_L(k*sqrt(d))/d;
		//B321
			case 23:
				return dampexpb*(r/8.)*(-8*x+r*(-12+36*x*x+12*r*x*(3-5*x*x)+r*r*(3-30*x*x+35*x*x*x*x)))*P_L(k*r)*P_L(k*sqrt(d))/pow2(d);
		//B322
			case 24:
				return dampexpb*(3*r/16.)*(x*x-1)*(-8*x+r*(-12+60*x*x+20*r*x*(3-7*x*x)+5*r*r*(1-14*x*x+21*x*x*x*x)))*P_L(k*r)*P_L(k*sqrt(d))/pow2(d);
		//B422 		//B422 // HAVE INCLUDED THE 'TYPO' IN 3RD TERM : 6rr -> 6rrx
			case 25:
				return dampexpb*(r/16.)*(8*x*(-3+5*x*x)-6*r*(3-30*x*x+35*x*x*x*x)+6*r*r*x*(15-70*x*x+63*x*x*x*x)+r*r*r*(5-21*x*x*(5-15*x*x+11*x*x*x*x)))*P_L(k*r)*P_L(k*sqrt(d))/pow2(d);
  	// C terms
			//C11
			case 26:
				return dampexpa*2.*((r*r*r)*(x*x-1)*(Dd_spt*F_spt*r*(r*x-1)-D_spt*Fd_spt*x*d)*P_L(k)*P_L(k*sqrt(d)))/ (Dd_spt*pow2(d));
			//C12
			case 27:
				return -dampexpa*D_spt*Fd_spt*((r*r*r*r)*pow2(x*x-1)*P_L(k)*P_L(k*sqrt(d)))/(Dd_spt*pow2(d));
			//C22
			case 28:
				return dampexpa*((r*r*r)*(x*x-1)*(2*Dd_spt*F_spt*r*(r*x-1)-D_spt*Fd_spt*(r+4*x+2*r*r*x-7*r*x*x))*P_L(k)*P_L(k*sqrt(d)))/ (Dd_spt*pow2(d));
			//C23
			case 29:
				return -dampexpa*D_spt*Fd_spt*((r*r*r*r)*pow2(x*x-1)*P_L(k)*P_L(k*sqrt(d)))/ (Dd_spt*pow2(d));
      //C33
			case 30:
				return dampexpa*D_spt*Fd_spt*((r*r*r)*(x*x-1)*(-2*x+r*(-1+3*x*x))*P_L(k)*P_L(k*sqrt(d)))/ (Dd_spt*pow2(d));
				//C11t
			case 31:
				return dampexpa*2*F_spt*(r*(x*x-1)*(-x+r*(-1+2*x*x))*P_L(k*r)*P_L(k*sqrt(d)))/ pow2(d);
				//C12t
			case 32:
				return -dampexpa*F_spt*(r*r*pow2(x*x-1)*P_L(k*r)*P_L(k*sqrt(d)))/ pow2(d);
				//C22t
			case 33:
				return dampexpa*r*(x*x-1)*(2*D_spt*Fd_spt*(-x+r*(-1+2*x*x)) + Dd_spt*F_spt*(-2*x+r*(-1+3*x*x)))*P_L(k*r)*P_L(k*sqrt(d))/ (Dd_spt*pow2(d));
				//C23t
			case 34:
				return -dampexpa*D_spt*Fd_spt*(r*r*pow2(x*x-1)*P_L(k*r)*P_L(k*sqrt(d)))/ (Dd_spt*pow2(d));
				//C33t
			case 35:
				return dampexpa*D_spt*Fd_spt*(r*(x*x-1)*(-2*x+r*(-1+3*x*x))*P_L(k*r)*P_L(k*sqrt(d)))/ (Dd_spt*pow2(d));
				//c11
			case 36:
				return  -dampexpa*2*r*(1-x*x)*(D_spt*Fd_spt*r*(r*x-1)-Dd_spt*F_spt*x*d)*P_L(k*r)*P_L(k)/(Dd_spt*d);
				//c12
			case 37:
				return -dampexpa*D_spt*Fd_spt*r*r*pow2(-1+x*x)*P_L(k*r)*P_L(k)/(Dd_spt*d);
				//c22
			case 38:
				return dampexpa*r*(-1+x*x)*(-2*Dd_spt*F_spt*x*d+D_spt*Fd_spt*(-2*x+2*r*r*x+3*r*(-1+x*x)))*P_L(k*r)*P_L(k)/(Dd_spt*d);
				//c33
			case 39:
				return dampexpa*D_spt*Fd_spt*r*(-1+x*x)*(-2*x+r*(-1+3*x*x))*P_L(k*r)*P_L(k)/(Dd_spt*d);
			default:
			warning("SPT: invalid indices, a = %d\n", a);
				return 0;
      }
}




/* function to select multipole or u-dependent Redshift PS */
// it gives either u^(2n) or factL(k, u(i.e.sigma_v), n, a, b, e, f) - see factL

// e = 0:  RSD PS
// e = 1 : Monopole
// e = 2 : Quadrupole
// e = 3 : Hexdecapole



// e = 3 : Hexdecapole


static real u0(real k, real u,real F0, int n, int a){
	if (a >= 1) {
		return factL(k,u,F0,1.,n,a,6); }
		else
			return pow(u,2*n);
}


static real midintabc(real bl, int a, const PowerSpectrum& P_L, real k, real u, real r) {
  real KMAX = QMAXp/k;
  real KMIN = QMINp/k;
	real YMIN = Max(XMIN, (1.+r*r-KMAX*KMAX)/2./r);
  real YMAX = ATS(k,r);
	real F0 =rempfr*fdgp_spt;
  real D1 = D_spt;
  real X1 = 1;
  real u01 = u0(k,u,F0,1,a);
  real u02 = u0(k,u,F0,2,a);
  real u03 = u0(k,u,F0,3,a);
  real u04 = u0(k,u,F0,4,a);

	//A term
		return     pow4(D1/dnorm_spt)*2.*(F0*u01*pow2(bl)*(Integrate(bind(ABC,1, cref(P_L), k, r, _1),  YMIN, YMAX ,  1e-3*P_L(k))
										   + Integrate(bind(ABC,5, cref(P_L), k,  r, _1),  YMIN, YMAX ,  1e-3*P_L(k))
										   + Integrate(bind(ABC,10, cref(P_L), k, r, _1),  YMIN, YMAX ,  1e-3*P_L(k)))

							+ F0*F0*u01*bl*(Integrate(bind(ABC,2, cref(P_L), k,  r, _1),  YMIN, YMAX , 1e-3*P_L(k))
										   + Integrate(bind(ABC,6, cref(P_L), k,  r, _1), YMIN, YMAX ,  1e-3*P_L(k))
										   + Integrate(bind(ABC,11, cref(P_L), k,  r, _1), YMIN, YMAX ,  1e-3*P_L(k)))

							+ F0*F0*u02*bl*(Integrate(bind(ABC,3, cref(P_L), k,  r, _1),   YMIN, YMAX , 1e-3*P_L(k))
										   + Integrate(bind(ABC,7, cref(P_L), k,  r, _1),   YMIN, YMAX ,  1e-3*P_L(k))
										   + Integrate(bind(ABC,12, cref(P_L), k,  r, _1),  YMIN, YMAX ,  1e-3*P_L(k)))

					     + F0*F0*F0*u02*(Integrate(bind(ABC,2, cref(P_L), k,  r, _1),  YMIN, YMAX ,  1e-3*P_L(k))
										   + Integrate(bind(ABC,8, cref(P_L), k,  r, _1),   YMIN, YMAX ,  1e-3*P_L(k))
										   + Integrate(bind(ABC,11, cref(P_L), k, r, _1), YMIN, YMAX ,  1e-3*P_L(k)))

						   + F0*F0*F0*u03*(Integrate(bind(ABC,4, cref(P_L), k, r, _1),   YMIN, YMAX , 1e-3*P_L(k))
										   + Integrate(bind(ABC,9, cref(P_L), k,  r, _1),  YMIN, YMAX ,  1e-3*P_L(k))
										   + Integrate(bind(ABC,13, cref(P_L), k,  r, _1),YMIN, YMAX ,  1e-3*P_L(k))))

	//B term
		  +      pow4(D1/dnorm_spt)*2*(F0*F0*pow2(bl)*u01*Integrate(bind(ABC,14,  cref(P_L), k, r, _1), YMIN, YMAX , 1e-3* P_L(k))

						  -F0*F0*F0*bl*u01*(Integrate(bind(ABC,15, cref(P_L), k, r, _1), YMIN, YMAX , 1e-3*P_L(k))
											    +Integrate(bind(ABC,16, cref(P_L),  k,  r, _1),YMIN, YMAX , 1e-3*P_L(k)))

						+F0*F0*F0*F0*u01*Integrate(bind(ABC,17, cref(P_L),  k, r, _1), YMIN, YMAX , 1e-3*P_L(k))

							  +F0*F0*pow2(bl)*u02*Integrate(bind(ABC,18,  cref(P_L), k,  r, _1),YMIN, YMAX ,  1e-3*P_L(k))

						 -F0*F0*F0*bl*u02*(Integrate(bind(ABC,19, cref(P_L), k,  r, _1), YMIN, YMAX , 1e-3*P_L(k))
										      + Integrate(bind(ABC,20, cref(P_L),  k,  r, _1), YMIN, YMAX ,  1e-3*P_L(k)))

					  +F0*F0*F0*F0*u02*Integrate(bind(ABC,21, cref(P_L),  k,  r, _1),YMIN, YMAX , 1e-3*P_L(k))

					    - F0*F0*F0*bl*u03*(Integrate(bind(ABC,22,  cref(P_L), k,  r, _1), YMIN, YMAX ,  1e-3*P_L(k))
											  + Integrate(bind(ABC,23,  cref(P_L), k,  r, _1), YMIN, YMAX , 1e-3*P_L(k)))

					  + F0*F0*F0*F0*u03*Integrate(bind(ABC,24,  cref(P_L), k,  r, _1),YMIN, YMAX ,  1e-3*P_L(k))

							+F0*F0*F0*F0*u04*Integrate(bind(ABC,25,  cref(P_L), k,  r, _1), YMIN, YMAX , 1e-3*P_L(k)))

	//C - nDGP analytic

		  +       X1*D1*D1*2/pow4(dnorm_spt)*(F0*pow2(bl)*u01*(Integrate(bind(ABC,26, cref(P_L), k, r, _1),   YMIN, YMAX ,  1e-3*P_L(k))
										+ Integrate(bind(ABC,31, cref(P_L), k, r, _1),   YMIN, YMAX ,  1e-3*P_L(k))
										+ Integrate(bind(ABC,36, cref(P_L), k, r, _1), YMIN, YMAX ,  1e-3*P_L(k))) +

							F0*F0*bl*u01*(Integrate(bind(ABC,27, cref(P_L), k, r, _1),   YMIN, YMAX ,  1e-3*P_L(k))
										+ Integrate(bind(ABC,32, cref(P_L), k, r, _1),   YMIN, YMAX ,  1e-3*P_L(k))
										+ Integrate(bind(ABC,37, cref(P_L), k, r, _1),  YMIN, YMAX ,  1e-3*P_L(k))) +

						   F0*F0*bl*u02*(Integrate(bind(ABC,28, cref(P_L), k, r, _1),   YMIN, YMAX ,  1e-3*P_L(k))
										+ Integrate(bind(ABC,33, cref(P_L), k, r, _1),   YMIN, YMAX ,  1e-3*P_L(k))
										+ Integrate(bind(ABC,38, cref(P_L), k,  r, _1), YMIN, YMAX ,  1e-3*P_L(k))) +

						F0*F0*F0*u02*(Integrate(bind(ABC,29, cref(P_L), k, r, _1),   YMIN, YMAX ,  1e-3*P_L(k))
										+ Integrate(bind(ABC,34, cref(P_L), k, r, _1),   YMIN, YMAX ,  1e-3*P_L(k))
										+ Integrate(bind(ABC,37, cref(P_L), k,  r, _1),  YMIN, YMAX ,  1e-3*P_L(k))) +

					   F0*F0*F0*u03*(Integrate(bind(ABC,30, cref(P_L), k, r, _1),   YMIN, YMAX ,  1e-3*P_L(k))
									   + Integrate(bind(ABC,35, cref(P_L), k, r, _1),   YMIN, YMAX ,  1e-3*P_L(k))
									   + Integrate(bind(ABC,39, cref(P_L), k,  r, _1),  YMIN, YMAX ,  1e-3*P_L(k))));

}



/* Non vanishing A, B and C terms */
//numerical calculations
// Notes for numerical AB terms : DEFINITION OF THETA=  DEL V / aH (and is accounted for in definition of Cross Bispectrum)
// Integrating them over angle
static real ABCn(double u0x[], int y, real bl, int a, const PowerSpectrum& P_L, real k, real u, real r) {
  double myresult = 0.;
	double temp_abc;
  double abc[10];
  // Set the limits of angular integration
  double KMAX = QMAXp/k;
  double KMIN = QMINp/k;
  double YMIN = Max(XMIN, (1.+r*r-KMAX*KMAX)/2./r);
  double YMAX = ATS(k,r);
  // Multipole factors
  real u01 = u0x[0];//u0(k,u,fl_spt,1,a);
  real u02 = u0x[1];//u0(k,u,fl_spt,2,a);
  real u03 = u0x[2];//u0(k,u,fl_spt,3,a);
  real u04 = u0x[3];//u0(k,u,fl_spt,4,a);
  // Find those corresponding limits in the GL Quad abscissae (See SpecialFunctions.h)
  int x_min = searchnearest(x128, YMIN);
  int x_max = searchnearest(x128, YMAX);
  for( int i = x_min; i <= x_max; i++ )
    				{
    		 		double x = x128[i];
    		 		double d = 1+ r*r - 2*r*x;
    		 		if(d < 1e-5){
    			 	temp_abc=0.;
    		 		}
    		  	else {
              // 1st order kernels F1/G1(k-p)
              abc[0] = F1kmp_nk[i*n2 + y];
              abc[1] = G1kmp_nk[i*n2 + y]/bl;
              // 1st order kernels F1/G1(p)
              abc[2] = F1p_nk[i*n2 + y];
              abc[3] = G1p_nk[i*n2 + y]/bl;
              //symmetrized 2nd order kernels for ps F2/G2(p,k-p)
              abc[4] = F2_nk[i*n2 + y];
              abc[5] = G2_nk[i*n2 + y]/bl;
              //symmetrized 2nd order kernels F2/G2(-p,k)
              abc[6] = F2B_nk[i*n2 + y];
              abc[7] = G2B_nk[i*n2 + y]/bl;
              //symmetrized 2nd order kernels F2/G2(-k,k-p)
              abc[8] = F2C_nk[i*n2 + y];
              abc[9] = G2C_nk[i*n2 + y]/bl;

              double dampexpa = exp(-pow2(F1_nk/dnorm_spt)*damp_term(k*r)/2.-pow2(F1_nk/dnorm_spt)*damp_term(k*sqrt(d))/2.-pow2(F1_nk/dnorm_spt)*damp_term(k)/2.);
              double dampexpb = exp(-pow2(F1_nk/dnorm_spt)*damp_term(k*sqrt(d))-pow2(F1_nk/dnorm_spt)*damp_term(k*r));

// A terms
//u^2
    		  temp_abc =  pow3(bl)*dampexpa/d*(-u01*(F1_nk*r*(-2*abc[8]*abc[1]*r*(r*x-1)+abc[9]*(2*abc[0]*x*d+abc[1]*r*(1-x*x)))*P_L(k)*P_L(k*sqrt(d))
            						+ abc[4]*r*(-2*abc[2]*abc[1]*r*(-1+r*x)+abc[3]*(2*abc[0]*x*d+abc[1]*r*(1-x*x)))*P_L(k*r)*P_L(k*sqrt(d))
            						+ F1_nk*r*(-2*abc[7]*abc[2]*r*(r*x-1)+abc[3]*(2*abc[6]*x*d+abc[7]*r*(1-x*x)))*P_L(k)*P_L(k*r))
//u^4
            -u02*(r*(2*abc[8]*G1_nk/bl*abc[1]*r*(r*x-1)+abc[9]*(G1_nk/bl*(-2*abc[0]*x*d+abc[1]*r*(x*x-1))+F1_nk*abc[1]*(-2*x+r*(-1+3*x*x))))*P_L(k)*P_L(k*sqrt(d))
                           +r*(abc[4]*abc[1]*abc[3]*(-2*x+r*(-1+3*x*x))+abc[5]*(2*abc[2]*abc[1]*r*(r*x-1)+abc[3]*(-2*abc[0]*x*d+abc[1]*r*(x*x-1))))*P_L(k*r)*P_L(k*sqrt(d))
                           +r*(abc[3]*abc[7]*F1_nk*(-2*x-r+3*x*x*r)+2*abc[2]*G1_nk/bl*abc[7]*r*(r*x-1)+G1_nk/bl*abc[3]*(-2*abc[6]*x*d+abc[7]*r*(x*x-1)))*P_L(k)*P_L(k*r))
//u^6
           -u03*(G1_nk/bl*abc[1]*abc[9]*r*(r+2*x-3*r*x*x)*P_L(k)*P_L(k*sqrt(d))
                        +abc[3]*abc[1]*abc[5]*r*(r+2*x-3*r*x*x)*P_L(k*r)*P_L(k*sqrt(d))
                        +abc[3]*G1_nk/bl*abc[7]*r*(r+2*x-3*r*x*x)*P_L(k*r)*P_L(k)))


// B terms
          +pow4(bl)*dampexpb*P_L(k*r)*P_L(k*sqrt(d))/(16.*d*d)*
//u^2
          (u01*(abc[1]*abc[3]*r*r*(x*x-1)*(2*abc[0]*d*(4*abc[2]+3*abc[3]*(x*x-1))
                         +abc[1]*r*r*(x*x-1)*(6*abc[2]+5*abc[3]*(x*x-1))))
//u^4
          +u02*(abc[1]*abc[3]*r*(-4*abc[0]*d*(3*abc[3]*(x*x-1)*(-2*x+r*(5*x*x-1))+abc[2]*(-4*x+r*(-2+6*x*x)))
                					   -3*abc[1]*r*(x*x-1)*(4*abc[2]*(2-6*r*x+r*r*(-1+5*x*x))
                						 +abc[3]*(x*x-1)*(6-30*r*x+5*r*r*(-1+7*x*x)))))
//u^6
          +u03*(abc[1]*abc[3]*r*(3*abc[3]*abc[1]*(x*x-1)*(-8*x+12*r*(5*x*x-1)-20*r*r*x*(7*x*x-3)+5*pow3(r)*(1-14*x*x+21*pow4(x)))
                              +2*abc[0]*abc[3]*(4*x*(3-5*x*x)+pow3(r)*(3-30*x*x+35*pow4(x))+r*(3-54*x*x+75*pow4(x))+r*r*(6*x+40*pow3(x)-70*pow5(x)))
                              +2*abc[2]*abc[1]*(-8*x+12*r*(3*x*x-1)+r*r*(36*x-60*pow3(x))+pow3(r)*(3-30*x*x+35*pow4(x)))))
//u^8
          +u04*(pow2(abc[1]*abc[3])*r*(8*x*(5*x*x-3)-6*r*(3-30*x*x+35*pow4(x))
                     				+6*r*r*x*(15-70*x*x+63*pow4(x))  		// HAVE INCLUDED THE 'TYPO' IN TERM : 6rr -> 6rrx
                     				+pow3(r)*(5-105*x*x+315*pow4(x)-231*pow6(x)))));
    				}
    				myresult += w128[i] * temp_abc;
    				}
          return  2*myresult;
}

/*Generalised A and B correction term */

/*Analytical A and B*/
// a = 1 : LCDM
// a = 2 : MG

// e = 0:  RSD PS
// e = 1 : Monopole
// e = 2 : Quadrupole
// e = 3 : Hexdecapole

//U = u for RSD PS
//U = sigma_v for multipoles

real RegPT::ABr(real k, real bl, real U, int a) const{
    int n3 = 300;
    real KMAX = QMAXp/k;
    real KMIN = QMINp/k;
    double y[n3];
    double integrand[n3];
    for (int i = 0; i<n3; i++){
    y[i] = KMIN * exp(i*log(KMAX/KMIN)/(n3-1.));
    integrand[i] = midintabc(bl, a,cref(P_l), k, U, y[i]);
    }
  double res = 0.;
    for( int i = 1; i < n3; ++ i ){
  res += 0.5 * (y[i] - y[i-1])*(integrand[i] + integrand[i-1]);
      }
    return  k*k*k/(4*M_PI*M_PI) * res;
    }


/*Numerical A and B */

// e = 0:  RSD PS
// e = 1 : Monopole
// e = 2 : Quadrupole
// e = 3 : Hexdecapole

//U = u for RSD PS
//U = sigma_v for multipoles

real RegPT::ABnr(double kmin, double kmax, real bl, real k, real U, int a) const{
	real KMAX = QMAXp/k;
  real KMIN = QMINp/k;
  int y1 = num_selec_mag(kmin, kmax, KMIN);
  int y2 = num_selec_mag(kmin, kmax, KMAX);
  double y[n2];
  double integrand[n2];
  double u0x[4];
  u0x[0]=u0(k,U,fl_spt,1,a);
  u0x[1]=u0(k,U,fl_spt,2,a);
  u0x[2]=u0(k,U,fl_spt,3,a);
  u0x[3]=u0(k,U,fl_spt,4,a);
  for (int i = y1; i<=y2; i++){
  y[i] = QMINp/kmax * exp(i*log(QMAXp*kmax/(QMINp*kmin))/(n2*1.-1.));
  integrand[i] = ABCn(u0x,i, bl, a, cref(P_l), k, U, y[i]);
  }
double res = 0.;
  for( int i = y1+1; i <= y2; ++ i ){
res += 0.5 * (y[i] - y[i-1])*(integrand[i] + integrand[i-1]);
}
  return  k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt) * res;
}


// Separation of numerical A B terms for interpolation

static void ABnr_selec(double abarray[], int y, const PowerSpectrum& P_L, real k, real r) {
  double myresult[14];
  for(int ii=0;ii<14;ii++){
    myresult[ii]=0.;
  }
	double temp_abs[14];
  double abc[10];
  // Set the limits of angular integration
  double KMAX = QMAXp/k;
  double KMIN = QMINp/k;
  double YMIN = Max(XMIN, (1.+r*r-KMAX*KMAX)/2./r);
  double YMAX = ATS(k,r);
  // Find those corresponding limits in the GL Quad abscissae (See SpecialFunctions.h)
  int x_min = searchnearest(x128, YMIN);
  int x_max = searchnearest(x128, YMAX);
  for( int i = x_min; i <= x_max; i++ )
    				{
    		 		double x = x128[i];
    		 		double d = 1+ r*r - 2*r*x;
    		 		if(d < 1e-5){
              for(int ii=0;ii<14;ii++){
                temp_abs[ii]=0.;
              }
    		 		}
    		  	else {
              // 1st order kernels F1/G1(k-p)
              abc[0] = F1kmp_nk[i*n2 + y];
              abc[1] = G1kmp_nk[i*n2 + y];
              // 1st order kernels F1/G1(p)
              abc[2] = F1p_nk[i*n2 + y];
              abc[3] = G1p_nk[i*n2 + y];
              //symmetrized 2nd order kernels for ps F2/G2(p,k-p)
              abc[4] = F2_nk[i*n2 + y];
              abc[5] = G2_nk[i*n2 + y];
              //symmetrized 2nd order kernels F2/G2(-p,k)
              abc[6] = F2B_nk[i*n2 + y];
              abc[7] = G2B_nk[i*n2 + y];
              //symmetrized 2nd order kernels F2/G2(-k,k-p)
              abc[8] = F2C_nk[i*n2 + y];
              abc[9] = G2C_nk[i*n2 + y];

              double dampexpa = exp(-pow2(F1_nk/dnorm_spt)*damp_term(k*r)/2.-pow2(F1_nk/dnorm_spt)*damp_term(k*sqrt(d))/2.-pow2(F1_nk/dnorm_spt)*damp_term(k)/2.);
              double dampexpb = exp(-pow2(F1_nk/dnorm_spt)*damp_term(k*sqrt(d))-pow2(F1_nk/dnorm_spt)*damp_term(k*r));
              double prefacb = dampexpb*P_L(k*r)*P_L(k*sqrt(d))/(16.*d*d);
//A terms
//u^2b^2
          temp_abs[0] =  -dampexpa/d*((F1_nk*r*(-2*abc[8]*abc[1]*r*(r*x-1)+abc[9]*2*abc[0]*x*d))*P_L(k)*P_L(k*sqrt(d))
                                    + abc[4]*r*(-2*abc[2]*abc[1]*r*(-1+r*x)+abc[3]*2*abc[0]*x*d)*P_L(k*r)*P_L(k*sqrt(d))
                                    +F1_nk*r*(-2*abc[7]*abc[2]*r*(r*x-1)+abc[3]*2*abc[6]*x*d)*P_L(k)*P_L(k*r));
//u^2b
          temp_abs[1] =  -dampexpa/d*((F1_nk*r*(abc[9]*abc[1]*r*(1-x*x))*P_L(k)*P_L(k*sqrt(d))
            						            + abc[4]*r*(abc[3]*abc[1]*r*(1-x*x))*P_L(k*r)*P_L(k*sqrt(d))
            						            + F1_nk*r*(abc[3]*abc[7]*r*(1-x*x))*P_L(k)*P_L(k*r)));
//u^4 b
          temp_abs[2] = -dampexpa/d*(r*(2*abc[8]*G1_nk*abc[1]*r*(r*x-1)+abc[9]*(G1_nk*(-2*abc[0]*x*d)+F1_nk*abc[1]*(-2*x+r*(-1+3*x*x))))*P_L(k)*P_L(k*sqrt(d))
                                     +r*(abc[4]*abc[1]*abc[3]*(-2*x+r*(-1+3*x*x))+abc[5]*(2*abc[2]*abc[1]*r*(r*x-1)+abc[3]*(-2*abc[0]*x*d)))*P_L(k*r)*P_L(k*sqrt(d))
                                     +r*(abc[3]*abc[7]*F1_nk*(-2*x-r+3*x*x*r)+2*abc[2]*G1_nk*abc[7]*r*(r*x-1)+G1_nk*abc[3]*(-2*abc[6]*x*d))*P_L(k)*P_L(k*r));
// u^4
          temp_abs[3] = -dampexpa/d*(r*(abc[9]*(G1_nk*(abc[1]*r*(x*x-1))))*P_L(k)*P_L(k*sqrt(d))
                                     +r*(abc[5]*(abc[3]*(abc[1]*r*(x*x-1))))*P_L(k*r)*P_L(k*sqrt(d))
                                     +r*(G1_nk*abc[3]*(abc[7]*r*(x*x-1)))*P_L(k)*P_L(k*r));

//u^6
          temp_abs[4] = -dampexpa/d*(G1_nk*abc[1]*abc[9]*r*(r+2*x-3*r*x*x)*P_L(k)*P_L(k*sqrt(d))
                                  +abc[3]*abc[1]*abc[5]*r*(r+2*x-3*r*x*x)*P_L(k*r)*P_L(k*sqrt(d))
                                  +abc[3]*G1_nk*abc[7]*r*(r+2*x-3*r*x*x)*P_L(k*r)*P_L(k));

//B terms
//u^2 b^2
          temp_abs[5]= prefacb*(abc[1]*abc[3]*r*r*(x*x-1)*(2*abc[0]*d*(4*abc[2])));
//u^2 b
          temp_abs[6]= prefacb*(abc[1]*abc[3]*r*r*(x*x-1)*(2*abc[0]*d*(3*abc[3]*(x*x-1))
                         +abc[1]*r*r*(x*x-1)*(6*abc[2])));
//u^2
          temp_abs[7]= prefacb*(abc[1]*abc[3]*r*r*(x*x-1)*(abc[1]*r*r*(x*x-1)*(5*abc[3]*(x*x-1))));
//u^4b^2
          temp_abs[8]= prefacb*(abc[1]*abc[3]*r*(-4*abc[0]*d*(abc[2]*(-4*x+r*(-2+6*x*x)))));
//u^4 b
          temp_abs[9]= prefacb*(abc[1]*abc[3]*r*(-4*abc[0]*d*(3*abc[3]*(x*x-1)*(-2*x+r*(5*x*x-1)))
                					   -3*abc[1]*r*(x*x-1)*(4*abc[2]*(2-6*r*x+r*r*(-1+5*x*x)))));
//u^4
          temp_abs[10]= prefacb*(abc[1]*abc[3]*r*(-3*abc[1]*r*(x*x-1)*(abc[3]*(x*x-1)*(6-30*r*x+5*r*r*(-1+7*x*x)))));

//u^6 b
          temp_abs[11] = prefacb*(abc[1]*abc[3]*r*(2*abc[0]*abc[3]*(4*x*(3-5*x*x)+pow3(r)*(3-30*x*x+35*pow4(x))+r*(3-54*x*x+75*pow4(x))+r*r*(6*x+40*pow3(x)-70*pow5(x)))
                                +2*abc[2]*abc[1]*(-8*x+12*r*(3*x*x-1)+r*r*(36*x-60*pow3(x))+pow3(r)*(3-30*x*x+35*pow4(x)))));

//u^6
          temp_abs[12]= prefacb*(abc[1]*abc[3]*r*(3*abc[3]*abc[1]*(x*x-1)*(-8*x+12*r*(5*x*x-1)-20*r*r*x*(7*x*x-3)+5*pow3(r)*(1-14*x*x+21*pow4(x)))));

//u^8
          temp_abs[13]= prefacb*(pow2(abc[1]*abc[3])*r*(8*x*(5*x*x-3)-6*r*(3-30*x*x+35*pow4(x))
                               				+6*r*r*x*(15-70*x*x+63*pow4(x))  		// HAVE INCLUDED THE 'TYPO' IN TERM : 6rr -> 6rrx
                               				+pow3(r)*(5-105*x*x+315*pow4(x)-231*pow6(x))));

    				}
            for(int ii=0;ii<14;ii++){
              myresult[ii] += w128[i]*temp_abs[ii];
            }
    				}
            for(int ii=0;ii<14;ii++){
              abarray[ii]=2*myresult[ii];
            }
}


void RegPT::ABr_selec(double myarray[], real k) const {
  	real KMAX = QMAXp/k;
    real KMIN = QMINp/k;
    int y1 = num_selec_mag(k, k, KMIN);
    int y2 = num_selec_mag(k, k, KMAX);
    double y[n2];
    double integrand[n2][14];
    double abarray[14];
    for (int i = y1; i<=y2; i++){
    y[i] = QMINp/k * exp(i*log(QMAXp/QMINp)/(n2*1.-1.));
    ABnr_selec(abarray,i, cref(P_l), k, y[i]);
    for(int ii=0;ii<14;ii++){
    integrand[i][ii] = abarray[ii];
      }
    }
  double res[14];
  for(int ii=0;ii<14;ii++){
    res[ii]=0.;
  }
    for( int i = y1+1; i <= y2; ++ i ){
      for(int ii=0;ii<14;ii++){
        res[ii] += 0.5 * (y[i] - y[i-1])*(integrand[i][ii] + integrand[i-1][ii]);
      }
  }
  for(int ii=0;ii<14;ii++){
    myarray[ii] = k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt) * res[ii];
  }
}
