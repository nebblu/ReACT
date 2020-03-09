#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/bind.hpp>
using boost::cref;

#include "Cosmology.h"
#include "PowerSpectrum.h"
#include "Quadrature.h"
#include "SPT.h"
#include "Spline.h"
#include "SpecialFunctions.h"
#include "LinearPS.h"


#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include <math.h>       /* pow */
#include <boost/math/tools/roots.hpp>
#include <random>


//Atsushi limit on k-angular integration
inline double ATS(real k, real r){
 real KMAX = QMAXp/k;
 real KMIN = QMINp/k;
 if(r>=0.5){
  return 1./(2.*r);
 }
 else {
 return  Min(0.99999999999999, (1.+r*r-KMIN*KMIN)/2./r);
 }
}

//Atsushi limit on k-angular integration
inline double ATSn(real k, real r){
 real KMAX = QMAXp/k;
 real KMIN = QMINp/k;
 if(r>=0.5){
  return 1./(2.*r);
 }
 else {
 return  Min(XMAX, (1.+r*r-KMIN*KMIN)/2./r);
 }
}



////////// PERTURBATION THEORY CONSTRUCTED STATISTICS ////////////

/*Load the linear power spectrum */

SPT::SPT(const Cosmology& C_, const PowerSpectrum& P_L_, real epsrel_)
: C(C_), P_L(P_L_)
{
    epsrel = epsrel_;
}


// normalization of linear growth
double sig_8;

// function to set normalization of linear growth to fiducial for data comparisons
void SPT::G(real dfid) const {
  sig_8 = dfid/D_spt;
       }

double rempf;
void SPT::remp(real f) const{
rempf = f/fl_spt;
}

void SPT::rempn(real f) const{
rempf = f/(-G1_nk/F1_nk);
}

void SPT::Gn(real dfid) const {
sig_8 = dfid/F1_nk;
     }

/* P_{ab}(k) 1-loop: analytical LCDM and nDGP */

/* NOTE THAT THE DEFINITIONS OF VEL DIV ARE THETA= DEL V / (aH) */

/* The division by Dl is my own normalisation of P_L - it requires z=0 always in spt2.cpp */
// 1: P_dd, 2: P_dt , 3: P_tt ; LCDM
// 4:  P_dd, 5: P_dt , 6: P_tt ; nDGP
real SPT::PLOOP(real k, int a) const {
    switch(a) {
        //LCDM
        case 1:
            return  pow2(sig_8*Dl_spt/dnorm_spt)*P_L(k) + pow4(sig_8)*(P13_dd(k)+P22_dd(k));
            break;
        case 2:
            return  -rempf*fl_spt*pow2(sig_8*Dl_spt/dnorm_spt)*P_L(k) + pow4(sig_8)*(P13_dt(k) + P22_dt(k));
            break;
        case 3:
            return  pow2(sig_8*rempf*fl_spt*Dl_spt/dnorm_spt)*P_L(k) +pow4(sig_8)*(P13_tt(k) + P22_tt(k));
            break;
        //DGP
        case 4:
            return   pow2(sig_8*D_spt/dnorm_spt)*P_L(k) + pow4(sig_8)*(P13D_dd(k)+P22D_dd(k));
            break;
        case 5:
            return  -fdgp_spt*pow2(sig_8*D_spt/dnorm_spt)*P_L(k) + pow4(sig_8)*(P13D_dt(k)+ P22D_dt(k));
            break;
        case 6:
            return   pow2(sig_8*fdgp_spt*D_spt/dnorm_spt)*P_L(k) + pow4(sig_8)*(P13D_tt(k)+P22D_tt(k));
            break;
        default:
            warning("SPT: invalid indices, a = %d\n", a);
            return 0;
    }
}


/* P_{ab}(k) 1-loop :  numerical for arbitrary model of gravity - kernel dependent */
real SPT::PLOOPn(double kmin, double kmax,  int a, real k ) const {
    switch(a) {
        case 1:
						return    pow2(sig_8)*pow2(F1_nk/dnorm_spt)*P_L(k) +  pow4(sig_8)*P13n(kmin, kmax, a, k) + pow4(sig_8)*P22n(kmin, kmax, a, k); //delta-delta
            break;
        case 2:
            return    pow2(sig_8)*rempf*F1_nk*G1_nk*P_L(k)/pow2(dnorm_spt) +  pow4(sig_8)*rempf*P13n(kmin, kmax, a, k) + pow4(sig_8)*rempf*P22n(kmin, kmax, a, k); //delta -theta
            break;
        case 3:
            return    pow2(sig_8)*pow2(rempf*G1_nk/dnorm_spt)*P_L(k) +  pow4(sig_8)*rempf*rempf*P13n(kmin, kmax, a, k) + pow4(sig_8)*rempf*rempf*P22n(kmin, kmax, a, k); // theta-theta
            break;
      default:
            warning("SPT: invalid indices, a = %d\n", a);
            return 0;
    }
}


/* P(k,u) RSD KAISER : */
// a =1 : LCDM
// a =2 : nDGP
// a =3 : numerical (arbitrary model)
//bl is linear  bias
real SPT::PRSD(real k, real u, real bl, int a) const {
	double F0,D0;
	switch (a) {
		case 1:
			F0=rempf*fl_spt;
			D0=sig_8*Dl_spt;
			break;
		case 2:
			F0=fdgp_spt;
			D0=D_spt;
			break;
		case 3:
			F0=-G1_nk/F1_nk;
			D0=F1_nk;
			break;
				  }

	return    pow2(D0*bl/dnorm_spt)*(P_L(k) + 2.*(F0/bl)*u*u*P_L(k) + pow2(F0/bl)*pow4(u)* P_L(k));
	 }


/* P(k,u) RSD TNS */
// a = 1 : LCDM
// a = 2 : nDGP
// u = k.z/|k|, z being a unit vector along the line of sight axis
//bl is linear  bias
real SPT::PTNS(real k, real u, real bl, real sigma_v, int a) const {
    switch(a) {
		case 1:
			return  DFOG(k,u,sigma_v,a)*(pow2(bl)*PLOOP(k,1) - 2.*u*u*bl*PLOOP(k,2) + pow4(u)* PLOOP(k,3) + pow4(sig_8)*AB(k, bl, u, a));
      break;
    case 2:
      return	DFOG(k,u,sigma_v,a)*(pow2(bl)*PLOOP(k,4) - 2.*u*u*bl*PLOOP(k,5) + pow4(u)*PLOOP(k,6) + AB(k, bl, u, a));
      break;
        default:
            warning("SPT: invalid indices, a = %d \n", a );
            return 0;
    }
}


/* P(k) RSD TNS: Arbitrary model   */
//bl is linear  bias
// sigma_v is velocity dispersion
real SPT::PTNSn(double kmin, double kmax, real bl, real sigma_v, real k, real u) const {
	return  DFOG(k, u, sigma_v,3 )*(pow2(bl)*PLOOPn(kmin,kmax,1,k) - 2.*u*u*bl*PLOOPn(kmin,kmax,2,k) + pow4(u)* PLOOPn(kmin,kmax,3,k) + ABn(kmin,kmax, bl, k, u, 0));
}


/*Multipoles*/
/*  LCDM/nDGP  Kaiser Multipoles with DFoG term (set sigma_v=0 for linear kaiser multipoles) */
//bl is linear  bias
// a = 1 : Monopole
// a = 2 : Quadrupole
// a = 3 : Hexdecapole
real SPT::KASM(real k, real bl, real sigma_v, int a) const {
	return   pow2(D_spt*bl/dnorm_spt)*(factL(k,sigma_v,fl_spt,1.,0,a,6) + 2.*fdgp_spt/bl*factL(k,sigma_v,fl_spt,1.,1,a,6) + pow2(fdgp_spt/bl)*factL(k,sigma_v,fl_spt,1.,2,a,6))*P_L(k);
}


/*  Arbitrary model Kaiser Multipoles with DFoG term (set sigma_v=0 for linear kaiser multipoles) */
//bl is linear  bias
// a = 1 : Monopole
// a = 2 : Quadrupole
// a = 3 : Hexdecapole
real SPT::KASMmg(real k, real bl, real sigma_v, int a) const {
	return   pow2(F1_nk*bl/dnorm_spt)*(factL(k,sigma_v,fl_spt,1.,0,a,6) *P_L(k) - 2.*(G1_nk/F1_nk/bl)*factL(k,sigma_v,fl_spt,1.,1,a,6) *P_L(k) + pow2(G1_nk/F1_nk/bl)*factL(k,sigma_v,fl_spt,1.,2,a,6) * P_L(k));
}


/*  LCDM  TNS  Multipoles */
//bl is linear bias
// a = 1 : Monopole
// a = 2 : Quadrupole
// a = 3 : Hexdecapole
real SPT::PTNSM(real k, real bl, real sigma_v, int a) const {
	return  pow2(bl)*factL(k,sigma_v,fl_spt,1.,0,a,6) *PLOOP(k,1) - 2.*factL(k,sigma_v,fl_spt,1.,1,a,6) *bl*PLOOP(k,2) + factL(k,sigma_v,fl_spt,1.,2,a,6) * PLOOP(k,3) + pow4(sig_8)*AB(k, bl, sigma_v, a);
    }

//qbias
real SPT::PTNSMq(real k, double barr[], real sigma_v, int a) const {
    double bk = barr[0]*sqrt((1+barr[1]*pow2(k))/(1+barr[2]*k));
  return  pow2(bk)*factL(k,sigma_v,fl_spt,1.,0,a,6) *PLOOP(k,1) - 2.*factL(k,sigma_v,fl_spt,1.,1,a,6) *bk*PLOOP(k,2) + factL(k,sigma_v,fl_spt,1.,2,a,6) * PLOOP(k,3) + pow4(sig_8)*AB(k, bk, sigma_v, a);
      }

/* Lag bias See Eq.23 1607.03150 for example*/
real SPT::PTNSMl(real k, double barr[], real sigma_v, int a) const {
  return  factL(k,sigma_v,fl_spt,1.,0,a,6)*(pow2(barr[0])*PLOOP(k,1)+Lag_bias(1,k,barr) + barr[2]) - 2.*factL(k,sigma_v,fl_spt,1.,1,a,6)*(barr[0]*PLOOP(k,2)+Lag_bias(2,k,barr)) + factL(k,sigma_v,fl_spt,1.,2,a,6)* PLOOP(k,3) + pow4(sig_8)*AB(k, barr[0], sigma_v, a);
      }

/*  nDGP TNS  Multipoles  */
// bl is linear bias
// a = 1 : Monopole
// a = 2 : Quadrupole
// a = 3 : Hexdecapole

real SPT::PTNSMnDGP(real k, real bl, real sigma_v, int a) const {
	return  pow2(bl)*factL(k,sigma_v,fl_spt,1.,0,a,6) *PLOOP(k,4) - 2.*factL(k,sigma_v,fl_spt,1.,1,a,6) *bl*PLOOP(k,5) + factL(k,sigma_v,fl_spt,1.,2,a,6) * PLOOP(k,6) + AB(k, bl, sigma_v, a);
}

real SPT::PTNSMnDGPq(real k, double barr[], real sigma_v, int a) const {
  double bk = barr[0]*sqrt((1+barr[1]*pow2(k))/(1+barr[2]*k));
	return  pow2(bk)*factL(k,sigma_v,fl_spt,1.,0,a,6) *PLOOP(k,4) - 2.*factL(k,sigma_v,fl_spt,1.,1,a,6) *bk*PLOOP(k,5) + factL(k,sigma_v,fl_spt,1.,2,a,6) * PLOOP(k,6) +  AB(k, bk, sigma_v, a);
}


/*  Arbitrary model TNS  Multipoles with DFoG term  */
//bl is linear bias
// a = 1 : Monopole
// a = 2 : Quadrupole
// a = 3 : Hexdecapole

real SPT::PTNSMmg(real k, double kmin, double kmax, int a,real bl, real sigma_v ) const {
	return  pow2(bl)*factL(k,sigma_v,fl_spt,1.,0,a,6)*PLOOPn(kmin,kmax,1,k)  - 2.*factL(k,sigma_v,fl_spt,1.,1,a,6) *bl*PLOOPn(kmin,kmax,2,k)  + factL(k,sigma_v,fl_spt,1.,2,a,6) * PLOOPn(kmin,kmax,3,k)+ pow4(sig_8)*ABn(kmin,kmax, bl, k, sigma_v, a);
}



/*Everything below are the functions used to construct the above spectra : */
//In order of appearance :
//1 LOOP CORRECTIONS : LCDM , nDGP, Numerical
// TNS A,B and C(nDGP only) correction terms - u integrated (Multipoles) and not


/*LCDM 1-loop terms */

/* LCDM P_{ab}^{(22)} */
real SPT::P22(real k, int a, int b) const {
    switch(a*b) {
        case 1:
            return P22_dd(k);
        case 2:
            return P22_dt(k);
        case 3:
            return P22_tt(k);
        default:
            warning("SPT: invalid indices, a = %d, b = %d\n", a,b);
            return 0;
    }
}

/* LCDM P_{ab}^{(13)} */
real SPT::P13(real k, int a, int b) const {
    switch(a*b) {
        case 1:
            return P13_dd(k);
        case 2:
            return P13_dt(k);
        case 3:
            return P13_tt(k);
        default:
        warning("SPT: invalid indices, a = %d, b = %d\n", a,b);
            return 0;
    }
}


/* P22- density, density */
static real f22_dd(const PowerSpectrum& P_L, real k, real r, real x) {
    real d = 1 + r*r - 2*r*x;
    if(d < 1e-5)
        return 0;
    else
        return  P_L(r*k)*P_L(k*sqrt(d)) * pow2(3*r + 7*x - 10*r*x*x) / (98.*pow2(d));
}

static real midintdd(const PowerSpectrum& P_L, real k, real r) {
    real KMAX = QMAXp/k;
    real KMIN = QMINp/k;
    real YMIN = Max(-1., (1.+r*r-KMAX*KMAX)/2./r);
    real YMAX = Min(0.9999999999999, (1.+r*r-KMIN*KMIN)/2./r);
return Integrate(bind(f22_dd, cref(P_L), k, r, _1), YMIN, YMAX ,  1e-4*P_L(k));
}

/* P_{\delta\delta}^{(22)} */
real SPT::P22_dd(real k) const {
 	real KMAX = QMAXp/k;
	real KMIN = QMINp/k;
    return k*k*k/(4*M_PI*M_PI) * pow4(Dl_spt/dnorm_spt) * Integrate<ExpSub>(bind(midintdd, cref(P_L), k, _1), KMIN, KMAX,epsrel);
}




/*P13 - DENSITY,DENSITY */
static real f13_dd(const PowerSpectrum& P_L, real k, real r) {
    real s;
    if(r < 1e-2)
        s = -168 + (928./5.)*pow2(r) - (4512./35.)*pow4(r) + (416./21.)*pow6(r);
    else if(fabs(r-1) < 1e-10)
        s = -88 + 8*(r-1);
    else if(r > 100)
        s = -488./5. + (96./5.)/pow2(r) - (160./21.)/pow4(r) - (1376./1155.)/pow6(r);
    else
        s = 12/pow2(r) - 158 + 100*pow2(r) - 42*pow4(r) + 3/pow3(r) * pow3(r*r - 1) * (7*r*r + 2) * log((1+r)/fabs(1-r));

    return P_L(k*r) * s;
}

/* P_{\delta\delta}^{(13)} */
real SPT::P13_dd(real k) const {
	real KMAX = QMAXp/k;
	real KMIN = QMINp/k;
    return k*k*k/(252.*4*M_PI*M_PI) *P_L(k) *  pow4(Dl_spt/dnorm_spt)* Integrate<ExpSub>(bind(f13_dd, cref(P_L), k, _1), KMIN, KMAX,epsrel);
}


/*P22 - DENSITY, VELOCITY */

static real f22_dt(const PowerSpectrum& P_L, real k, real r, real x) {
    real d = 1 + r*r - 2*r*x;
    if(d < 1e-5)
        return 0;
    else
        return P_L(k*r) *P_L(k*sqrt(d)) * (3*r + 7*x - 10*r*x*x)*(7*x - r - 6*r*x*x) / pow2(d);
}

static real midintdt(const PowerSpectrum& P_L, real k, real r) {
    real KMAX = QMAXp/k;
    real KMIN = QMINp/k;
    real YMIN = Max(-1., (1.+r*r-KMAX*KMAX)/2./r);
    real YMAX = Min(0.9999999999999, (1.+r*r-KMIN*KMIN)/2./r);
return Integrate(bind(f22_dt, cref(P_L), k, r, _1), YMIN, YMAX ,  1e-4*P_L(k));
}


/* P_{\delta\theta}^{(22)} */
real SPT::P22_dt(real k) const {
	real KMAX = QMAXp/k;
	real KMIN = QMINp/k;
    return  k*k*k/(98.*4.*M_PI*M_PI) *pow4(Dl_spt/dnorm_spt)*(-rempf*fl_spt)* Integrate<ExpSub>(bind(midintdt, cref(P_L), k, _1), KMIN, KMAX, epsrel);
}


/*P13- DENSITY,VELOCITY */
static real f13_dt(const PowerSpectrum& P_L, real k, real r) {
    real s;
    if(r < 1e-2)
        s = -168 + (416./5.)*pow2(r) - (2976./35.)*pow4(r) + (224./15.)*pow6(r);
    else if(fabs(r-1) < 1e-10)
        s = -152 - 56*(r-1);
    else if(r > 100)
        s = -200 + (2208./35.)/pow2(r) - (1312./105.)/pow4(r) - (1888./1155.)/pow6(r);
    else
        s = 24./pow2(r) - 202. + 56.*pow2(r) - 30.*pow4(r) + 3./pow3(r) * pow3(r*r - 1.) * (5.*r*r + 4.) * log((1.+r)/fabs(1.-r));

    return P_L(k*r) * s;
}


/* P_{\delta\theta}^{(13)} */
real SPT::P13_dt(real k) const {
   	real KMAX = QMAXp/k;
	  real KMIN = QMINp/k;
    return k*k*k/(252.*4.*M_PI*M_PI) *P_L(k) *(-rempf*fl_spt)*pow4(Dl_spt/dnorm_spt)* Integrate<ExpSub>(bind(f13_dt, cref(P_L), k, _1), KMIN, KMAX, epsrel);
}


/*P22- VELOCITY, VELOCITY */
static real f22_tt(const PowerSpectrum& P_L, real k, real r, real x) {
     real d = 1 + r*r - 2*r*x;
    if(d < 1e-5)
        return 0;
    else
        return P_L(k*r) * P_L(k*sqrt(d)) * pow2(7*x - r - 6*r*x*x) / pow2(d);
}


static real midintt(const PowerSpectrum& P_L, real k, real r) {
    real KMAX = QMAXp/k;
    real KMIN = QMINp/k;
    real YMIN = Max(-1., (1.+r*r-KMAX*KMAX)/2./r);
    real YMAX = Min(0.9999999999999, (1.+r*r-KMIN*KMIN)/2./r);
return Integrate(bind(f22_tt, cref(P_L), k, r, _1), YMIN, YMAX ,  1e-4*P_L(k));
}

 /* P_{\theta\theta}^{(22)} */
 real SPT::P22_tt(real k) const {
	real KMAX = QMAXp/k;
	real KMIN = QMINp/k;
	 return k*k*k/(98*4*M_PI*M_PI) *  pow4(Dl_spt/dnorm_spt)*pow2(rempf*fl_spt)*Integrate<ExpSub>(bind(midintt, cref(P_L), k, _1), KMIN, KMAX, epsrel);
}

/* P13- VELOCITY, VELOCITY*/
static real f13_tt(const PowerSpectrum& P_L, real k, real r) {
    real s;
    if(r < 1e-2)
        s = -56 - (32./5.)*pow2(r) - (96./7.)*pow4(r) + (352./105.)*pow6(r);
    else if(fabs(r-1) < 1e-10)

        s = -72 - 40*(r-1);

    else if(r > 100)
        s = -504./5. + (1248./35.)/pow2(r) - (608./105.)/pow4(r) - (160./231.)/pow6(r);
    else
        s = 12/pow2(r) - 82 + 4*pow2(r) - 6*pow4(r) + 3/pow3(r) * pow3(r*r - 1) * (r*r + 2) * log((1+r)/fabs(1-r));

    return P_L(k*r) * s;
}

/* P_{\theta\theta}^{(13)} */
real SPT::P13_tt(real k) const {
   	real KMAX = QMAXp/k;
	real KMIN = QMINp/k;
    return k*k*k/(84*4*M_PI*M_PI) *P_L(k) *pow4(Dl_spt/dnorm_spt)*pow2(rempf*fl_spt) * Integrate<ExpSub>(bind(f13_tt, cref(P_L), k, _1), KMIN, KMAX, epsrel);
}

/*DGP 1-LOOP TERMS */

/* DGP P_{ab}^{(22)} */
real SPT::P22D(real k, int a) const {
    switch(a) {
        case 1:
			return P22D_dd(k);
        case 2:
			return P22D_dt(k);
        case 3:
      return P22D_tt(k);
        default:
            warning("SPT: invalid indices, a = %d\n", a);
            return 0;
    }
}

/* DGP P_{ab}^{(13)} */
real SPT::P13D(real k, int a) const {
    switch(a) {
        case 1:
            return  P13D_dd(k);
        case 2:
            return  P13D_dt(k);
        case 3:
            return  P13D_tt(k);
        default:
            warning("SPT: invalid indices, a = %d\n", a);
            return 0;
    }
}

static real Df22_dd(const PowerSpectrum& P_L, real k, real r, real x) {
    real d = 1 + r*r - 2*r*x;
    if(d < 1e-5)
        return 0;
    else
        return P_L(r*k) * P_L(k*sqrt(d)) * 2 * r*r* (pow4(D_spt)*pow2((3*r + 7*x - 10*r*x*x)/(14.*r*d)) + pow2(F_spt) * pow2((x*x-1)/d) + pow2(D_spt)*F_spt*(3*r + 7*x - 10*r*x*x) * (1-x*x) / (7.*r*pow2(d)));}


static real Dmidintdd(const PowerSpectrum& P_L, real k, real r) {
    real KMAX = QMAXp/k;
    real KMIN = QMINp/k;
    real YMIN = Max(-1., (1.+r*r-KMAX*KMAX)/2./r);
    real YMAX = Min(0.9999999999999, (1.+r*r-KMIN*KMIN)/2./r);
return Integrate(bind(Df22_dd, cref(P_L), k, r, _1), YMIN, YMAX ,  1e-4*P_L(k));
}


/* P_{\delta\delta}^{(22)} */
real SPT::P22D_dd(real k) const {
 	real KMAX = QMAXp/k;
	real KMIN = QMINp/k;
    return k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt)* Integrate<ExpSub>(bind(Dmidintdd, cref(P_L), k, _1), KMIN, KMAX, epsrel);
}



static real Df13_dd(const PowerSpectrum& P_L, real k, real r) {
    real s;
	 if(r < 1e-2)
		s = pow4(D_spt)/252.*(-168 + (928./5.)*pow2(r) - (4512./35.)*pow4(r) + (416./21.)*pow6(r)) +  (D_spt*Cx_spt/6.) * (-96./5. * r*r + 96./35. * pow4(r) + 32./105.*pow6(r)) + (D_spt*I_spt/6.) * (32.*r*r - 96./5. * pow4(r) + 96./35. * pow6(r)) + (D_spt*KL_spt/12.) * ( 256./5. * r*r - 768./35. * pow4(r) + 256./105. * pow6(r)) ;

	else if(fabs(r-1) < 1e-10)

		s = pow4(D_spt)/252.*(-88 + 8*(r-1)) +  D_spt*Cx_spt/6. * (-16 - 24*(r-1)) + D_spt*I_spt/6. *(16 + 8*(r-1)) +  D_spt*KL_spt/12. *(32 + 32*(r-1)) ;

	else if(r>100)
		s = pow4(D_spt)/252.*(-488./5. + (96./5.)/pow2(r) - (160./21.)/pow4(r) - (1376./1155.)/pow6(r)) + D_spt*Cx_spt/6. * (-32 + 96./(5.*r*r) - 96./(35.*pow4(r)) - 32./(105.*pow6(r)))+ D_spt*I_spt/6. *(96./5. - 96./(35.*pow2(r)) - 32./(105.*pow4(r)) - 32./(385.*pow6(r))) + (D_spt*KL_spt/12.) * (256./5. - 768./(35.*pow2(r)) + 256./(105.*pow4(r)) + 256./(1155.*pow6(r))) ;

	else
   s = pow4(D_spt)/252.*(12./pow2(r) - 158. + 100.*pow2(r) - 42.*pow4(r) + 3./pow3(r) * pow3(r*r - 1) * (7*r*r + 2) * log((1+r)/fabs(1-r))) + D_spt*Cx_spt*(1./pow2(r) - 16./6. - pow2(r) + 1./(2.*pow3(r)) * pow3(r*r - 1) * log((1+r)/fabs(1-r))) + D_spt*I_spt*(1. + 16./6. * pow2(r) - pow4(r) + 1./(2.*r) * pow3(r*r - 1) * log((1+r)/fabs(1-r))) + D_spt*KL_spt/(12. * pow3(r)) * ( -6*r + 22*pow3(r) + 22*pow5(r) - 6*pow7(r) +  3*pow4(r*r - 1) * log((1+r)/fabs(1-r))) ;

    return P_L(k*r) * s;
}

/* P_{\delta\delta}^{(13)} */
real SPT::P13D_dd(real k) const {
	real KMAX = QMAXp/k;
	real KMIN = QMINp/k;
    return  k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt) *P_L(k) * Integrate<ExpSub>(bind(Df13_dd, cref(P_L), k, _1), KMIN, KMAX, epsrel);
}


static real Df22_dt(const PowerSpectrum& P_L, real k, real r, real x) {
    real d = 1 + r*r - 2*r*x;
    if(d < 1e-5)
        return 0;
    else
	 return  P_L(k*r) * P_L(k*sqrt(d)) * r * r * (fdgp_spt*pow4(D_spt)*(3*r + 7*x - 10*r*x*x)*(r - 7*x + 6*x*x*r) / (98.*r*r*pow2(d)) - F_spt*Fd_spt/H_spt * 2 * pow2(x*x-1) / pow2(d) +  Fd_spt*pow2(D_spt)/H_spt * (3*r + 7*x - 10*r*x*x) * (x*x-1) / (7.*r*pow2(d)) - F_spt*fdgp_spt*pow2(D_spt)*(r - 7*x + 6*r*x*x)*(x*x-1)/(7.*r*pow2(d)));
}


static real Dmidintdt(const PowerSpectrum& P_L, real k, real r) {
    real KMAX = QMAXp/k;
    real KMIN = QMINp/k;
    real YMIN = Max(-1., (1.+r*r-KMAX*KMAX)/2./r);
    real YMAX = Min(0.9999999999, (1.+r*r-KMIN*KMIN)/2./r);
return Integrate(bind(Df22_dt, cref(P_L), k, r, _1), YMIN, YMAX ,  1e-4*P_L(k));
}


/* P_{\delta\theta}^{(22)} */
real SPT::P22D_dt(real k) const {
        real KMAX = QMAXp/k;
    	real KMIN = QMINp/k;
    return k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt) * Integrate<ExpSub>(bind(Dmidintdt, cref(P_L), k, _1), KMIN, KMAX, epsrel);
}


static real Df13_dt(const PowerSpectrum& P_L, real k, real r) {
real s;
if(r < 1e-2)
s = -fdgp_spt*pow4(D_spt)/252.*(-168 + (416./5.)*pow2(r) - (2976./35.)*pow4(r) + (224./15.)*pow6(r)) - (D_spt*Cdx_spt + Cx_spt*fdgp_spt*D_spt*H_spt)/(12.*H_spt) * (-96./5. * r*r + 96./35. * pow4(r) + 32./105.*pow6(r)) - (fdgp_spt*D_spt*H_spt*I_spt + Id_spt*D_spt -Fd_spt*pow2(D_spt))/(H_spt*12.)* (32*r*r - 96./5. * pow4(r) + 96./35. * pow6(r)) - (D_spt*KLd_spt + D_spt*fdgp_spt*H_spt*KL_spt)/(H_spt*24.) * ( 256./5. * r*r - 768./35. * pow4(r) + 256./105. * pow6(r)) ;

else if(fabs(r-1) < 1e-10)

s = -fdgp_spt*pow4(D_spt)/252.*(-152 - 56*(r-1)) - (D_spt*Cdx_spt + Cx_spt*fdgp_spt*D_spt*H_spt)/(12*H_spt) * (-16 - 24*(r-1)) - (fdgp_spt*D_spt*H_spt*I_spt + Id_spt*D_spt -Fd_spt*pow2(D_spt))/(H_spt*12.)*(16 + 8*(r-1)) - (D_spt*KLd_spt + D_spt*fdgp_spt*H_spt*KL_spt)/(H_spt*24.) *(32 + 32*(r-1)) ;

else if(r>100)
s = -fdgp_spt*pow4(D_spt)/252.*(-200 + (2208./35.)/pow2(r) - (1312./105.)/pow4(r) - (1888./1155.)/pow6(r)) -  (D_spt*Cdx_spt + Cx_spt*fdgp_spt*D_spt*H_spt)/(12.*H_spt) *( -32 + 96./(5.*r*r) - 96./(35.*pow4(r)) - 32./(105.*pow6(r))) - (fdgp_spt*D_spt*H_spt*I_spt + Id_spt*D_spt -Fd_spt*pow2(D_spt))/(H_spt*12.)* (96./5. - 96./(35.*pow2(r)) - 32./(105.*pow4(r)) - 32./(385.*pow6(r))) - (D_spt*KLd_spt + D_spt*fdgp_spt*H_spt*KL_spt)/(12.*H_spt) * (256./5. - 768./(35.*pow2(r)) + 256./(105.*pow4(r)) + 256./(1155.*pow6(r))) ;

else
s = -fdgp_spt*pow4(D_spt)/252.*(24/pow2(r) - 202 + 56*pow2(r) - 30*pow4(r) + 3/pow3(r) * pow3(r*r - 1) * (5*r*r + 4) * log((1+r)/fabs(1-r))) - (D_spt*Cdx_spt + Cx_spt*fdgp_spt*D_spt*H_spt)/(H_spt*12.*pow3(r))*(6*r - 16*pow3(r) - 6*pow5(r) + 3*pow3(r*r-1)*log((1+r)/fabs(1-r))) - (fdgp_spt*D_spt*H_spt*I_spt + Id_spt*D_spt -Fd_spt*pow2(D_spt))/(H_spt*12.*r)*(6*r + 16*pow3(r) - 6*pow5(r) + 3*pow3(r*r-1)*log((1+r)/fabs(1-r))) - (D_spt*KLd_spt + D_spt*fdgp_spt*H_spt*KL_spt)/(H_spt*24.*pow3(r))*(-6*r + 22*pow3(r) + 22*pow5(r) - 6*pow7(r) + 3*pow4(r*r-1)*log((1+r)/fabs(1-r)));

return  P_L(k*r) * s;
}

/* P_{\delta\theta}^{(13)} */
real SPT::P13D_dt(real k) const {
   	real KMAX = QMAXp/k;
	  real KMIN = QMINp/k;
    return k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt) *P_L(k) * Integrate<ExpSub>(bind(Df13_dt, cref(P_L), k, _1), KMIN, KMAX,  epsrel);
}


static real Df22_tt(const PowerSpectrum& P_L, real k, real r, real x) {
    real d = 1 + r*r - 2*r*x;
    if(d < 1e-5)
        return 0;
    else
        return  P_L(k*r) * P_L(k*sqrt(d)) * r * r * 2 * (pow2(fdgp_spt*pow2(D_spt)) * pow2((7*x - r - 6*r*x*x)/(14*r*d)) + pow2(Fd_spt/H_spt)*pow2((x*x-1)/d) + (fdgp_spt*pow2(D_spt)*Fd_spt/H_spt)*(r-7*x+6*x*x*r)*(x*x-1)/(7*r*pow2(d)));}


static real Dmidintt(const PowerSpectrum& P_L, real k, real r) {
    real KMAX = QMAXp/k;
    real KMIN = QMINp/k;
    real YMIN = Max(-1., (1.+r*r-KMAX*KMAX)/2./r);
    real YMAX = Min(0.9999999999999, (1.+r*r-KMIN*KMIN)/2./r);
return Integrate(bind(Df22_tt, cref(P_L), k, r, _1), YMIN, YMAX ,  1e-4*P_L(k));
}

/* P_{\theta\theta}^{(22)} */
real SPT::P22D_tt(real k) const {
	real KMAX = QMAXp/k;
	real KMIN = QMINp/k;
    return k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt) * Integrate<ExpSub>(bind(Dmidintt, cref(P_L), k, _1), KMIN, KMAX, epsrel);
}


static real Df13_tt(const PowerSpectrum& P_L, real k, real r) {
real s;

if(r < 1e-2)

s =  pow2(fdgp_spt*pow2(D_spt))/84.*(-56 - (32./5.)*pow2(r) - (96./7.)*pow4(r) + (352./105.)*pow6(r)) +  fdgp_spt*D_spt*Cdx_spt/(H_spt*6.) * (-96./5. * r*r + 96./35. * pow4(r) + 32./105.*pow6(r)) + fdgp_spt*D_spt*(Id_spt-Fd_spt*D_spt)/(6.*H_spt) * (32*r*r - 96./5. * pow4(r) + 96./35. * pow6(r)) + fdgp_spt*D_spt*KLd_spt/(H_spt*12.) * ( 256./5.*r*r - 768./35. * pow4(r) + 256./105. * pow6(r)) ;

else if(fabs(r-1) < 1e-10)

s =  pow2(fdgp_spt*pow2(D_spt))/84.*(-72 - 40*(r-1)) +   fdgp_spt*D_spt*Cdx_spt/(H_spt*6.) * (-16 - 24*(r-1)) + fdgp_spt*D_spt*(Id_spt-Fd_spt*D_spt)/(6.*H_spt) *(16 + 8*(r-1)) +  fdgp_spt*D_spt*KLd_spt/(H_spt*12.) *(32 + 32*(r-1)) ;

else if(r>100)
s =  pow2(fdgp_spt*pow2(D_spt))/84.*(-504./5. + (1248./35.)/pow2(r) - (608./105.)/pow4(r) - (160./231.)/pow6(r)) + fdgp_spt*D_spt*Cdx_spt/(H_spt*6.)*(-32 + 96./(5.*r*r) - 96./(35.*pow4(r)) - 32./(105.*pow6(r)))+ fdgp_spt*D_spt*(Id_spt-Fd_spt*D_spt)/(6.*H_spt) *(96./5. - 96./(35.*pow2(r)) - 32./(105.*pow4(r)) - 32./(385.*pow6(r)) ) + fdgp_spt*D_spt*KLd_spt/(H_spt*12.)  * (256./5. - 768./(35.*pow2(r)) + 256./(105.*pow4(r)) + 256./(1155.*pow6(r))) ;

else
s =  pow2(fdgp_spt*pow2(D_spt))/84.*(12./pow2(r) - 82 + 4*pow2(r) - 6*pow4(r) + 3/pow3(r) * pow3(r*r - 1) * (r*r + 2) * log((1+r)/fabs(1-r))) + fdgp_spt*D_spt*Cdx_spt/H_spt * 1./(6.*pow3(r))*(6*r - 16*pow3(r) - 6*pow5(r) + 3*pow3(r*r-1)* log((1+r)/fabs(1-r))) + fdgp_spt*D_spt*(Id_spt-Fd_spt*D_spt)/(H_spt*6.*r)*(6*r + 16*pow3(r) - 6*pow5(r) + 3*pow3(r*r-1)*log((1+r)/fabs(1-r))) + fdgp_spt*D_spt*(KLd_spt)/(H_spt*12.*pow3(r))*(-6*r + 22*pow3(r) + 22*pow5(r) - 6*pow7(r) + 3*pow4(r*r-1)* log((1+r)/fabs(1-r))) ;

return P_L(k*r) * s;}

/* P_{\theta\theta}^{(13)} */

real SPT::P13D_tt(real k) const {
   	real KMAX = QMAXp/k;
	real KMIN = QMINp/k;
    return k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt) *P_L(k) * Integrate<ExpSub>(bind(Df13_tt, cref(P_L), k, _1), KMIN, KMAX, epsrel);
}


/* Higher order bias terms for P_dd. See Eq.24 and 25 of 1607.03150 for example*/

// P_b2
static real lbias_selec(const PowerSpectrum& P_L, int a,  real bias[], real k, real r, real x) {
    real d = 1 + r*r - 2*r*x;
    if(d < 1e-5)
        return 0;
    else{
      double d1 = sqrt(d);
      double sker = 2./3.-(1.-pow2((x-r)/d1));
      double dker = -2./7.*(1.-x*x);
      double plkr = P_L(r*k);
      double bs2 = -4./7.*(bias[0]-1.);
      double b3n = 32./315.*(bias[0]-1.);
      double G2ker = G2eds(k*r,k*d1,(x-r)/d1);
      double F2ker = F2eds(k*r,k*d1,(x-r)/d1);
      switch (a) {
      case 1:
      return  plkr*P_L(k*d1)*pow2(r)*
              (2.*bias[0]*bias[1]*F2ker // b2
              +2.*bias[0]*bs2*F2ker*sker //bs2
              +pow2(bias[1])/2. //b22 a
              +bias[1]*bs2*sker //b2s2 a
              +pow2(bs2)*pow2(sker)/2.) //bs22 a

            - pow2(r*plkr)*(pow2(bias[1])/2. //b22 b
                            +2.*bias[1]*bs2/3. // b2s2 b
                            +pow2(bs2)*2./9.) // bs22 b

            +P_L(k)*b3n*bias[0]*plkr*2.*105./16.*pow2(r)*(dker*sker+8./63.); // b3 terms

     case 2 :
     return   -rempf*fl_spt*plkr*P_L(k*d1)*pow2(r)*(bias[1]*G2ker //b2t
                                                 +bs2*G2ker*sker) // bs2t
              -rempf*fl_spt*P_L(k)*b3n*plkr*105./16.*pow2(r)*(dker*sker+8./63.);
            }
}
}

static real midpb2(const PowerSpectrum& P_L, int a, real bias[], real k, real r) {
    real KMAX = QMAXp/k;
    real KMIN = QMINp/k;
    real YMIN = Max(-1., (1.+r*r-KMAX*KMAX)/2./r);
    real YMAX = Min(0.9999999999999, (1.+r*r-KMIN*KMIN)/2./r);
return Integrate(bind(lbias_selec, cref(P_L), a, bias, k, r, _1), YMIN, YMAX , 1e-3);
}

//P_bs2
real SPT::Lag_bias(int a, real k, real bias[]) const {
 	real KMAX = QMAXp/k;
	real KMIN = QMINp/k;
    return k*k*k/(4*M_PI*M_PI) * pow4(sig_8*Dl_spt/dnorm_spt) * Integrate<ExpSub>(bind(midpb2, cref(P_L), a, bias, k, _1), KMIN, KMAX,epsrel);
}


// NUMERICAL KERNEL 1-LOOP TERMS
//Integrating over angle - P22 numerical:
static double F2F(int y, int a,  const PowerSpectrum& P_L, double k, double r){
	double temp_ps;
	double myresult=0.;
	switch (a) {
			case 1:
			for( int i = 0; i < n1; i++ )
				{
        double d = 1+ r*r - 2*r*x128[i];
		 		if(d < 1e-5){
			 	temp_ps=0.;
		 		}
		  	else {
		  	temp_ps = P_L(k*sqrt(d)) * pow2(F2_nk[i*n2 + y]);
				}
        myresult += w128[i] * temp_ps;
				}
        return  2. * r * r * myresult;
        break;
			case 2:
      myresult = 0.;
			for( int i = 0; i < n1; i++ )
				{
        double d = 1+ r*r - 2*r*x128[i];
				if(d < 1e-5){
				temp_ps=0.;
				}
				else {
				temp_ps = P_L(k*sqrt(d)) * G2_nk[i*n2 + y]*F2_nk[i*n2 + y];
				}
        myresult += w128[i] * temp_ps;
				}
        return  2. * r * r * myresult;
        break;
			case 3:
			for( int i = 0; i < n1; i++ )
				{
				double d = 1+ r*r - 2*r*x128[i];
				if(d < 1e-5){
				temp_ps=0.;
				}
				else {
				temp_ps = P_L(k*sqrt(d)) * pow2(G2_nk[i*n2 + y]);
				}
        myresult += w128[i] * temp_ps;
				}
        return 2. * r * r * myresult;
        break;
			}
	}

  //Integrating over angle -P13 numerical:
  static double F3F(int y, int a,  const PowerSpectrum& P_L, double k, double r){
    double myresult = 0.;
    // Set the limits of angular integrations
  	switch (a) {
  			case 1:
  			 for( int i = 0; i < n1; i++ )
  				{
  				myresult += w128[i] * F1_nk * F3_nk[i*n2 + y];
  				}
        return   6. * r * r * myresult;
        break;
  			case 2:
  			 for( int i = 0; i < n1; i++ )
  				{
  					myresult += w128[i] * (G1_nk*F3_nk[i*n2 + y] + F1_nk*G3_nk[i*n2 + y]);
  				}
        return    3. * r * r * myresult;
        break;
  			case 3:
  				for( int i = 0; i < n1; i++ )
  				{
  					myresult += w128[i] * G1_nk * G3_nk[i*n2 + y];
  				}
          return   6. * r * r * myresult;
          break;

  		}
    }

/*Selection function for numerical kernel magnitude */

inline int num_selec_mag(double kmin, double kmax, double y){
      double YMIN =QMINp/kmax;
      double YMAX = QMAXp/kmin;
    	int mag_int;
  //  	mag_int = (y-YMIN)*(n2*1.-1)/(YMAX-YMIN); // linear
    //mag_int = sqrt((y-QMINp/kmax)*kmin/QMAXp)*(n2-1); //quadratic
    mag_int = (int)round((n2-1)*log(y/YMIN)/log(YMAX/YMIN)); //exponential
    		return mag_int;
}

/* P^{(22)} */
// a = 1 : delta delta
// a = 2 : delta theta
// a = 3 : theta theta
real SPT::P22n(double kmin, double kmax,  int a, real k) const {
    	real KMAX = QMAXp/k;
    	real KMIN = QMINp/k;
      int y1 = num_selec_mag(kmin, kmax, KMIN);
      int y2 = num_selec_mag(kmin, kmax, KMAX);
      double y[n2];
      double integrand[n2];
      for (int i = y1; i<=y2; i++){
      y[i] = QMINp/kmax * exp(i*log(QMAXp*kmax/(QMINp*kmin))/(n2*1.-1.)); // exponential sampling
  //    y[i] = i*1./(n2-1.)*(QMAXp/kmin-QMINp/kmax)+QMINp/kmax; // linear sampling
      integrand[i] = P_L(k*y[i]) * F2F(i, a, cref(P_L), k, y[i]);
}
double res = 0.;
  for( int i = y1+1; i <= y2; ++ i ){
res += 0.5 * (y[i] - y[i-1])*(integrand[i] + integrand[i-1]);
}
return  k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt) * res;
}

real SPT::P13n(double kmin, double kmax,  int a, real k) const {
  real KMAX = QMAXp/k;
  real KMIN = QMINp/k;
  int y1 = num_selec_mag(kmin, kmax, KMIN);
  int y2 = num_selec_mag(kmin, kmax, KMAX);
  double y[n2];
  double integrand[n2];
  for (int i = y1; i<=y2; i++){
  y[i] = QMINp/kmax * exp(i*log(QMAXp*kmax/(QMINp*kmin))/(n2*1.-1.));
//  y[i] = i*1./(n2-1.)*(QMAXp/kmin-QMINp/kmax)+QMINp/kmax; // linear sampling
  integrand[i] = P_L(k*y[i]) * F3F(i, a, cref(P_L), k, y[i]);
  }
double res = 0.;
  for( int i = y1+1; i <= y2; ++ i ){
res += 0.5 * (y[i] - y[i-1])*(integrand[i] + integrand[i-1]);
}
  return  k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt) *  P_L(k) * res;
}

/// REVAMP P22n and P13n integration
static double ploopn2_mgdd( const PowerSpectrum& P_L, double vars[], double k, double r, double x){
  double kargs[4],kv[3],xv[3], p22,p13,d;
  IOW iow;
        kv[0] = k;
        kv[1] = k*r;
        kv[2] = kv[1];
        xv[0] = -0.9999999999999999; // critical number for de solving - singularity encountered otherwise
        xv[1] = x;
        xv[2] = -x;
        kargs[0] = sqrt(kv[1]*kv[1]+kv[0]*kv[0]-2.*kv[1]*kv[0]*xv[1]);
        kargs[2] = sqrt(kv[1]*kv[1]+2.*kv[1]*kv[2]*xv[0]+kv[2]*kv[2]);
        kargs[1] = sqrt(kv[2]*kv[2]+2.*kv[2]*kv[0]*xv[2]+kv[0]*kv[0]);
        kargs[3] = sqrt(kv[0]*kv[0]+2.*kv[0]*kv[1]*xv[1]+kv[1]*kv[1]);
        iow.initn2(vars[4],kv,xv,kargs,vars[0],vars[1],vars[2],vars[3]);
        p22 = pow2(F2_nk[0]);
        p13 = F1_nk*F3_nk[0];

    return pow2(r)*2.*P_L(k*r)*(P_L(kargs[0])*p22 + 3.*P_L(k)*p13);
}
static double ploopn2_mgdt( const PowerSpectrum& P_L, double vars[], double k, double r, double x){
  double kargs[4],kv[3],xv[3], p22,p13,d;
  IOW iow;
        kv[0] = k;
        kv[1] = k*r;
        kv[2] = kv[1];
        xv[0] = -0.99999999;
        xv[1] = x;
        xv[2] = -x;
        kargs[0] = sqrt(kv[1]*kv[1]+kv[0]*kv[0]-2.*kv[1]*kv[0]*xv[1]);
        kargs[2] = sqrt(kv[1]*kv[1]+2.*kv[1]*kv[2]*xv[0]+kv[2]*kv[2]);
        kargs[1] = sqrt(kv[2]*kv[2]+2.*kv[2]*kv[0]*xv[2]+kv[0]*kv[0]);
        kargs[3] = sqrt(kv[0]*kv[0]+2.*kv[0]*kv[1]*xv[1]+kv[1]*kv[1]);
        iow.initn2(vars[4],kv,xv,kargs,vars[0],vars[1],vars[2],vars[3]);
        p22 = G2_nk[0]*F2_nk[0];
        p13 = 0.5*(G1_nk*F3_nk[0] + F1_nk*G3_nk[0]);
    return pow2(r)*2.*P_L(k*r)*(P_L(kargs[0])*p22 + 3.*P_L(k)*p13);
}

static double ploopn2_mgtt( const PowerSpectrum& P_L, double vars[], double k, double r, double x){
  double kargs[4],kv[3],xv[3], p22,p13,d;
  IOW iow;
        kv[0] = k;
        kv[1] = k*r;
        kv[2] = kv[1];
        xv[0] = -0.99999999;
        xv[1] = x;
        xv[2] = -x;
        kargs[0] = sqrt(kv[1]*kv[1]+kv[0]*kv[0]-2.*kv[1]*kv[0]*xv[1]);
        kargs[2] = sqrt(kv[1]*kv[1]+2.*kv[1]*kv[2]*xv[0]+kv[2]*kv[2]);
        kargs[1] = sqrt(kv[2]*kv[2]+2.*kv[2]*kv[0]*xv[2]+kv[0]*kv[0]);
        kargs[3] = sqrt(kv[0]*kv[0]+2.*kv[0]*kv[1]*xv[1]+kv[1]*kv[1]);
        iow.initn2(vars[4],kv,xv,kargs,vars[0],vars[1],vars[2],vars[3]);
        p22 = pow2(G2_nk[0]);
        p13 = G1_nk * G3_nk[0];
    return pow2(r)*2.*P_L(k*r)*(P_L(kargs[0])*p22 + 3.*P_L(k)*p13);
}
/// .... for IDE
static double ploopn2_idedd( const PowerSpectrum& P_L, double vars[], double k, double r, double x){
  double kargs[4],kv[3],xv[3], p22,p13,d;
  IOW iow;
        int p4 = vars[5];
        kv[0] = k;
        kv[1] = k*r;
        kv[2] = kv[1];
        xv[0] = -0.99999999;
        xv[1] = x;
        xv[2] = -x;
        kargs[0] = sqrt(kv[1]*kv[1]+kv[0]*kv[0]-2.*kv[1]*kv[0]*xv[1]);
        kargs[2] = sqrt(kv[1]*kv[1]+2.*kv[1]*kv[2]*xv[0]+kv[2]*kv[2]);
        kargs[1] = sqrt(kv[2]*kv[2]+2.*kv[2]*kv[0]*xv[2]+kv[0]*kv[0]);
        kargs[3] = sqrt(kv[0]*kv[0]+2.*kv[0]*kv[1]*xv[1]+kv[1]*kv[1]);
        iow.initn2_ide(vars[4],kv,xv,kargs,vars[0],vars[1],vars[2],vars[3],p4);
        p22 = pow2(F2_nk[0]);
        p13 = F1_nk*F3_nk[0];
    return pow2(r)*2.*P_L(k*r)*(P_L(kargs[0])*p22 + 3.*P_L(k)*p13);
}
static double ploopn2_idedt( const PowerSpectrum& P_L, double vars[], double k, double r, double x){
  double kargs[4],kv[3],xv[3], p22,p13,d;
  IOW iow;
        int p4 = vars[5];
        kv[0] = k;
        kv[1] = k*r;
        kv[2] = kv[1];
        xv[0] = -0.99999999;
        xv[1] = x;
        xv[2] = -x;
        kargs[0] = sqrt(kv[1]*kv[1]+kv[0]*kv[0]-2.*kv[1]*kv[0]*xv[1]);
        kargs[2] = sqrt(kv[1]*kv[1]+2.*kv[1]*kv[2]*xv[0]+kv[2]*kv[2]);
        kargs[1] = sqrt(kv[2]*kv[2]+2.*kv[2]*kv[0]*xv[2]+kv[0]*kv[0]);
        kargs[3] = sqrt(kv[0]*kv[0]+2.*kv[0]*kv[1]*xv[1]+kv[1]*kv[1]);
        iow.initn2_ide(vars[4],kv,xv,kargs,vars[0],vars[1],vars[2],vars[3],p4);
        p22 = G2_nk[0]*F2_nk[0];
        p13 = 0.5*(G1_nk*F3_nk[0] + F1_nk*G3_nk[0]);
    return pow2(r)*2.*P_L(k*r)*(P_L(kargs[0])*p22 + 3.*P_L(k)*p13);
}

static double ploopn2_idett( const PowerSpectrum& P_L, double vars[], double k, double r, double x){
  double kargs[4],kv[3],xv[3], p22,p13,d;
  IOW iow;
        int p4 = vars[5];
        kv[0] = k;
        kv[1] = k*r;
        kv[2] = kv[1];
        xv[0] = -0.99999999;
        xv[1] = x;
        xv[2] = -x;
        kargs[0] = sqrt(kv[1]*kv[1]+kv[0]*kv[0]-2.*kv[1]*kv[0]*xv[1]);
        kargs[2] = sqrt(kv[1]*kv[1]+2.*kv[1]*kv[2]*xv[0]+kv[2]*kv[2]);
        kargs[1] = sqrt(kv[2]*kv[2]+2.*kv[2]*kv[0]*xv[2]+kv[0]*kv[0]);
        kargs[3] = sqrt(kv[0]*kv[0]+2.*kv[0]*kv[1]*xv[1]+kv[1]*kv[1]);
        iow.initn2_ide(vars[4],kv,xv,kargs,vars[0],vars[1],vars[2],vars[3],p4);
        p22 = pow2(G2_nk[0]);
        p13 = G1_nk * G3_nk[0];
    return pow2(r)*2.*P_L(k*r)*(P_L(kargs[0])*p22 + 3.*P_L(k)*p13);
}

// Choose a {1,...6}: P_dd,P_dt, P_tt (MG) ,P_dd,P_dt, P_tt (IDE)

double SPT::PLOOPn2(int a, double vars[], double k, double err) const{
  IOW iow;
double loop, tree;
double prefac = k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt);
real KMAX = QMAXp/k;
real KMIN = QMINp/k;
double c[2] = {KMIN,-1.};
double d[2] = {KMAX, 1.};
switch (a) {
  case 0:
  //  loop = prefac*Integrate<2>(bind(ploopn2_mgdd,cref(P_L),vars,k,_1,_2), c, d, err,1e-2);
    iow.initn_lin(0, vars[4], k, vars[0],vars[1], 1.,1.);
    tree = pow2(F1_nk/dnorm_spt)*P_L(k);
    return tree;
  case 1:
    loop = prefac*Integrate<2>(bind(ploopn2_mgdd,cref(P_L),vars,k,_1,_2), c, d, err,1e-4);
    tree = pow2(F1_nk/dnorm_spt)*P_L(k);
    return loop+tree;
  case 2:
    loop = prefac*Integrate<2>(bind(ploopn2_mgdt,cref(P_L),vars,k,_1,_2), c, d, err,1e-2);
    tree = (G1_nk*F1_nk)/pow2(dnorm_spt)*P_L(k);
    return loop+tree;
  case 3:
    loop = prefac*Integrate<2>(bind(ploopn2_mgtt,cref(P_L),vars,k,_1,_2), c, d, err,1e-2);
    tree = pow2(G1_nk/dnorm_spt)*P_L(k);
    return loop+tree;
  case 4:
    iow.initnorm_ide(vars[4],vars[0],vars[1], vars[2],vars[3],1);
    tree = pow2(D_spt/dnorm_spt)*P_L(k);
    return tree;
  case 5:
    iow.initnorm_ide(vars[4],vars[0],vars[1], vars[2],vars[3],1);
    loop = k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt)*Integrate<2>(bind(ploopn2_idedd,cref(P_L),vars,k,_1,_2), c, d, err,1e-2);
    tree = pow2(F1_nk/dnorm_spt)*P_L(k);
    return loop+tree;
  case 6:
    iow.initnorm_ide(vars[4],vars[0],vars[1], vars[2],vars[3],1);
    loop = k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt)*Integrate<2>(bind(ploopn2_idedt,cref(P_L),vars,k,_1,_2), c, d, err,1e-2);
    tree = (G1_nk*F1_nk)/pow2(dnorm_spt)*P_L(k);
    return loop+tree;
  case 7:
    iow.initnorm_ide(vars[4],vars[0],vars[1], vars[2],vars[3],1);
    loop = k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt)*Integrate<2>(bind(ploopn2_idett,cref(P_L),vars,k,_1,_2), c, d, err,1e-2);
    tree = pow2(G1_nk/dnorm_spt)*P_L(k);
    return loop+tree;
}}

/// Should extend this ^ to TNS : AB term needs special initialisation because of x-integral limits

////////////////////////////////////////////////
///////HALO-FIT FOR NON-SCALE DEP MODELS : takashi et al 2012 fits  1208.2701///////////
///////////////////////////////////////////////

double wintcambs[3];
static double wintcamb(const PowerSpectrum& P_L, double r){
 int n1 = 3000;
 double sum1 = 0.;
 double sum2 = 0.;
 double sum3 = 0.;
 double anorm = 1./(2.*pow2(M_PI));
 double t,k, pkterm,x,x2,w1,w2,w3,fac;
 for(int i = 1; i<=n1 ; i++){
   t = (i-0.5)/n1;
   k = -1. + 1./t;
   pkterm = pow2(D_spt/dnorm_spt)*P_L(k)*pow3(k)*anorm;
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
wintcambs[0] = sqrt(sum1);
wintcambs[1] = -sum2/sum1;
wintcambs[2] = -sum2*sum2/sum1/sum1 - sum3/sum1;
}

double phpars[13];
void parscamb(const PowerSpectrum& P_L)
{
  double xlogr1 = -2.;
  double xlogr2 = 3.5;
  double rmid = 0.;
  double diff = 1.;
  for (;;){
    rmid = (xlogr2 + xlogr1)/2.;
    rmid = pow(10,rmid);
    wintcamb(cref(P_L),rmid);
    diff = wintcambs[0]-1.;
    if (fabs(diff)<=0.0001) {
      phpars[0] = 1./rmid;
      phpars[11] = -3.-wintcambs[1];
      phpars[12] = -wintcambs[2];
      break;
    }
   else if(diff>0.0001){
     xlogr1 = log10(rmid);
   }
   else if (diff < -0.0001){
     xlogr2 = log10(rmid);
   }
    if (xlogr2<-1.9999){
      phpars[0] = 1000.;
      break;
    }
    else if(xlogr2>3.4999){
      phpars[0] = 1000.;
      break;
    }
  }
}

void SPT::phinit(double scalef, double omega0)const{
    parscamb(cref(P_L));

    double neff = phpars[11];
    double curv = phpars[12];

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

     phpars[1] = pow(10.,ane); // an
     phpars[2] = pow(10.,bne); // bn
     phpars[3] = pow(10.,cne); //cn

     phpars[4] = 0.1971 - 0.0843*neff + 0.8460*curv; //gan
     phpars[5] = fabs(6.0835 + 1.3373*neff - 0.1959*neff2- 5.5274*curv); // alpha
     phpars[6] = 2.0379 - 0.7354*neff + 0.3157*neff2 + 1.2490*neff3 + 0.3980*neff4- 0.1682*curv; //beta

     phpars[7] = pow(10.,nue); // nun

     if (fabs(1-omgm) > 0.01) {
       double f1a = pow(omgm,-0.0732);
       double f2a = pow(omgm,-0.1423);
       double f3a = pow(omgm,0.0725);

       double f1b = pow(omgm,-0.0307);
       double f2b = pow(omgm,-0.0585);
       double f3b = pow(omgm,0.0743);

        double frac = omgl/(1. - omgm);
       phpars[8] = frac*f1b + (1.-frac)*f1a;
       phpars[9] = frac*f2b + (1.-frac)*f2a;
       phpars[10] = frac*f3b + (1.-frac)*f3a;

        }
      else{
        phpars[8]=1.;
        phpars[9]=1.;
        phpars[10]=1.;
      }

  }

double SPT::PHALO(double k) const{
  if (phpars[0] == 1000. || k<=0.005) {
    return pow2(D_spt/dnorm_spt)*P_L(k);
  }
  else{
  double kcub = pow3(k)/2./pow2(M_PI);
  double deltahp = phpars[1]*pow(k/phpars[0],phpars[8]*3.)/(1.+phpars[2]*pow(k/phpars[0],phpars[9])+pow(phpars[3]*phpars[10]*k/phpars[0],3.-phpars[4]));
  double deltah = deltahp/(1.+ phpars[7]/pow2(k/phpars[0]));
  double deltaplin = kcub*pow2(D_spt/dnorm_spt)*P_L(k);
  double deltaq = deltaplin*pow(1.+ deltaplin,phpars[6])/(1.+phpars[5]*deltaplin)*exp(-(k/phpars[0])/4.-pow2(k/phpars[0])/8.);
  double deltanl = deltaq + deltah;

  return  deltanl/kcub;
}
}
// MG HALOFIT FOR HU-SAWICKI MODEL : 1312.1291 ///
// double SPT::MGPHALO(double k, double mg1) const{
//   double kcub = pow3(k)/2./pow2(M_PI);
//   double deltahp = phpars[1]*pow(k/phpars[0],phpars[8]*3.)/(1.+phpars[2]*pow(k/phpars[0],phpars[9])+pow(phpars[3]*phpars[10]*k/phpars[0],3.-phpars[4]));
//
//   double deltah = deltahp/(1.+ phpars[7]/pow2(k/phpars[0]));
//
//   double amg = phpars[1] + mg1*(1.943697 + 7.776061*phpars[11] + 3.186278*pow2(phpars[11]) + 6.916149*phpars[12]);
//   double bmg = phpars[2] + mg1*(0.999088 + 8.480852*phpars[11] + 3.644990*pow2(phpars[11]) + 9.519407*phpars[12]);
//   double cmg = phpars[3] + mg1*(1.934338 + 2.511626*phpars[11] + 0.792323*pow2(phpars[11]) + 0.337545*phpars[12]);
//   double gmg = phpars[4] + mg1*(2.971117 - 1.702803*phpars[11] - 1.284630*pow2(phpars[11]) - 6.797889*phpars[12]);
//   double deltahppmg = amg*pow(k/phpars[0],phpars[8]*3.)/(1.+ bmg*pow(k/phpars[0],phpars[9])+pow(cmg*phpars[10]*k/phpars[0],3.-gmg));
//
//   double dfac =
//   double ximg = exp(dfac*(-10.6564 - 0.995708*phpars[11] + 1.169303*pow2(phpars[11]) + 17.519593*phpars[12]))
//   double mumg = mg1*(1.440371 + 1.819927*phpars[11] + 0.564780*pow2(phpars[11]) + 0.274286*phpars[12]);
//   double numg =  phpars[7] + mg1*(-2282.5327 - 2135.1213*phpars[11] - 2258.1919*pow2(phpars[11]) - 2378.1342*phpars[12]);
//
//   double deltah = deltahpmg*ximg/(1.+ mumg/(k/phpars[0]) +numg/pow2(k/phpars[0]));
//
//
//   double deltaplin = kcub*pow2(D_spt/dnorm_spt)*P_L(k);
//   double deltaplinmg = deltaplin*(1.+ mg1*(-0.832105 - 0.238632*phpars[11] + 0.427827*phpars[12]));
//   double betamg = phpars[6] + mg1*(-0.318559 + 2.963588*phpars[11] + 1.551244*pow2(phpars[11]) + 1.150983*phpars[12]);
//   double alphamg = phpars[5] + mg1*(-3.367256 + 3.888473*phpars[11] + 2.294713*pow2(phpars[11]) + 8.821165*phpars[12])
//   double deltaq = deltaplin*pow(1.+ deltaplinmg,betamg)/(1.+alphamg*deltaplinmg)*exp(-(k/phpars[0])/4.-pow2(k/phpars[0])/8.);
//
//   double deltanl = deltaq + deltah;
//
//   return  deltanl/kcub;
// }



////////////////////////////////////////////////
///////REDSHIFT SPACE MODEL-TNS TERMS///////////
///////////////////////////////////////////////

/* Non vanishing A, B and C terms */

// Notes for analytic ABC terms : DEFINITION OF THETA= - DEL V / aHf (and is accounted for in definition of Cross Bispectrum)
// See http://arxiv.org/pdf/1006.0699v1.pdf for derivation

// NOTE ON PARAMETERS: x =k1.k/(k^2*r)


static real ABC(int a, const PowerSpectrum& P_L, real k, real r, real x) {
    real d = 1. + r*r - 2.*r*x;
    if(d < 1e-5)
        return 0;
    else
		switch(a) {
		// A terms
				//A11
			case 1:
				return (-r*r*r/7.)*(x+6.*x*x*x+r*r*x*(-3.+10.*x*x)+r*(-3.+x*x-12.*x*x*x*x))*P_L(k) * P_L(k*sqrt(d)) / pow2(d);
        break;
      	//A12
			case 2:
				return  (r*r*r*r/14.)*(x*x-1.)*(-1.+7.*r*x-6.*x*x) * P_L(k) * P_L(k*sqrt(d)) / pow2(d);
        break;
      	//A22
			case 3:
				return  (r*r*r/14.)*(r*r*x*(13.-41.*x*x)-4.*(x+6.*x*x*x)+r*(5.+9.*x*x+42.*x*x*x*x))*P_L(k)*P_L(k*sqrt(d))/ pow2(d);
        break;
      	//A33
			case 4:
				return (r*r*r/14.)*(1.-7.*r*x+6.*x*x)*(-2.*x+r*(-1.+3.*x*x))*P_L(k)*P_L(k*sqrt(d))/ pow2(d);
        break;
      	//A11t
			case 5:
				return  (1./7.)*(x+r-2.*r*x*x)*(3.*r+7.*x-10.*r*x*x)*P_L(k*r)*P_L(k*sqrt(d))/pow2(d);
        break;
      	//A12t
			case 6:
				return (r/14.)*(x*x-1.)*(3.*r+7.*x-10.*r*x*x)*P_L(k*r)*P_L(k*sqrt(d))/pow2(d);
        break;
      	//A22t
			case 7:
				return (1./14.)*(28.*x*x+r*x*(25.-81.*x*x)+r*r*(1.-27.*x*x+54.*x*x*x*x))*P_L(k*r)*P_L(k*sqrt(d))/pow2(d);
        break;
      	//A23t
			case 8:
				return (r/14.)*(-x*x+1.)*(r-7.*x+6.*r*x*x)*P_L(k*r)*P_L(k*sqrt(d))/pow2(d);
        break;
		//A33t
			case 9:
				return (1./14.)*(-2.*x-r+3.*r*x*x)*(r-7.*x+6.*r*x*x)*P_L(k*r)*P_L(k*sqrt(d))/pow2(d);
        break;
		//a11
			case 10:
				return (-7.*x*x+pow3(r)*x*(-3.+10.*x*x)+3.*r*(x+6.*pow3(x))+r*r*(6.-19.*x*x-8.*pow4(x)))*P_L(k)*P_L(k*r)/(7.*d);
        break;
		//a12
			case 11:
				return r*(-1.+x*x)*(6.*r-7.*(1.+r*r)*x + 8.*r*x*x)*P_L(k)*P_L(k*r)/(14.*d);
        break;
			//a22
			case 12:
				return (-28.*x*x+r*r*r*x*(-13.+41.*x*x)+r*x*(11.+73.*x*x)-2.*r*r*(-9.+31.*x*x+20.*x*x*x*x))*P_L(k)*P_L(k*r)/(14.*d);
        break;
		//a33
			case 13:
				return (7.*x + r*(-6.+7.*r*x-8.*x*x))*(-2.*x+r*(-1.+3.*x*x))*P_L(k)*P_L(k*r)/(14.*d);
        break;

		//B terms
				//B111
			case 14:
				return (r*r/2.)*(x*x-1.)*P_L(k*r)*P_L(k*sqrt(d))/d;
        break;
	//B112
			case 15:
				return (3.*r*r/8.)*pow2(x*x-1.)*P_L(k*r)*P_L(k*sqrt(d))/d;
        break;
		//B121
			case 16:
				return (3.*r*r*r*r/8.)*pow2(x*x-1.)*P_L(k*r)*P_L(k*sqrt(d))/pow2(d);
        break;
			//B122
			case 17:
				return (5.*r*r*r*r/16.)*pow3(x*x-1.)*P_L(k*r)*P_L(k*sqrt(d))/pow2(d);
        break;
			//B211
			case 18:
				return (r/2.)*(r+2.*x-3.*r*x*x)*P_L(k*r)*P_L(k*sqrt(d))/d;
        break;
		//B212
			case 19:
				return (-3.*r/4.)*(x*x-1.)*(-r-2.*x+5.*r*x*x)*P_L(k*r)*P_L(k*sqrt(d))/d;
        break;
			//B221
			case 20:
				return (3.*r*r/4.)*(x*x-1.)*(-2.+r*r+6.*r*x-5.*r*r*x*x)*P_L(k*r)*P_L(k*sqrt(d))/pow2(d);
        break;
			//B222
			case 21:
				return (-3.*r*r/16.)*pow2(x*x-1.)*(6.-30.*r*x-5.*r*r+35.*r*r*x*x)*P_L(k*r)*P_L(k*sqrt(d))/pow2(d);
        break;
		//B312
			case 22:
				return (r/8.)*(4.*x*(3-5.*x*x)+r*(3.-30.*x*x+35.*x*x*x*x))*P_L(k*r)*P_L(k*sqrt(d))/d;
        break;
		//B321
			case 23:
				return (r/8.)*(-8.*x+r*(-12.+36.*x*x+12.*r*x*(3.-5.*x*x)+r*r*(3.-30.*x*x+35.*x*x*x*x)))*P_L(k*r)*P_L(k*sqrt(d))/pow2(d);
        break;
		//B322
			case 24:
				return (3.*r/16.)*(x*x-1.)*(-8.*x+r*(-12.+60.*x*x+20.*r*x*(3.-7.*x*x)+5.*r*r*(1.-14.*x*x+21.*x*x*x*x)))*P_L(k*r)*P_L(k*sqrt(d))/pow2(d);
        break;
		//B422 // HAVE INCLUDED THE 'TYPO' IN 3RD TERM : 6rr -> 6rrx
			case 25:
				return (r/16.)*(8.*x*(-3.+5.*x*x)-6.*r*(3.-30.*x*x+35.*x*x*x*x)+6.*r*r*x*(15.-70.*x*x+63.*x*x*x*x)+r*r*r*(5.-21.*x*x*(5.-15.*x*x+11.*x*x*x*x)))*P_L(k*r)*P_L(k*sqrt(d))/pow2(d);
		// C terms
    break;
			//C11
			case 26:
				return 2.*((r*r*r)*(x*x-1.)*(Dd_spt*F_spt*r*(r*x-1.)-D_spt*Fd_spt*x*d)*P_L(k)*P_L(k*sqrt(d)))/ (Dd_spt*pow2(d));
        break;
			//C12
			case 27:
				return -D_spt*Fd_spt*((r*r*r*r)*pow2(x*x-1.)*P_L(k)*P_L(k*sqrt(d)))/(Dd_spt*pow2(d));
        break;
			//C22
			case 28:
				return ((r*r*r)*(x*x-1.)*(2.*Dd_spt*F_spt*r*(r*x-1.)-D_spt*Fd_spt*(r+4.*x+2.*r*r*x-7.*r*x*x))*P_L(k)*P_L(k*sqrt(d)))/ (Dd_spt*pow2(d));
        break;
			//C23
			case 29:
				return -D_spt*Fd_spt*((r*r*r*r)*pow2(x*x-1.)*P_L(k)*P_L(k*sqrt(d)))/ (Dd_spt*pow2(d));
        break;
				//C33
			case 30:
				return D_spt*Fd_spt*((r*r*r)*(x*x-1.)*(-2.*x+r*(-1.+3.*x*x))*P_L(k)*P_L(k*sqrt(d)))/ (Dd_spt*pow2(d));
        break;
				//C11t
			case 31:
				return 2.*F_spt*(r*(x*x-1.)*(-x+r*(-1.+2.*x*x))*P_L(k*r)*P_L(k*sqrt(d)))/ pow2(d);
        break;
				//C12t
			case 32:
				return -F_spt*(r*r*pow2(x*x-1.)*P_L(k*r)*P_L(k*sqrt(d)))/ pow2(d);
        break;
				//C22t
			case 33:
				return r*(x*x-1.)*(2.*D_spt*Fd_spt*(-x+r*(-1.+2.*x*x)) + Dd_spt*F_spt*(-2.*x+r*(-1.+3.*x*x)))*P_L(k*r)*P_L(k*sqrt(d))/ (Dd_spt*pow2(d));
        break;
				//C23t
			case 34:
				return -D_spt*Fd_spt*(r*r*pow2(x*x-1.)*P_L(k*r)*P_L(k*sqrt(d)))/ (Dd_spt*pow2(d));
        break;
				//C33t
			case 35:
				return D_spt*Fd_spt*(r*(x*x-1.)*(-2.*x+r*(-1.+3.*x*x))*P_L(k*r)*P_L(k*sqrt(d)))/ (Dd_spt*pow2(d));
        break;
				//c11
			case 36:
				return  -2.*r*(1.-x*x)*(D_spt*Fd_spt*r*(r*x-1.)-Dd_spt*F_spt*x*d)*P_L(k*r)*P_L(k)/(Dd_spt*d);
        break;
				//c12
			case 37:
				return -D_spt*Fd_spt*r*r*pow2(-1.+x*x)*P_L(k*r)*P_L(k)/(Dd_spt*d);
        break;
				//c22
			case 38:
				return r*(-1.+x*x)*(-2.*Dd_spt*F_spt*x*d+D_spt*Fd_spt*(-2.*x+2.*r*r*x+3.*r*(-1.+x*x)))*P_L(k*r)*P_L(k)/(Dd_spt*d);
        break;
				//c33
			case 39:
				return D_spt*Fd_spt*r*(-1.+x*x)*(-2.*x+r*(-1.+3.*x*x))*P_L(k*r)*P_L(k)/(Dd_spt*d);
        break;
			default:
			warning("SPT: invalid indices, a = %d\n", a);
				return 0;
}
}




/* function to select multipole or u-dependent Redshift PS */
// it gives either u^(2n) or see factL in SpecialFunctions.cpp


// e = 0:  RSD PS
// e = 1 : Monopole
// e = 2 : Quadrupole
// e = 3 : Hexdecapole

static real u0(real k, real u, real F0, int n, int a){
	if (a >= 1) {
		return factL(k,u,F0,1.,n,a,6); }
		else
			return pow(u,2*n);
}

static real midintabc(double u0[],real bl, int a, const PowerSpectrum& P_L, real k, real u, real r) {
  real KMAX = QMAXp/k;
  real KMIN = QMINp/k;
	real YMIN = Max(-1., (1.+r*r-KMAX*KMAX)/2./r);
  real YMAX = ATS(k,r);
	real F0,D1,X1;
			F0 = rempf*fdgp_spt;
			D1 = D_spt;
			X1 = 1.;
      // real u01 = u0(k,u,F0,1,a);
      // real u02 = u0(k,u,F0,2,a);
      // real u03 = u0(k,u,F0,3,a);
      // real u04 = u0(k,u,F0,4,a);
      real u01 = u0[0];//u0(k,u,1.,1,a);
      real u02 = u0[1];//u0(k,u,1.,2,a);
      real u03 = u0[2];//u0(k,u,1.,3,a);
      real u04 = u0[3];//u0(k,u,1.,4,a);

	//A term
		return     pow4(D1/dnorm_spt)*2.*(F0*u01*pow2(bl)*(Integrate(bind(ABC, 1, cref(P_L), k, r, _1),  YMIN, YMAX ,  1e-3*P_L(k))
										   + Integrate(bind(ABC, 5, cref(P_L), k,  r, _1),  YMIN, YMAX ,  1e-3*P_L(k))
										   + Integrate(bind(ABC, 10, cref(P_L), k, r, _1),  YMIN, YMAX ,  1e-3*P_L(k)))

							+ F0*F0*u01*bl*(Integrate(bind(ABC, 2, cref(P_L), k,  r, _1),  YMIN, YMAX , 1e-3*P_L(k))
										   + Integrate(bind(ABC, 6, cref(P_L), k,  r, _1), YMIN, YMAX ,  1e-3*P_L(k))
										   + Integrate(bind(ABC, 11, cref(P_L), k,  r, _1), YMIN, YMAX ,  1e-3*P_L(k)))

							+ F0*F0*u02*bl*(Integrate(bind(ABC, 3, cref(P_L), k,  r, _1),   YMIN, YMAX , 1e-3*P_L(k))
										   + Integrate(bind(ABC, 7, cref(P_L), k,  r, _1),   YMIN, YMAX ,  1e-3*P_L(k))
										   + Integrate(bind(ABC, 12, cref(P_L), k,  r, _1),  YMIN, YMAX ,  1e-3*P_L(k)))

					     + F0*F0*F0*u02*(Integrate(bind(ABC, 2, cref(P_L), k,  r, _1),  YMIN, YMAX ,  1e-3*P_L(k))
										   + Integrate(bind(ABC, 8, cref(P_L), k,  r, _1),   YMIN, YMAX ,  1e-3*P_L(k))
										   + Integrate(bind(ABC, 11, cref(P_L), k, r, _1), YMIN, YMAX ,  1e-3*P_L(k)))

						   + F0*F0*F0*u03*(Integrate(bind(ABC, 4, cref(P_L), k, r, _1),   YMIN, YMAX , 1e-3*P_L(k))
										   + Integrate(bind(ABC, 9, cref(P_L), k,  r, _1),  YMIN, YMAX ,  1e-3*P_L(k))
										   + Integrate(bind(ABC, 13, cref(P_L), k,  r, _1),YMIN, YMAX ,  1e-3*P_L(k))))

	//B term
		  +      pow4(D1/dnorm_spt)*2.*(F0*F0*pow2(bl)*u01*Integrate(bind(ABC, 14,  cref(P_L), k, r, _1), YMIN, YMAX , 1e-3* P_L(k))

						  -F0*F0*F0*bl*u01*(Integrate(bind(ABC, 15, cref(P_L), k, r, _1), YMIN, YMAX , 1e-3*P_L(k))
											    +Integrate(bind(ABC, 16, cref(P_L),  k,  r, _1),YMIN, YMAX , 1e-3*P_L(k)))

						+F0*F0*F0*F0*u01*Integrate(bind(ABC, 17, cref(P_L),  k, r, _1), YMIN, YMAX , 1e-3*P_L(k))

							  +F0*F0*pow2(bl)*u02*Integrate(bind(ABC, 18,  cref(P_L), k,  r, _1),YMIN, YMAX ,  1e-3*P_L(k))

						 -F0*F0*F0*bl*u02*(Integrate(bind(ABC, 19, cref(P_L), k,  r, _1), YMIN, YMAX , 1e-3*P_L(k))
										      + Integrate(bind(ABC, 20, cref(P_L),  k,  r, _1), YMIN, YMAX ,  1e-3*P_L(k)))

					  +F0*F0*F0*F0*u02*Integrate(bind(ABC, 21, cref(P_L),  k,  r, _1),YMIN, YMAX , 1e-3*P_L(k))

					    - F0*F0*F0*bl*u03*(Integrate(bind(ABC, 22,  cref(P_L), k,  r, _1), YMIN, YMAX ,  1e-3*P_L(k))
											  + Integrate(bind(ABC, 23,  cref(P_L), k,  r, _1), YMIN, YMAX , 1e-3*P_L(k)))

					  + F0*F0*F0*F0*u03*Integrate(bind(ABC, 24,  cref(P_L), k,  r, _1),YMIN, YMAX ,  1e-3*P_L(k))

							+F0*F0*F0*F0*u04*Integrate(bind(ABC, 25,  cref(P_L), k,  r, _1), YMIN, YMAX , 1e-3*P_L(k)))

	//C - nDGP analytic

		  +       X1*D1*D1*2./pow4(dnorm_spt)*(F0*pow2(bl)*u01*(Integrate(bind(ABC, 26, cref(P_L), k, r, _1),   YMIN, YMAX ,  1e-3*P_L(k))
										+ Integrate(bind(ABC, 31, cref(P_L), k, r, _1),   YMIN, YMAX ,  1e-3*P_L(k))
										+ Integrate(bind(ABC, 36, cref(P_L), k, r, _1), YMIN, YMAX ,  1e-3*P_L(k))) +

							F0*F0*bl*u01*(Integrate(bind(ABC, 27, cref(P_L), k, r, _1),   YMIN, YMAX ,  1e-3*P_L(k))
										+ Integrate(bind(ABC, 32, cref(P_L), k, r, _1),   YMIN, YMAX ,  1e-3*P_L(k))
										+ Integrate(bind(ABC, 37, cref(P_L), k, r, _1),  YMIN, YMAX ,  1e-3*P_L(k))) +

						   F0*F0*bl*u02*(Integrate(bind(ABC, 28, cref(P_L), k, r, _1),   YMIN, YMAX ,  1e-3*P_L(k))
										+ Integrate(bind(ABC, 33, cref(P_L), k, r, _1),   YMIN, YMAX ,  1e-3*P_L(k))
										+ Integrate(bind(ABC, 38, cref(P_L), k,  r, _1), YMIN, YMAX ,  1e-3*P_L(k))) +

						F0*F0*F0*u02*(Integrate(bind(ABC, 29, cref(P_L), k, r, _1),   YMIN, YMAX ,  1e-3*P_L(k))
										+ Integrate(bind(ABC, 34, cref(P_L), k, r, _1),   YMIN, YMAX ,  1e-3*P_L(k))
										+ Integrate(bind(ABC, 37, cref(P_L), k,  r, _1),  YMIN, YMAX ,  1e-3*P_L(k))) +

					   F0*F0*F0*u03*(Integrate(bind(ABC, 30, cref(P_L), k, r, _1),   YMIN, YMAX ,  1e-3*P_L(k))
									   + Integrate(bind(ABC, 35, cref(P_L), k, r, _1),   YMIN, YMAX ,  1e-3*P_L(k))
									   + Integrate(bind(ABC, 39, cref(P_L), k,  r, _1),  YMIN, YMAX ,  1e-3*P_L(k))));

}



/* Non vanishing A, B and C terms */
//numerical calculations
// Notes for numerical AB terms : DEFINITION OF THETA=  DEL V / aH (and is accounted for in definition of Cross Bispectrum)
// Integrating them over angle
static real ABCn(double u0[], int y, real bl, int a, const PowerSpectrum& P_L, real k, real u, real r) {
  double myresult = 0.;
	double temp_abc;
  double abc[10];
  // Set the limits of angular integration
  double KMAX = QMAXp/k;
  double KMIN = QMINp/k;
  double YMIN = Max(XMIN, (1.+r*r-KMAX*KMAX)/2./r);
  double YMAX = ATSn(k,r);
  // Multipole factors
  real u01 = u0[0];//u0(k,u,1.,1,a);
  real u02 = u0[1];//u0(k,u,1.,2,a);
  real u03 = u0[2];//u0(k,u,1.,3,a);
  real u04 = u0[3];//u0(k,u,1.,4,a);
  // Find those corresponding limits in the GL Quad abscissae (See SpecialFunctions.h)
  int x_min = searchnearest(x128, YMIN);
  int x_max = searchnearest(x128, YMAX);
  for( int i = x_min; i <= x_max; i++ )
    				{
    		 		double x = x128[i];
    		 		double d = 1.+ r*r - 2.*r*x;
    		 		if(d < 1e-5){
    			 	temp_abc=0.;
    		 		}
    		  	else {
              // 1st order kernels F1/G1(k-p)
              abc[0] = F1kmp_nk[i*n2 + y];
              abc[1] = rempf*G1kmp_nk[i*n2 + y]/bl;
              // 1st order kernels F1/G1(p)
              abc[2] = F1p_nk[i*n2 + y];
              abc[3] = rempf*G1p_nk[i*n2 + y]/bl;
              //symmetrized 2nd order kernels for ps F2/G2(p,k-p)
              abc[4] = F2_nk[i*n2 + y];
              abc[5] = rempf*G2_nk[i*n2 + y]/bl;
              //symmetrized 2nd order kernels F2/G2(-p,k)
              abc[6] = F2B_nk[i*n2 + y];
              abc[7] = rempf*G2B_nk[i*n2 + y]/bl;
              //symmetrized 2nd order kernels F2/G2(-k,k-p)
              abc[8] = F2C_nk[i*n2 + y];
              abc[9] = rempf*G2C_nk[i*n2 + y]/bl;


// A terms
//u^2
    		  temp_abc =  pow3(bl)/d*(-u01*(F1_nk*r*(-2.*abc[8]*abc[1]*r*(r*x-1.)+abc[9]*(2.*abc[0]*x*d+abc[1]*r*(1.-x*x)))*P_L(k)*P_L(k*sqrt(d))
            						+ abc[4]*r*(-2.*abc[2]*abc[1]*r*(-1.+r*x)+abc[3]*(2.*abc[0]*x*d+abc[1]*r*(1.-x*x)))*P_L(k*r)*P_L(k*sqrt(d))
            						+ F1_nk*r*(-2.*abc[7]*abc[2]*r*(r*x-1.)+abc[3]*(2.*abc[6]*x*d+abc[7]*r*(1.-x*x)))*P_L(k)*P_L(k*r))
//u^4
            -u02*(r*(2.*abc[8]*rempf*G1_nk/bl*abc[1]*r*(r*x-1.)+abc[9]*(rempf*G1_nk/bl*(-2.*abc[0]*x*d+abc[1]*r*(x*x-1.))+F1_nk*abc[1]*(-2.*x+r*(-1.+3*x*x))))*P_L(k)*P_L(k*sqrt(d))
                           +r*(abc[4]*abc[1]*abc[3]*(-2.*x+r*(-1.+3.*x*x))+abc[5]*(2.*abc[2]*abc[1]*r*(r*x-1.)+abc[3]*(-2.*abc[0]*x*d+abc[1]*r*(x*x-1))))*P_L(k*r)*P_L(k*sqrt(d))
                           +r*(abc[3]*abc[7]*F1_nk*(-2.*x-r+3.*x*x*r)+2.*abc[2]*rempf*G1_nk/bl*abc[7]*r*(r*x-1.)+rempf*G1_nk/bl*abc[3]*(-2.*abc[6]*x*d+abc[7]*r*(x*x-1)))*P_L(k)*P_L(k*r))
//u^6
           -u03*(rempf*G1_nk/bl*abc[1]*abc[9]*r*(r+2.*x-3.*r*x*x)*P_L(k)*P_L(k*sqrt(d))
                        +abc[3]*abc[1]*abc[5]*r*(r+2.*x-3.*r*x*x)*P_L(k*r)*P_L(k*sqrt(d))
                        +abc[3]*rempf*G1_nk/bl*abc[7]*r*(r+2.*x-3.*r*x*x)*P_L(k*r)*P_L(k)))


// B terms
          +pow4(bl)*P_L(k*r)*P_L(k*sqrt(d))/(16.*d*d)*
//u^2
          (u01*(abc[1]*abc[3]*r*r*(x*x-1.)*(2.*abc[0]*d*(4.*abc[2]+3.*abc[3]*(x*x-1.))
                         +abc[1]*r*r*(x*x-1.)*(6.*abc[2]+5.*abc[3]*(x*x-1))))
//u^4
          +u02*(abc[1]*abc[3]*r*(-4.*abc[0]*d*(3.*abc[3]*(x*x-1.)*(-2.*x+r*(5.*x*x-1.))+abc[2]*(-4*x+r*(-2.+6.*x*x)))
                					   -3.*abc[1]*r*(x*x-1.)*(4.*abc[2]*(2.-6.*r*x+r*r*(-1.+5.*x*x))
                						 +abc[3]*(x*x-1)*(6.-30.*r*x+5.*r*r*(-1.+7.*x*x)))))
//u^6
          +u03*(abc[1]*abc[3]*r*(3.*abc[3]*abc[1]*(x*x-1.)*(-8.*x+12.*r*(5.*x*x-1.)-20.*r*r*x*(7.*x*x-3.)+5.*pow3(r)*(1.-14.*x*x+21.*pow4(x)))
                              +2.*abc[0]*abc[3]*(4.*x*(3.-5.*x*x)+pow3(r)*(3.-30.*x*x+35.*pow4(x))+r*(3.-54.*x*x+75.*pow4(x))+r*r*(6.*x+40.*pow3(x)-70.*pow5(x)))
                              +2.*abc[2]*abc[1]*(-8.*x+12.*r*(3.*x*x-1.)+r*r*(36.*x-60.*pow3(x))+pow3(r)*(3.-30.*x*x+35.*pow4(x)))))
//u^8
          +u04*(pow2(abc[1]*abc[3])*r*(8.*x*(5.*x*x-3.)-6.*r*(3.-30.*x*x+35.*pow4(x))
                     				+6.*r*r*x*(15.-70.*x*x+63.*pow4(x)) 		 // HAVE INCLUDED THE 'TYPO' IN TERM : 6rr -> 6rrx
                     				+pow3(r)*(5.-105.*x*x+315.*pow4(x)-231.*pow6(x)))));
    				}
    				myresult += w128[i] * temp_abc;
    				}
          return  2.*myresult;
}


/*Generalised A and B correction term */

/*Analytical A and B*/
// a = 1 : LCDM
// a = 2 : MG

// e = 0:  RSD PS
// e = 1 : Monopole
// e = 2 : Quadrupole
// e = 3 : Hexdecapole

//U = u for 2d PS
//U = sigma_v for multipoles

// For A or B or C seperated : just remove as necessary.
real SPT::AB(real k, real bl, real U, int a) const{
      int n3 = 300;
      real KMAX = QMAXp/k;
      real KMIN = QMINp/k;
      double y[n3];
      double integrand[n3];
      double u0x[4];
      u0x[0]=u0(k,U,1.,1,a);
      u0x[1]=u0(k,U,1.,2,a);
      u0x[2]=u0(k,U,1.,3,a);
      u0x[3]=u0(k,U,1.,4,a);
      for (int i = 0; i<n3; i++){
      y[i] = KMIN * exp(i*log(KMAX/KMIN)/(n3-1.));
      integrand[i] = midintabc(u0x,bl, a ,cref(P_L), k, U, y[i]);
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

//U = u for 2d PS
//U = sigma_v for multipoles

real SPT::ABn(double kmin, double kmax, real bl, real k, real U, int a) const{
	real KMAX = QMAXp/k;
  real KMIN = QMINp/k;
  int y1 = num_selec_mag(kmin, kmax, KMIN);
  int y2 = num_selec_mag(kmin, kmax, KMAX);
  double y[n2];
  double integrand[n2];
  double u0x[4];
  u0x[0]=u0(k,U,1.,1,a);
  u0x[1]=u0(k,U,1.,2,a);
  u0x[2]=u0(k,U,1.,3,a);
  u0x[3]=u0(k,U,1.,4,a);
  for (int i = y1; i<=y2; i++){
  y[i] = QMINp/kmax * exp(i*log(QMAXp*kmax/(QMINp*kmin))/(n2*1.-1.));
  integrand[i] = ABCn(u0x,i, bl, a, cref(P_L), k, U, y[i]);
  }
double res = 0.;
  for( int i = y1+1; i <= y2; ++ i ){
res += 0.5 * (y[i] - y[i-1])*(integrand[i] + integrand[i-1]);
}
  return  k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt) * res;
}


// Used for scale dep interpolation scheme

static void ABns_selec(double abarray[], int y, const PowerSpectrum& P_L, real k, real r) {
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
  double YMAX = ATSn(k,r);
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

              double dampexpa = 1.;
              double dampexpb = 1.;
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


void SPT::ABs_selec(double myarray[], real k) const {
  	real KMAX = QMAXp/k;
    real KMIN = QMINp/k;
    int y1 = num_selec_mag(k, k, KMIN);
    int y2 = num_selec_mag(k, k, KMAX);
    double y[n2];
    double integrand[n2][14];
    double abarray[14];
    for (int i = y1; i<=y2; i++){
    y[i] = KMIN * exp(i*log(KMAX/KMIN)/(n2*1.-1.));
    ABns_selec(abarray,i, cref(P_L), k, y[i]);
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
