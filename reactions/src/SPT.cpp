#if HAVE_CONFIG_H
# include <config.h>
#endif


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
#include <random>
#include <functional>

using std::cref;
using std::bind;


////////// PERTURBATION THEORY CONSTRUCTED STATISTICS ////////////

/*Load the linear power spectrum */

SPT::SPT(const Cosmology& C_, const PowerSpectrum& P_L_, real epsrel_)
: C(C_), P_L(P_L_)
{
    epsrel = epsrel_;
}

// Function to vary growth rate f as a free parameter instead of predicted from theory.
double rempf;

//LCDM
void SPT::remp(real f) const{
rempf = f/fl_spt;
}

//nDGP
void SPT::rempdgp(real f) const{
rempf = f/fdgp_spt;
}

// numerical  --- TODO: not in RSD spectra yet --- should be easy to implement
void SPT::rempn(real f) const{
rempf = f/(-G1_nk/F1_nk);
}



///////////////////////////////ANALYTIC PREDICTIONS IN REAL SPACE  (EDS APPROXIMATED) /////////////////////////////////

/* P_{ab}(k) 1-loop: analytical LCDM and nDGP */
/* NOTE THAT THE DEFINITIONS OF VEL DIV ARE THETA= DEL V / (aH) */

/* The division by dnorm_spt is my own normalisation of P_L - it requires z=0 always in cosmology class */
// 1: P_dd, 2: P_dt , 3: P_tt ; LCDM
// 4:  P_dd, 5: P_dt , 6: P_tt ; nDGP
real SPT::PLOOP(real k, int a) const {
    switch(a) {
        //LCDM
        case 1:
            return  pow2(Dl_spt/dnorm_spt)*P_L(k) + P13_dd(k)+P22_dd(k);
            break;
        case 2:
            return  -rempf*fl_spt*pow2(Dl_spt/dnorm_spt)*P_L(k) + P13_dt(k) + P22_dt(k);
            break;
        case 3:
            return  pow2(rempf*fl_spt*Dl_spt/dnorm_spt)*P_L(k) + P13_tt(k) + P22_tt(k);
            break;
        //DGP
        case 4:
            return   pow2(D_spt/dnorm_spt)*P_L(k) + P13D_dd(k)+P22D_dd(k);
            break;
        case 5:
            return  -fdgp_spt*pow2(D_spt/dnorm_spt)*P_L(k) + P13D_dt(k)+ P22D_dt(k);
            break;
        case 6:
            return   pow2(fdgp_spt*D_spt/dnorm_spt)*P_L(k) + P13D_tt(k)+P22D_tt(k);
            break;
        default:
            warning("SPT: invalid indices, a = %d\n", a);
            return 0;
    }
}


/////////////////////////////// PREDICTIONS IN REDSHIFT SPACE - Analytic /////////////////////////////////

/* P(k,u) RSD KAISER : */
// a =1 : LCDM
// a =2 : nDGP
// u = k.z/|k|, z being a unit vector along the line of sight axis
// bl is linear  bias
real SPT::PRSD(real k, real u, real bl, int a) const {
	double F0,D0;
	switch (a) {
		case 1:
			F0=rempf*fl_spt;
			D0= Dl_spt;
			break;
		case 2:
			F0=rempf*fdgp_spt;
			D0=D_spt;
			break;
    default:
      warning("SPT: invalid indices, a = %d\n", a);
      return 0;
				  }
	return    pow2(D0*bl/dnorm_spt)*(P_L(k) + 2.*(F0/bl)*u*u*P_L(k) + pow2(F0/bl)*pow4(u)* P_L(k));
	 }


/* P(k,u) ANALYTIC RSD TNS */
// a = 1 : LCDM
// a = 2 : nDGP
// u = k.z/|k|, z being a unit vector along the line of sight axis
//bl is linear  bias
real SPT::PTNS(real k, real u, real bl, real sigma_v, int a) const {
    switch(a) {
		case 1:
			return  DFOG(k,u,sigma_v,a)*(pow2(bl)*PLOOP(k,1) - 2.*u*u*bl*PLOOP(k,2) + pow4(u)* PLOOP(k,3) + AB_mu(k, bl, u, a));
      break;
    case 2:
      return	DFOG(k,u,sigma_v,a)*(pow2(bl)*PLOOP(k,4) - 2.*u*u*bl*PLOOP(k,5) + pow4(u)*PLOOP(k,6) + AB_mu(k, bl, u, a));
      break;
        default:
            warning("SPT: invalid indices, a = %d \n", a );
            return 0;
    }
}


/*Multipoles*/
/*  LCDM/nDGP  Kaiser Multipoles with DFoG term (set sigma_v=0 for linear kaiser multipoles) */
//bl is linear  bias
// a = 1 : Monopole
// a = 2 : Quadrupole
// a = 3 : Hexdecapole
real SPT::KASM(real k, real bl, real sigma_v, int a) const {
	return   pow2(D_spt*bl/dnorm_spt)*(factL(k,sigma_v,1.,1.,0,a,7) + 2.*rempf*fdgp_spt/bl*factL(k,sigma_v,1.,1.,1,a,7) + pow2(rempf*fdgp_spt/bl)*factL(k,sigma_v,1.,1.,2,a,7))*P_L(k);
}


/*  LCDM  TNS  Multipoles */
//bl is linear bias
// a = 1 : Monopole
// a = 2 : Quadrupole
// a = 3 : Hexdecapole
real SPT::PTNSM(real k, real bl, real sigma_v, int a) const {
	return  pow2(bl)*factL(k,sigma_v,1.,1.,0,a,7) *PLOOP(k,1) - 2.*factL(k,sigma_v,1.,1.,1,a,7)*bl*PLOOP(k,2) + factL(k,sigma_v,1.,1.,2,a,7) * PLOOP(k,3) + AB(k, bl, sigma_v, a);
    }


/* biased tracers */
//Q-bias
real SPT::PTNSMq(real k, double barr[], real sigma_v, int a) const {
    double bk = barr[0]*sqrt((1+barr[1]*pow2(k))/(1+barr[2]*k));
  return  pow2(bk)*factL(k,sigma_v,1.,1.,0,a,7) *PLOOP(k,1) - 2.*factL(k,sigma_v,1.,1.,1,a,7) *bk*PLOOP(k,2) + factL(k,sigma_v,1.,1.,2,a,7) * PLOOP(k,3) + AB(k, bk, sigma_v, a);
      }

/* Lagrangian bias of Roy and MacDonald See Eq.23 1607.03150 for example*/
real SPT::PTNSMl(real k, double barr[], real sigma_v, int a) const {
  return  factL(k,sigma_v,1.,1.,0,a,7)*(pow2(barr[0])*PLOOP(k,1)+Lag_bias(1,k,barr) + barr[2]) - 2.*factL(k,sigma_v,1.,1.,1,a,7)*(barr[0]*PLOOP(k,2)+Lag_bias(2,k,barr)) + factL(k,sigma_v,1.,1.,2,a,7)* PLOOP(k,3) +  AB(k, barr[0], sigma_v, a);
      }

/*  analytic nDGP TNS  Multipoles  */
// bl is linear bias
// a = 1 : Monopole
// a = 2 : Quadrupole
// a = 3 : Hexdecapole

real SPT::PTNSMnDGP(real k, real bl, real sigma_v, int a) const {
	return  pow2(bl)*factL(k,sigma_v,1.,1.,0,a,7) *PLOOP(k,4) - 2.*factL(k,sigma_v,1.,1.,1,a,7) *bl*PLOOP(k,5) + factL(k,sigma_v,1.,1.,2,a,7) * PLOOP(k,6) + AB(k, bl, sigma_v, a);
}

// Q-bias for nDGP
real SPT::PTNSMnDGPq(real k, double barr[], real sigma_v, int a) const {
  double bk = barr[0]*sqrt((1+barr[1]*pow2(k))/(1+barr[2]*k));
	return  pow2(bk)*factL(k,sigma_v,1.,1.,0,a,7) *PLOOP(k,4) - 2.*factL(k,sigma_v,1.,1.,1,a,7) *bk*PLOOP(k,5) + factL(k,sigma_v,1.,1.,2,a,7) * PLOOP(k,6) +  AB(k, bk, sigma_v, a);
}




/*All functions below are the components used to construct the above spectra :
 In order of appearance :

 1 LOOP CORRECTIONS : LCDM , nDGP, Numerical
 TNS A,B and nDGP correction terms
 modified gravity TNS and Kaiser with linear, Q-bias and Lag bias (GR modelling only for Lag Bias)*/


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
return Integrate(bind(f22_dd, cref(P_L), k, r, std::placeholders::_1), YMIN, YMAX ,  1e-4*P_L(k));
}

/* P_{\delta\delta}^{(22)} */
real SPT::P22_dd(real k) const {
 	real KMAX = QMAXp/k;
	real KMIN = QMINp/k;
    return k*k*k/(4*M_PI*M_PI) * pow4(Dl_spt/dnorm_spt) * Integrate<ExpSub>(bind(midintdd, cref(P_L), k, std::placeholders::_1), KMIN, KMAX,epsrel);
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
    return k*k*k/(252.*4*M_PI*M_PI) *P_L(k) *  pow4(Dl_spt/dnorm_spt)* Integrate<ExpSub>(bind(f13_dd, cref(P_L), k, std::placeholders::_1), KMIN, KMAX,epsrel);
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
return Integrate(bind(f22_dt, cref(P_L), k, r, std::placeholders::_1), YMIN, YMAX ,  1e-4*P_L(k));
}


/* P_{\delta\theta}^{(22)} */
real SPT::P22_dt(real k) const {
	real KMAX = QMAXp/k;
	real KMIN = QMINp/k;
    return  k*k*k/(98.*4.*M_PI*M_PI) *pow4(Dl_spt/dnorm_spt)*(-rempf*fl_spt)* Integrate<ExpSub>(bind(midintdt, cref(P_L), k, std::placeholders::_1), KMIN, KMAX, epsrel);
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
    return k*k*k/(252.*4.*M_PI*M_PI) *P_L(k) *(-rempf*fl_spt)*pow4(Dl_spt/dnorm_spt)* Integrate<ExpSub>(bind(f13_dt, cref(P_L), k, std::placeholders::_1), KMIN, KMAX, epsrel);
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
return Integrate(bind(f22_tt, cref(P_L), k, r, std::placeholders::_1), YMIN, YMAX ,  1e-4*P_L(k));
}

 /* P_{\theta\theta}^{(22)} */
 real SPT::P22_tt(real k) const {
	real KMAX = QMAXp/k;
	real KMIN = QMINp/k;
	 return k*k*k/(98*4*M_PI*M_PI) *  pow4(Dl_spt/dnorm_spt)*pow2(rempf*fl_spt)*Integrate<ExpSub>(bind(midintt, cref(P_L), k, std::placeholders::_1), KMIN, KMAX, epsrel);
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
    return k*k*k/(84*4*M_PI*M_PI) *P_L(k) *pow4(Dl_spt/dnorm_spt)*pow2(rempf*fl_spt) * Integrate<ExpSub>(bind(f13_tt, cref(P_L), k, std::placeholders::_1), KMIN, KMAX, epsrel);
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
return Integrate(bind(Df22_dd, cref(P_L), k, r, std::placeholders::_1), YMIN, YMAX ,  1e-4*P_L(k));
}


/* P_{\delta\delta}^{(22)} */
real SPT::P22D_dd(real k) const {
 	real KMAX = QMAXp/k;
	real KMIN = QMINp/k;
    return k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt)* Integrate<ExpSub>(bind(Dmidintdd, cref(P_L), k, std::placeholders::_1), KMIN, KMAX, epsrel);
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
    return  k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt) *P_L(k) * Integrate<ExpSub>(bind(Df13_dd, cref(P_L), k, std::placeholders::_1), KMIN, KMAX, epsrel);
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
return Integrate(bind(Df22_dt, cref(P_L), k, r, std::placeholders::_1), YMIN, YMAX ,  1e-4*P_L(k));
}


/* P_{\delta\theta}^{(22)} */
real SPT::P22D_dt(real k) const {
        real KMAX = QMAXp/k;
    	real KMIN = QMINp/k;
    return k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt) * Integrate<ExpSub>(bind(Dmidintdt, cref(P_L), k, std::placeholders::_1), KMIN, KMAX, epsrel);
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
    return k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt) *P_L(k) * Integrate<ExpSub>(bind(Df13_dt, cref(P_L), k, std::placeholders::_1), KMIN, KMAX,  epsrel);
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
return Integrate(bind(Df22_tt, cref(P_L), k, r, std::placeholders::_1), YMIN, YMAX ,  1e-4*P_L(k));
}

/* P_{\theta\theta}^{(22)} */
real SPT::P22D_tt(real k) const {
	real KMAX = QMAXp/k;
	real KMIN = QMINp/k;
    return k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt) * Integrate<ExpSub>(bind(Dmidintt, cref(P_L), k, std::placeholders::_1), KMIN, KMAX, epsrel);
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

return P_L(k*r) * s;
}


/* P_{\theta\theta}^{(13)} */

real SPT::P13D_tt(real k) const {
   	real KMAX = QMAXp/k;
	real KMIN = QMINp/k;
    return k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt) *P_L(k) * Integrate<ExpSub>(bind(Df13_tt, cref(P_L), k, std::placeholders::_1), KMIN, KMAX, epsrel);
}




/* Higher order bias terms in Roy and Mac Donald model - GR only -  See Eq.24 and 25 of 1607.03150 for example*/

// P_b2
static real lbias_selec(const PowerSpectrum& P_L, int a,  real bias[], real k, real r, real x) {
    real d = 1. + r*r - 2.*r*x;
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
    default:
    warning("SPT: invalid indices, a = %d \n", a);
        return 0;
            }
}
}

static real midpb2(const PowerSpectrum& P_L, int a, real bias[], real k, real r) {
    real KMAX = QMAXp/k;
    real KMIN = QMINp/k;
    real YMIN = Max(-1., (1.+r*r-KMAX*KMAX)/2./r);
    real YMAX = Min(0.999999999, (1.+r*r-KMIN*KMIN)/2./r);
return Integrate(bind(lbias_selec, cref(P_L), a, bias, k, r, std::placeholders::_1), YMIN, YMAX , 1e-3);
}

//P_bs2
real SPT::Lag_bias(int a, real k, real bias[]) const {
 	real KMAX = QMAXp/k;
	real KMIN = QMINp/k;
    return k*k*k/(4*M_PI*M_PI) * pow4(Dl_spt/dnorm_spt) * Integrate<ExpSub>(bind(midpb2, cref(P_L), a, bias, k, std::placeholders::_1), KMIN, KMAX,epsrel);
}



/////// Numerical 1-loop spectra /////////////

static double ploopn2_mgdd( const PowerSpectrum& P_L, double vars[], int model,  double k, double r, double x){
  double kargs[4],kv[3],xv[3], p22,p13;
  // tolerance for ode solver (see SpecialFunctions.cpp, initn2). This encounters a singularity if k'.-k' = -1 exactly.
  double tol = 1e-8;
  // The integrated |k'| is parametrised as k*r
  IOW iow;
        kv[0] = k;
        kv[1] = k*r;
        kv[2] = kv[1];
        xv[0] = -1. + tol;
        xv[1] = x;
        xv[2] = -x;
        kargs[0] = sqrt(kv[1]*kv[1]+kv[0]*kv[0]-2.*kv[1]*kv[0]*xv[1]);
        kargs[2] = sqrt(kv[1]*kv[1]+2.*kv[1]*kv[2]*xv[0]+kv[2]*kv[2]);
        kargs[1] = sqrt(kv[2]*kv[2]+2.*kv[2]*kv[0]*xv[2]+kv[0]*kv[0]);
        kargs[3] = sqrt(kv[0]*kv[0]+2.*kv[0]*kv[1]*xv[1]+kv[1]*kv[1]);
        iow.initn2(vars[0],kv,xv,kargs,vars[1],vars[2],vars[3],vars[4],vars[6],model);
        p22 = pow2(F2_nk);
        p13 = F1_nk*F3_nk;

    return pow2(r)*2.*P_L(k*r)*(P_L(kargs[0])*p22 + 3.*P_L(k)*p13);
}


static double ploopn2_mgdt( const PowerSpectrum& P_L, double vars[], int model, double k, double r, double x){
  double kargs[4],kv[3],xv[3], p22,p13;
  // tolerance for ode solver (see SpecialFunctions.cpp, initn2). This encounters a singularity if k'.-k' = -1 exactly.
  double tol = 1e-8;
  // The integrated |k'| is parametrised as k*r
  IOW iow;
        kv[0] = k;
        kv[1] = k*r;
        kv[2] = kv[1];
        xv[0] = -1. + tol ;
        xv[1] = x;
        xv[2] = -x;
        kargs[0] = sqrt(kv[1]*kv[1]+kv[0]*kv[0]-2.*kv[1]*kv[0]*xv[1]);
        kargs[2] = sqrt(kv[1]*kv[1]+2.*kv[1]*kv[2]*xv[0]+kv[2]*kv[2]);
        kargs[1] = sqrt(kv[2]*kv[2]+2.*kv[2]*kv[0]*xv[2]+kv[0]*kv[0]);
        kargs[3] = sqrt(kv[0]*kv[0]+2.*kv[0]*kv[1]*xv[1]+kv[1]*kv[1]);
        iow.initn2(vars[0],kv,xv,kargs,vars[1],vars[2],vars[3],vars[4],vars[6],model);
        p22 = G2_nk*F2_nk;
        p13 = 0.5*(G1_nk*F3_nk + F1_nk*G3_nk);
    return pow2(r)*2.*P_L(k*r)*(P_L(kargs[0])*p22 + 3.*P_L(k)*p13);
}

static double ploopn2_mgtt(const PowerSpectrum& P_L, double vars[], int model, double k, double r, double x){
  double kargs[4],kv[3],xv[3], p22,p13;
  // tolerance for ode solver (see SpecialFunctions.cpp, initn2). This encounters a singularity if k'.-k' = -1 exactly.
  double tol = 1e-8;
  // The integrated |k'| is parametrised as k*r
  IOW iow;
        kv[0] = k;
        kv[1] = k*r;
        kv[2] = kv[1];
        xv[0] = -1. + tol;
        xv[1] = x;
        xv[2] = -x;
        kargs[0] = sqrt(kv[1]*kv[1]+kv[0]*kv[0]-2.*kv[1]*kv[0]*xv[1]);
        kargs[2] = sqrt(kv[1]*kv[1]+2.*kv[1]*kv[2]*xv[0]+kv[2]*kv[2]);
        kargs[1] = sqrt(kv[2]*kv[2]+2.*kv[2]*kv[0]*xv[2]+kv[0]*kv[0]);
        kargs[3] = sqrt(kv[0]*kv[0]+2.*kv[0]*kv[1]*xv[1]+kv[1]*kv[1]);
        iow.initn2(vars[0],kv,xv,kargs,vars[1],vars[2],vars[3],vars[4],vars[6],model);
        p22 = pow2(G2_nk);
        p13 = G1_nk * G3_nk;
    return pow2(r)*2.*P_L(k*r)*(P_L(kargs[0])*p22 + 3.*P_L(k)*p13);
}


// pseudo 1-loop matter power spectrum (1812.05594)
static double ploopn2_mgdd_pseudo( const PowerSpectrum& P_L, double vars[], int model, double k, double r, double x){
  double kargs[4],kv[3],xv[3], p22,p13;
  // tolerance for ode solver (see SpecialFunctions.cpp, initn2). This encounters a singularity if k'.-k' = -1 exactly.
  double tol = 1e-8;
  // The integrated |k'| is parametrised as k*r
  IOW iow;
        kv[0] = k;
        kv[1] = k*r;
        kv[2] = kv[1];
        xv[0] = -1. + tol;
        xv[1] = x;
        xv[2] = -x;
        kargs[0] = sqrt(kv[1]*kv[1]+kv[0]*kv[0]-2.*kv[1]*kv[0]*xv[1]);
        kargs[2] = sqrt(kv[1]*kv[1]+2.*kv[1]*kv[2]*xv[0]+kv[2]*kv[2]);
        kargs[1] = sqrt(kv[2]*kv[2]+2.*kv[2]*kv[0]*xv[2]+kv[0]*kv[0]);
        kargs[3] = sqrt(kv[0]*kv[0]+2.*kv[0]*kv[1]*xv[1]+kv[1]*kv[1]);
        iow.initn2_pseudo(vars[0],kv,xv,kargs,vars[1],vars[2],vars[3],vars[4],vars[6],model);
        p22 = pow2(F2_nk);
        p13 = F1_nk*F3_nk;

	// for full pseudo we need to multiply by the modified growth / LCDM growth
        double normalisation = pow4(F1_nk);

        iow.initn_lin(vars[0], k*r, vars[1],vars[2], vars[3],vars[4],model);
        double pnorm = pow2(F1_nk);
        iow.initn_lin(vars[0], k, vars[1],vars[2], vars[3],vars[4],model);
        double knorm = pow2(F1_nk);
        iow.initn_lin(vars[0], kargs[0], vars[1],vars[2], vars[3],vars[4],model);
        double kmpnorm = pow2(F1_nk);

    return pow2(r)*2.*P_L(k*r)*pnorm*(P_L(kargs[0])*p22*kmpnorm + 3.*P_L(k)*p13*knorm) / normalisation;
}


// Choose a {0,...,4}: P_linear, P_dd,P_dt, P_tt (MG), P_dd pseudo (see HALO.cpp)
// vars: 0 =  scale factor, 1= omega_m(z=0), 2 = mg param , 3 = mg param, 4 = mg param,
double SPT::PLOOPn2(int a, double vars[], int model, double k, double err) const{
  IOW iow;
double loop, tree;
double prefac = k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt);
real KMAX = QMAXp/k;
real KMIN = QMINp/k;
double c[2] = {KMIN,-1.};
double d[2] = {KMAX, 1.};
switch (a) {
  case 0:
    iow.initn_lin(vars[0], k, vars[1],vars[2], vars[3],vars[4],model);
    tree = pow2(F1_nk/dnorm_spt)*P_L(k);
    return tree;
  case 1:
    loop = prefac*Integrate<2>(bind(ploopn2_mgdd,cref(P_L),vars,model,k,std::placeholders::_1,std::placeholders::_2), c, d, err);
    tree = pow2(F1_nk/dnorm_spt)*P_L(k);
    return loop+tree;
  case 2:
    loop = prefac*Integrate<2>(bind(ploopn2_mgdt,cref(P_L),vars,model,k,std::placeholders::_1,std::placeholders::_2), c, d, err,1e-2);
    tree = (G1_nk*F1_nk)/pow2(dnorm_spt)*P_L(k);
    return loop+tree;
  case 3:
    loop = prefac*Integrate<2>(bind(ploopn2_mgtt,cref(P_L),vars,model,k,std::placeholders::_1,std::placeholders::_2), c, d, err,1e-2);
    tree = pow2(G1_nk/dnorm_spt)*P_L(k);
    return loop+tree;
  case 4:
    loop = prefac*Integrate<2>(bind(ploopn2_mgdd_pseudo,cref(P_L),vars,model,k,std::placeholders::_1,std::placeholders::_2), c, d, err);
    iow.initn_lin(vars[0], k, vars[1],vars[2], vars[3],vars[4],model);
    tree = pow2(F1_nk/dnorm_spt)*P_L(k);
    return loop+tree;
    default:
    warning("SPT: invalid indices, a = %d \n", a);
        return 0;
}}



/* Massive neutrinos 1-loop real and pseudo spectra as described in https://arxiv.org/abs/1902.10692*/


static double ploopn2_mgdd_nu( const PowerSpectrum& P_L, double vars[], int model, double k, double r, double x){
  double kargs[4],kv[3],xv[3], p22,p13;
  // tolerance for ode solver (see SpecialFunctions.cpp, initn2). This encounters a singularity if k'.-k' = -1 exactly.
  double tol = 1e-8;
  // The integrated |k'| is parametrised as k*r
  IOW iow;
        kv[0] = k;
        kv[1] = k*r;
        kv[2] = kv[1];
        xv[0] = -1. + tol;
        xv[1] = x;
        xv[2] = -x;
        kargs[0] = sqrt(kv[1]*kv[1]+kv[0]*kv[0]-2.*kv[1]*kv[0]*xv[1]);
        kargs[2] = sqrt(kv[1]*kv[1]+2.*kv[1]*kv[2]*xv[0]+kv[2]*kv[2]);
        kargs[1] = sqrt(kv[2]*kv[2]+2.*kv[2]*kv[0]*xv[2]+kv[0]*kv[0]);
        kargs[3] = sqrt(kv[0]*kv[0]+2.*kv[0]*kv[1]*xv[1]+kv[1]*kv[1]);
        iow.initn2(vars[0],kv,xv,kargs,vars[1],vars[2],vars[3],vars[4],vars[6],model);
        p22 = pow2(F2_nk);
        p13 = F3_nk;

    return pow2(r)*2.*(P_L(k*r)/pow2(F1p_nk))*( (P_L(kargs[0])/pow2(F1kmp_nk))*p22 + 3.*(P_L(k)/F1_nk)*p13 );
}


static double ploopn2_mgdd_pseudo_nu( const PowerSpectrum& P_L, double vars[], int model, double k, double r, double x){
  double kargs[4],kv[3],xv[3], p22,p13;
  // tolerance for ode solver (see SpecialFunctions.cpp, initn2). This encounters a singularity if k'.-k' = -1 exactly.
  double tol = 1e-8;
  // The integrated |k'| is parametrised as k*r
  IOW iow;
        kv[0] = k;
        kv[1] = k*r;
        kv[2] = kv[1];
        xv[0] = -1. + tol;
        xv[1] = x;
        xv[2] = -x;
        kargs[0] = sqrt(kv[1]*kv[1]+kv[0]*kv[0]-2.*kv[1]*kv[0]*xv[1]);
        kargs[2] = sqrt(kv[1]*kv[1]+2.*kv[1]*kv[2]*xv[0]+kv[2]*kv[2]);
        kargs[1] = sqrt(kv[2]*kv[2]+2.*kv[2]*kv[0]*xv[2]+kv[0]*kv[0]);
        kargs[3] = sqrt(kv[0]*kv[0]+2.*kv[0]*kv[1]*xv[1]+kv[1]*kv[1]);
	// can replace this initialisation to the F2 and F3 EdS kernels for quicker run (only in full pseudo computation, not in unscreened approx)
        iow.initn2_pseudo(vars[0],kv,xv,kargs,vars[1],vars[2],vars[3],vars[4],vars[6],model);
        p22 = pow2(F2_nk);
        p13 = F3_nk;

    return pow2(r)*2.*(P_L(k*r)/pow2(F1p_nk))*( (P_L(kargs[0])/pow2(F1kmp_nk))*p22 + 3.*(P_L(k)/F1_nk)*p13 );
}

// Same as PLOOPn2 but for massive neutrino case
double SPT::PLOOPn2_nu(int a, double vars[], int model, double k, double err) const{
  IOW iow;
double loop, tree;
double prefac = k*k*k/(4*M_PI*M_PI);
real KMAX = QMAXp/k;
real KMIN = QMINp/k;
double c[2] = {KMIN,-1.};
double d[2] = {KMAX, 1.};
switch (a) {
  case 0:
    iow.initn_lin(vars[0], k, vars[1],vars[2], vars[3],vars[4],model);
    tree = P_L(k);
    return tree;
  case 1:
    loop = prefac*Integrate<2>(bind(ploopn2_mgdd_nu,cref(P_L),vars,model,k,std::placeholders::_1,std::placeholders::_2), c, d, err);
    tree = P_L(k);
    return loop+tree;
  case 2:
    loop = prefac*Integrate<2>(bind(ploopn2_mgdd_pseudo_nu,cref(P_L),vars,model,k,std::placeholders::_1,std::placeholders::_2), c, d, err);
    tree = P_L(k);
    return loop+tree;
    default:
    warning("SPT: invalid indices, a = %d \n", a);
        return 0;
}}



// P_loop(k0,z) array initialisation for real and pseudo spectra used in reactions: 1812.05594 (HALO.cpp)
// The pseudo spectrum simply omits screening terms (gamma_2 = gamma_3 = 0) (1606.02520)
// ploopr - real spectrum array
// ploopp - pseudo spectrum array
// redshifts - list of redshifts to compute spectra at
// noz - number of redshifts
// vars:  1 = omega0,  2 = mg param
// k0 - scale at which to compute Spectra

void SPT::ploop_init(double ploopr[], double ploopp[], double redshifts[], int noz, double vars[], int model, double k0){
  IOW iow;

  double res[noz],resp[noz]; // stores loop integral results for real and pseudo spectra
  const int loop_N = 50; // steps in trap rule for loop integral in integrated wave vector p. 50 ensures sub percent accuracy but may be optimised.
  double KMAX = QMAXp/k0; // max r = p/k (QMAXp = 30 - see SpecialFunctions.h)
  double KMIN = QMINp/k0; // min r = p/k (QMINp = 1e-4 - ...)
  double myresult[noz][loop_N];  // stores result for angular integration
  double myresultp[noz][loop_N]; // stores pseudo result
  double k2val[loop_N]; // k array
  double mykernelarray[noz][20]; // SPT kernels stored to this array

  // tolerance for ode solver (see SpecialFunctions.cpp, initn2). This encounters a singularity if k'.-k' = -1 exactly.
  double tol = 1e-8;

// set these to zero because of issue with cosmosis producing nans
  for(int zi=0; zi<noz; zi++){
     res[zi] = 0;
     resp[zi]=0;
  for(int ki=0;ki<loop_N; ki++){
	myresult[zi][ki]=0.;
        myresultp[zi][ki]=0.;
     }}

 // 2D loop integrals
 for(int k2i = 0; k2i< loop_N; k2i++){

        k2val[k2i] = KMIN * exp(k2i*log(KMAX/KMIN)/(loop_N*1.-1.));

        double kargs[4],kv[3],xv[3], temp_ps[noz][32], temp_psp[noz][32]; // some intermediary arrays

//integrate over angle using GL quadrature
	for( int i = 0; i < 32; i++ ){
          kv[0] = k0;
          kv[1] = k0*k2val[k2i];
          kv[2] = kv[1];
          xv[0] = -1. + tol;
          xv[1] = x32[i];
          xv[2] = -x32[i];
          kargs[0] = sqrt(kv[1]*kv[1]+kv[0]*kv[0]-2.*kv[1]*kv[0]*xv[1]);
          kargs[2] = sqrt(kv[1]*kv[1]+2.*kv[1]*kv[2]*xv[0]+kv[2]*kv[2]);
          kargs[1] = sqrt(kv[2]*kv[2]+2.*kv[2]*kv[0]*xv[2]+kv[0]*kv[0]);
          kargs[3] = sqrt(kv[0]*kv[0]+2.*kv[0]*kv[1]*xv[1]+kv[1]*kv[1]);
          // solve differential equations numerically
          iow.initn3(redshifts, noz, kv,xv,kargs,vars[1],vars[2],vars[3],vars[4],mykernelarray,model);

 	  double p22[noz],p13[noz], p22p[noz],p13p[noz];
          // assign all redshifts
          for(int zi = 0; zi<noz; zi++){
            p22[zi] = pow2(mykernelarray[zi][2]);
            p13[zi] = mykernelarray[zi][0]*mykernelarray[zi][4];
            p22p[zi] = pow2(mykernelarray[zi][8]);
            p13p[zi] = mykernelarray[zi][6]*mykernelarray[zi][10];
	   if(kargs[0]<1e-4){
		kargs[0] =1e-4;
	    }
	  temp_ps[zi][i] = pow2(k2val[k2i])*2.*P_L(kv[1])*(P_L(kargs[0])*p22[zi] + 3.*P_L(k0)*p13[zi]);
          temp_psp[zi][i] = pow2(k2val[k2i])*2.*P_L(kv[1])*(P_L(kargs[0])*p22p[zi] + 3.*P_L(k0)*p13p[zi]);
	 }

  // save angular integral per redshift
        for(int zii = 0; zii<noz; zii++){
        myresult[zii][k2i] += w32[i] * temp_ps[zii][i];
        myresultp[zii][k2i] += w32[i] * temp_psp[zii][i];
		}
	      }
	        }

// integrate over r = k2/k using trap rule
for(int zi = 0; zi<noz; zi++){
  for( int i = 1; i < loop_N; i++){
       res[zi] += 0.5 * (k2val[i] - k2val[i-1])*(myresult[zi][i] + myresult[zi][i-1]);
       resp[zi] += 0.5 * (k2val[i] - k2val[i-1])*(myresultp[zi][i] + myresultp[zi][i-1]);
	 }
      }

// assign values of 1-loop spectra to arrays
  for(int zi=0; zi<noz; zi++){
    double tree = pow2(mykernelarray[zi][0]/dnorm_spt)*P_L(k0);
    double loop = pow3(k0)/(4*M_PI*M_PI) * res[zi]/pow4(dnorm_spt);
    double loopp = pow3(k0)/(4*M_PI*M_PI) * resp[zi]/pow4(dnorm_spt);
    ploopr[zi] = tree + loop;
    ploopp[zi] = tree + loopp;
}
}



// P_loop(k0,z) array initialisation for real and pseudo spectra used in reactions including massive neutrinos: 2105.12114
// pkz0 is the P_{L,cb}(k,z) at 'redshifts' and for all integrated k0
// pkz1 is the P_{L,cb} at 'redshifts' and for all integrated k1
// pkz2 is the P_{L,cb} at 'redshifts' and for all integrated k2
// pkz0p is the P_{L,m}(k,z) at 'redshifts' and for all integrated k0
// pkz1p is the P_{L,m} at 'redshifts' and for all integrated k1
// pkz2p is the P_{L,m} at 'redshifts' and for all integrated k2
// ploopr - real spectrum array to be filled by ploop_init_nu
// ploopp - pseudo spectrum array to be filled by ploop_init_nu
// redshifts - list of redshifts to compute spectra at
// noz - number of redshifts
// vars:  1 = omega0,  2 = mg param
// k0 - scale at which to compute Spectra


void SPT::ploop_init_nu(double pkz0[], double pkz1[][50], double pkz2[][50*32], double pkz0p[], double pkz1p[][50], double pkz2p[][50*32], double ploopr[], double ploopp[], double redshifts[], int noz, double vars[], int model, double k0){
  IOW iow;

  double res[noz],resp[noz]; // stores loop integral results for real and pseudo spectra
  const int loop_N = 50; // steps in trap rule for loop integral in integrated wave vector p. 50 ensures sub percent accuracy but may be optimised.
  double KMAX = QMAXp/k0; // max r = p/k (QMAXp = 30 - see SpecialFunctions.h)
  double KMIN = QMINp/k0; // min r = p/k (QMINp = 1e-4 - ...)
  double myresult[noz][loop_N];  // stores result for angular integration
  double myresultp[noz][loop_N]; // stores pseudo result
  double k2val[loop_N]; // k array
  double mykernelarray[noz][20]; // SPT kernels stored to this array

  // tolerance for ode solver (see SpecialFunctions.cpp, initn2). This encounters a singularity if k'.-k' = -1 exactly.
  double tol = 1e-8;

// set these to zero because of issue with cosmosis producing nans
  for(int zi=0; zi<noz; zi++){
     res[zi] = 0;
     resp[zi]=0;
  for(int ki=0;ki<loop_N; ki++){
	myresult[zi][ki]=0.;
        myresultp[zi][ki]=0.;
     }}

 // 2D loop integrals
 for(int k2i = 0; k2i< loop_N; k2i++){

        k2val[k2i] = KMIN * exp(k2i*log(KMAX/KMIN)/(loop_N*1.-1.));

        double kargs[4],kv[3],xv[3], temp_ps[noz][32], temp_psp[noz][32]; // some intermediary arrays

//integrate over angle using GL quadrature
	for( int i = 0; i < 32; i++ ){
          kv[0] = k0;
          kv[1] = k0*k2val[k2i];
          kv[2] = kv[1];
          xv[0] = -1. + tol;
          xv[1] = x32[i];
          xv[2] = -x32[i];
          kargs[0] = sqrt(kv[1]*kv[1]+kv[0]*kv[0]-2.*kv[1]*kv[0]*xv[1]);
          kargs[2] = sqrt(kv[1]*kv[1]+2.*kv[1]*kv[2]*xv[0]+kv[2]*kv[2]);
          kargs[1] = sqrt(kv[2]*kv[2]+2.*kv[2]*kv[0]*xv[2]+kv[0]*kv[0]);
          kargs[3] = sqrt(kv[0]*kv[0]+2.*kv[0]*kv[1]*xv[1]+kv[1]*kv[1]);
          // solve differential equations numerically
          iow.initn3(redshifts, noz, kv,xv,kargs,vars[1],vars[2],vars[3],vars[4],mykernelarray, vars[6], model);

        //  printf("%d \n", i);

 	  double p22[noz],p13[noz], p22p[noz],p13p[noz];
          // assign all redshifts
          for(int zi = 0; zi<noz; zi++){
            p22[zi] = pow2(mykernelarray[zi][2]); // mykernelarray[zi][2] is F2
            p13[zi] = mykernelarray[zi][0]*mykernelarray[zi][4]; // mykernelarray[zi][0] is F1, mykernelarray[zi][4] is F3
            p22p[zi] = pow2(mykernelarray[zi][8]); // mykernelarray[zi][8] is F2_noscr
            p13p[zi] = mykernelarray[zi][6]*mykernelarray[zi][10]; // mykernelarray[zi][6] is F1_noscr, mykernelarray[zi][10] is F3_noscr

    if(kargs[0]<1e-4){
		    kargs[0] =1e-4;
	    }
	      temp_ps[zi][i] = pow2(k2val[k2i])*2.*(pkz1[zi][k2i]/pow2(mykernelarray[zi][14]))*((pkz2[zi][i*loop_N+k2i]/pow2(mykernelarray[zi][12]))*p22[zi] + 3.*(pkz0[zi]/pow2(mykernelarray[zi][0]))*p13[zi]);
        temp_psp[zi][i] = pow2(k2val[k2i])*2.*(pkz1p[zi][k2i]/pow2(mykernelarray[zi][18]))*((pkz2p[zi][i*loop_N+k2i]/pow2(mykernelarray[zi][16]))*p22p[zi] + 3.*(pkz0p[zi]/pow2(mykernelarray[zi][6]))*p13p[zi]);
	 }

  // save angular integral per redshift
        for(int zii = 0; zii<noz; zii++){
        myresult[zii][k2i] += w32[i] * temp_ps[zii][i];
        myresultp[zii][k2i] += w32[i] * temp_psp[zii][i];
		}
	      }
	        }

// integrate over r = k2/k using trap rule
for(int zi = 0; zi<noz; zi++){
  for( int i = 1; i < loop_N; i++){
       res[zi] += 0.5 * (k2val[i] - k2val[i-1])*(myresult[zi][i] + myresult[zi][i-1]);
       resp[zi] += 0.5 * (k2val[i] - k2val[i-1])*(myresultp[zi][i] + myresultp[zi][i-1]);
	 }
      }

  // assign values of 1-loop spectra to arrays
    for(int zi=0; zi<noz; zi++){
      double tree = pkz0[zi];
      double treep = pkz0p[zi];
      double loop = pow3(k0)/(4*M_PI*M_PI) * res[zi];
      double loopp = pow3(k0)/(4*M_PI*M_PI) * resp[zi];
      ploopr[zi] = tree + loop;
      ploopp[zi] = treep + loopp;
  }
  }


////////////////////////////////////////////////
///////HALO-FIT FOR NON-SCALE DEP MODELS : takashi et al 2012 fits  1208.2701///////////
///////////////////////////////////////////////

double wintcambs[3];
static void wintcamb(const PowerSpectrum& P_L, double r){
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

// Initialises components

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

// halofit spectrum
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





////////////////////////////////////////////////
///////REDSHIFT SPACE MODEL-TNS Correction TERMS///////////
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

  	// C terms (DGP 'C' not TNS C)
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

//Atsushi limit on k-angular integration
inline double ATS(real k, real r){
 real KMIN = QMINp/k;
 if(r>=0.5){
  return 1./(2.*r);
 }
 else {
 return  Min(0.999999999, (1.+r*r-KMIN*KMIN)/2./r);
 }
}

static real midintabc(double u0[],real bl, int a, const PowerSpectrum& P_L, real k, real r) {
  real KMAX = QMAXp/k;
	real YMIN = Max(-1., (1.+r*r-KMAX*KMAX)/2./r);
  real YMAX = ATS(k,r);

	real F0,D1,X1;
			F0 = rempf*fdgp_spt;
			D1 = D_spt;
			X1 = 1.;
      real u01 = u0[0];
      real u02 = u0[1];
      real u03 = u0[2];
      real u04 = u0[3];

	//A term
		return   pow4(D1/dnorm_spt)*2.*(F0*u01*pow2(bl)*(Integrate(bind(ABC, 1, cref(P_L), k, r, std::placeholders::_1),  YMIN, YMAX ,  1e-3*P_L(k))
										   + Integrate(bind(ABC, 5, cref(P_L), k,  r, std::placeholders::_1),  YMIN, YMAX ,  1e-3*P_L(k))
										   + Integrate(bind(ABC, 10, cref(P_L), k, r, std::placeholders::_1),  YMIN, YMAX ,  1e-3*P_L(k)))

							+ F0*F0*u01*bl*(Integrate(bind(ABC, 2, cref(P_L), k,  r, std::placeholders::_1),  YMIN, YMAX , 1e-3*P_L(k))
										   + Integrate(bind(ABC, 6, cref(P_L), k,  r, std::placeholders::_1), YMIN, YMAX ,  1e-3*P_L(k))
										   + Integrate(bind(ABC, 11, cref(P_L), k,  r, std::placeholders::_1), YMIN, YMAX ,  1e-3*P_L(k)))

							+  F0*F0*u02*bl*(Integrate(bind(ABC, 3, cref(P_L), k,  r, std::placeholders::_1),   YMIN, YMAX , 1e-3*P_L(k))
										   + Integrate(bind(ABC, 7, cref(P_L), k,  r, std::placeholders::_1),   YMIN, YMAX ,  1e-3*P_L(k))
										   + Integrate(bind(ABC, 12, cref(P_L), k,  r, std::placeholders::_1),  YMIN, YMAX ,  1e-3*P_L(k)))

					     + F0*F0*F0*u02*(Integrate(bind(ABC, 2, cref(P_L), k,  r, std::placeholders::_1),  YMIN, YMAX ,  1e-3*P_L(k))
										   + Integrate(bind(ABC, 8, cref(P_L), k,  r, std::placeholders::_1),   YMIN, YMAX ,  1e-3*P_L(k))
										   + Integrate(bind(ABC, 11, cref(P_L), k, r, std::placeholders::_1), YMIN, YMAX ,  1e-3*P_L(k)))

						   + F0*F0*F0*u03*(Integrate(bind(ABC, 4, cref(P_L), k, r, std::placeholders::_1),   YMIN, YMAX , 1e-3*P_L(k))
										   + Integrate(bind(ABC, 9, cref(P_L), k,  r, std::placeholders::_1),  YMIN, YMAX ,  1e-3*P_L(k))
										   + Integrate(bind(ABC, 13, cref(P_L), k,  r, std::placeholders::_1),YMIN, YMAX ,  1e-3*P_L(k))))

	//B term
		  +      pow4(D1/dnorm_spt)*2.*(F0*F0*pow2(bl)*u01*Integrate(bind(ABC, 14,  cref(P_L), k, r, std::placeholders::_1), YMIN, YMAX , 1e-3* P_L(k))

						  -F0*F0*F0*bl*u01*(Integrate(bind(ABC, 15, cref(P_L), k, r, std::placeholders::_1), YMIN, YMAX , 1e-3*P_L(k))
											    +Integrate(bind(ABC, 16, cref(P_L),  k,  r, std::placeholders::_1),YMIN, YMAX , 1e-3*P_L(k)))

						+F0*F0*F0*F0*u01*Integrate(bind(ABC, 17, cref(P_L),  k, r, std::placeholders::_1), YMIN, YMAX , 1e-3*P_L(k))

							  +F0*F0*pow2(bl)*u02*Integrate(bind(ABC, 18,  cref(P_L), k,  r, std::placeholders::_1),YMIN, YMAX ,  1e-3*P_L(k))

						 -F0*F0*F0*bl*u02*(Integrate(bind(ABC, 19, cref(P_L), k,  r, std::placeholders::_1), YMIN, YMAX , 1e-3*P_L(k))
										      + Integrate(bind(ABC, 20, cref(P_L),  k,  r, std::placeholders::_1), YMIN, YMAX ,  1e-3*P_L(k)))

					  +F0*F0*F0*F0*u02*Integrate(bind(ABC, 21, cref(P_L),  k,  r, std::placeholders::_1),YMIN, YMAX , 1e-3*P_L(k))

					    - F0*F0*F0*bl*u03*(Integrate(bind(ABC, 22,  cref(P_L), k,  r, std::placeholders::_1), YMIN, YMAX ,  1e-3*P_L(k))
											  + Integrate(bind(ABC, 23,  cref(P_L), k,  r, std::placeholders::_1), YMIN, YMAX , 1e-3*P_L(k)))

					  + F0*F0*F0*F0*u03*Integrate(bind(ABC, 24,  cref(P_L), k,  r, std::placeholders::_1),YMIN, YMAX ,  1e-3*P_L(k))

							+F0*F0*F0*F0*u04*Integrate(bind(ABC, 25,  cref(P_L), k,  r, std::placeholders::_1), YMIN, YMAX , 1e-3*P_L(k)))

	//C - nDGP analytic

		  +       X1*D1*D1*2./pow4(dnorm_spt)*(F0*pow2(bl)*u01*(Integrate(bind(ABC, 26, cref(P_L), k, r, std::placeholders::_1),   YMIN, YMAX ,  1e-3*P_L(k))
										+ Integrate(bind(ABC, 31, cref(P_L), k, r, std::placeholders::_1),   YMIN, YMAX ,  1e-3*P_L(k))
										+ Integrate(bind(ABC, 36, cref(P_L), k, r, std::placeholders::_1), YMIN, YMAX ,  1e-3*P_L(k))) +

							F0*F0*bl*u01*(Integrate(bind(ABC, 27, cref(P_L), k, r, std::placeholders::_1),   YMIN, YMAX ,  1e-3*P_L(k))
										+ Integrate(bind(ABC, 32, cref(P_L), k, r, std::placeholders::_1),   YMIN, YMAX ,  1e-3*P_L(k))
										+ Integrate(bind(ABC, 37, cref(P_L), k, r, std::placeholders::_1),  YMIN, YMAX ,  1e-3*P_L(k))) +

						   F0*F0*bl*u02*(Integrate(bind(ABC, 28, cref(P_L), k, r, std::placeholders::_1),   YMIN, YMAX ,  1e-3*P_L(k))
										+ Integrate(bind(ABC, 33, cref(P_L), k, r, std::placeholders::_1),   YMIN, YMAX ,  1e-3*P_L(k))
										+ Integrate(bind(ABC, 38, cref(P_L), k,  r, std::placeholders::_1), YMIN, YMAX ,  1e-3*P_L(k))) +

						F0*F0*F0*u02*(Integrate(bind(ABC, 29, cref(P_L), k, r, std::placeholders::_1),   YMIN, YMAX ,  1e-3*P_L(k))
										+ Integrate(bind(ABC, 34, cref(P_L), k, r, std::placeholders::_1),   YMIN, YMAX ,  1e-3*P_L(k))
										+ Integrate(bind(ABC, 37, cref(P_L), k,  r, std::placeholders::_1),  YMIN, YMAX ,  1e-3*P_L(k))) +

					   F0*F0*F0*u03*(Integrate(bind(ABC, 30, cref(P_L), k, r, std::placeholders::_1),   YMIN, YMAX ,  1e-3*P_L(k))
									   + Integrate(bind(ABC, 35, cref(P_L), k, r, std::placeholders::_1),   YMIN, YMAX ,  1e-3*P_L(k))
									   + Integrate(bind(ABC, 39, cref(P_L), k,  r, std::placeholders::_1),  YMIN, YMAX ,  1e-3*P_L(k))));

}


/* A and B correction term */

/*Analytical A and B*/
// a = 1 : Monopole
// a = 2 : Quadrupole
// a = 3 : Hexdecapole

// For A or B or C seperated : just remove as necessary.
real SPT::AB(real k, real bl, real sigmav, int a) const{
      int n3 = 300;
      real KMAX = QMAXp/k;
      real KMIN = QMINp/k;
      double y[n3];
      double integrand[n3];
      double u0x[4];
      u0x[0]=factL(k, sigmav, 1., 1., 1, a, 7);
      u0x[1]=factL(k, sigmav, 1., 1., 2, a, 7);
      u0x[2]=factL(k, sigmav, 1., 1., 3, a, 7);
      u0x[3]=factL(k, sigmav, 1., 1., 4, a, 7);

      for (int i = 0; i<n3; i++){
      y[i] = KMIN * exp(i*log(KMAX/KMIN)/(n3-1.));
      integrand[i] = midintabc(u0x, bl, a ,cref(P_L), k, y[i]);
      }
    double res = 0.;
      for( int i = 1; i < n3; ++ i ){
    res += 0.5 * (y[i] - y[i-1])*(integrand[i] + integrand[i-1]);
        }
      return  k*k*k/(4*M_PI*M_PI) * res;
      }


real SPT::AB_mu(real k, real bl, real u, int a) const{
      int n3 = 300;
      real KMAX = QMAXp/k;
      real KMIN = QMINp/k;
      double y[n3];
      double integrand[n3];
      double u0x[4];
      u0x[0]= pow2(u);
      u0x[1]= pow2(u0x[0]);
      u0x[2]= u0x[0]*u0x[1];
      u0x[3]= u0x[2]*u0x[2];
      for (int i = 0; i<n3; i++){
      y[i] = KMIN * exp(i*log(KMAX/KMIN)/(n3-1.));
      integrand[i] = midintabc(u0x, bl, a ,cref(P_L), k, y[i]);
      }
    double res = 0.;
      for( int i = 1; i < n3; ++ i ){
    res += 0.5 * (y[i] - y[i-1])*(integrand[i] + integrand[i-1]);
        }
      return  k*k*k/(4*M_PI*M_PI) * res;
      }




/*NUMERICAL TNS WITH VARIOUS BIAS PRESCRIPTIONS */


//// Integrands for numerical TNS //////


// Q-bias
static double ptns_qb( const PowerSpectrum& P_L, double u0[], double vars[], int model, double bk, double k, double r, double x){
  double kargs[4],kv[3],xv[3], abct, pdd, pdd22, pdd13, pdt, pdt22, pdt13, ptt, ptt22, ptt13, prefac, plk, plkmp, plkr,d,r2,r3,x2,x3,x4;
  prefac = k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt);
  d = 1.+ r*r - 2.*r*x;
  // tolerance for ode solver (see SpecialFunctions.cpp, initn2). This encounters a singularity if k'.-k' = -1 exactly.
  double tol = 1e-8;


  if(d < 1e-5){
      return 0;
    }

    // The integrated |k'| is parametrised as k*r
  IOW iow;
        kv[0] = k;
        kv[1] = k*r;
        kv[2] = kv[1];
        xv[0] = -1. + tol ;
        xv[1] = x;
        xv[2] = -x;
        kargs[0] = sqrt(kv[1]*kv[1]+kv[0]*kv[0]-2.*kv[1]*kv[0]*xv[1]);
        kargs[2] = sqrt(kv[1]*kv[1]+2.*kv[1]*kv[2]*xv[0]+kv[2]*kv[2]);
        kargs[1] = sqrt(kv[2]*kv[2]+2.*kv[2]*kv[0]*xv[2]+kv[0]*kv[0]);
        kargs[3] = sqrt(kv[0]*kv[0]+2.*kv[0]*kv[1]*xv[1]+kv[1]*kv[1]);
        iow.initn_rsd(vars[0],kv,xv,kargs,vars[1],vars[2],vars[3],vars[4],model);

        pdd22 = pow2(F2_nk);
        pdd13 = F1_nk*F3_nk;
        pdt22 = G2_nk*F2_nk;
        pdt13 = 0.5*(G1_nk*F3_nk + F1_nk*G3_nk);
        ptt22 = pow2(G2_nk);
        ptt13 = G1_nk * G3_nk;

        plk = P_L(k);
        plkmp = P_L(kargs[0]);
        plkr = P_L(k*r);
        r2 = pow2(r);
        r3 = r2*r;
        x2 = pow2(x);
        x3 = x2*x;
        x4 = x3*x;


// 1-loop spectra
         pdd = bk*bk*u0[4]*(prefac*r2*2.*plkr*(plkmp*pdd22+ 3.*plk*pdd13));
         pdt = bk*u0[0]*(prefac*r2*2.*plkr*(plkmp*pdt22 + 3.*plk*pdt13));
         ptt = u0[1]*(prefac*r2*2.*plkr*(plkmp*ptt22 + 3.*plk*ptt13));


// ab terms
        double abc[10];
        // 1st order kernels F1/G1(k-p)
        abc[0] = F1kmp_nk;
        abc[1] = G1kmp_nk/bk;
        // 1st order kernels F1/G1(p)
        abc[2] = F1p_nk;
        abc[3] = G1p_nk/bk;
        //symmetrized 2nd order kernels for ps F2/G2(p,k-p)
        abc[4] = F2_nk;
        abc[5] = G2_nk/bk;
        //symmetrized 2nd order kernels F2/G2(-p,k)
        abc[6] = F2B_nk;
        abc[7] = G2B_nk/bk;
        //symmetrized 2nd order kernels F2/G2(-k,k-p)
        abc[8] = F2C_nk;
        abc[9] = G2C_nk/bk;


        // A terms
        //u^2
        abct =  pow3(bk)/d*(
                  -u0[0]*(F1_nk*r*(-2.*abc[8]*abc[1]*r*(r*x-1.)+abc[9]*(2.*abc[0]*x*d+abc[1]*r*(1.-x2)))*plk*plkmp // 10% discrep with analytic
                  + abc[4]*r*(-2.*abc[2]*abc[1]*r*(-1.+r*x)+abc[3]*(2.*abc[0]*x*d+abc[1]*r*(1.-x2)))*plkr*plkmp // < 4%
                  + F1_nk*r*(-2.*abc[7]*abc[2]*r*(r*x-1.)+abc[3]*(2.*abc[6]*x*d+abc[7]*r*(1.-x2)))*plk*plkr) // < 3%

        //u^4
        - u0[1]*(r*(2.*abc[8]*G1_nk/bk*abc[1]*r*(r*x-1.)+abc[9]*(G1_nk/bk*(-2.*abc[0]*x*d+abc[1]*r*(x2-1.))+F1_nk*abc[1]*(-2.*x+r*(-1.+3*x2))))*plk*plkmp // < 10% discrep with analytic up to k=0.5
                +r*(abc[4]*abc[1]*abc[3]*(-2.*x+r*(-1.+3.*x2))+abc[5]*(2.*abc[2]*abc[1]*r*(r*x-1.)+abc[3]*(-2.*abc[0]*x*d+abc[1]*r*(x2-1.))))*plkr*plkmp  // < 4% discrep up to k=0.5
                +r*(abc[3]*abc[7]*F1_nk*(-2.*x-r+3.*x2*r)+2.*abc[2]*G1_nk/bk*abc[7]*r*(r*x-1.)+G1_nk/bk*abc[3]*(-2.*abc[6]*x*d+abc[7]*r*(x2-1.)))*plk*plkr) // < 4 % discrep

        //u^6
        - u0[2]*(G1_nk/bk*abc[1]*abc[9]*r*(r+2.*x-3.*r*x2)*plk*plkmp
                  +abc[3]*abc[1]*abc[5]*r*(r+2.*x-3.*r*x2)*plkr*plkmp
                  +abc[3]*G1_nk/bk*abc[7]*r*(r+2.*x-3.*r*x2)*plkr*plk))

        // B terms
        + pow4(bk)*plkr*plkmp/(16.*d*d)*
        //u^2
        (u0[0]*(abc[1]*abc[3]*r2*(x2-1.)*(2.*abc[0]*d*(4.*abc[2]+3.*abc[3]*(x2-1.))
                   +abc[1]*r2*(x2-1.)*(6.*abc[2]+5.*abc[3]*(x2-1.))))
        //u^4
        +u0[1]*(abc[1]*abc[3]*r*(-4.*abc[0]*d*(3.*abc[3]*(x2-1.)*(-2.*x+r*(5.*x2-1.))+abc[2]*(-4*x+r*(-2.+6.*x2)))
                       -3.*abc[1]*r*(x2-1.)*(4.*abc[2]*(2.-6.*r*x+r*r*(-1.+5.*x2))
                       +abc[3]*(x2-1)*(6.-30.*r*x+5.*r*r*(-1.+7.*x2)))))
        //u^6
        +u0[2]*(abc[1]*abc[3]*r*(3.*abc[3]*abc[1]*(x2-1.)*(-8.*x+12.*r*(5.*x2-1.)-20.*r2*x*(7.*x2-3.)+5.*r3*(1.-14.*x2+21.*x4))
                        +2.*abc[0]*abc[3]*(4.*x*(3.-5.*x2)+r3*(3.-30.*x2+35.*x4)+r*(3.-54.*x2+75.*x4)+r2*(6.*x+40.*x3-70.*x4*x))
                        +2.*abc[2]*abc[1]*(-8.*x+12.*r*(3.*x2-1.)+r2*(36.*x-60.*x3)+r3*(3.-30.*x2+35.*x4))))
        //u^8
        +u0[3]*(pow2(abc[1]*abc[3])*r*(8.*x*(5.*x2-3.)-6.*r*(3.-30.*x2+35.*x4)
                      +6.*r2*x*(15.-70.*x2+63.*x4)
                      +r3*(5.-105.*x2+315.*x4-231.*x3*x3))));


      return  pdd - 2.*pdt + ptt + prefac*abct;
    }


//// Lagrangian bias of Roy and Mac Donald

static double ptns_lagb( const PowerSpectrum& P_L, double u0[], double vars[], int model, double bias[], double k, double r, double x){
    double kargs[4],kv[3],xv[3], abct, pdd, pdd22, pdd13, pdt, pdt22, pdt13, ptt, ptt22, ptt13, prefac, plk, plkmp, plkr,d,r2,r3,x2,x3,x4;
    prefac = k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt);
    d = 1.+ r*r - 2.*r*x;
    // tolerance for ode solver (see SpecialFunctions.cpp, initn2). This encounters a singularity if k'.-k' = -1 exactly.
    double tol = 1e-8;

    if(d < 1e-5){
        return 0;
      }

      // The integrated |k'| is parametrised as k*r
    IOW iow;
          kv[0] = k;
          kv[1] = k*r;
          kv[2] = kv[1];
          xv[0] = -1. + tol;
          xv[1] = x;
          xv[2] = -x;
          kargs[0] = sqrt(kv[1]*kv[1]+kv[0]*kv[0]-2.*kv[1]*kv[0]*xv[1]);
          kargs[2] = sqrt(kv[1]*kv[1]+2.*kv[1]*kv[2]*xv[0]+kv[2]*kv[2]);
          kargs[1] = sqrt(kv[2]*kv[2]+2.*kv[2]*kv[0]*xv[2]+kv[0]*kv[0]);
          kargs[3] = sqrt(kv[0]*kv[0]+2.*kv[0]*kv[1]*xv[1]+kv[1]*kv[1]);
          iow.initn_rsd(vars[0],kv,xv,kargs,vars[1],vars[2],vars[3],vars[4],model);

          pdd22 = pow2(F2_nk);
          pdd13 = F1_nk*F3_nk;
          pdt22 = G2_nk*F2_nk;
          pdt13 = 0.5*(G1_nk*F3_nk + F1_nk*G3_nk);
          ptt22 = pow2(G2_nk);
          ptt13 = G1_nk * G3_nk;
          plk = P_L(k);
          plkmp = P_L(kargs[0]);
          plkr = P_L(k*r);
          r2 = pow2(r);
          r3 = r2*r;
          x2 = pow2(x);
          x3 = x2*x;
          x4 = x3*x;

// 1-loop spectra
       pdd = (prefac*r2*2.*plkr*(plkmp*pdd22+ 3.*plk*pdd13));
       pdt = (prefac*r2*2.*plkr*(plkmp*pdt22 + 3.*plk*pdt13));
       ptt = (prefac*r2*2.*plkr*(plkmp*ptt22 + 3.*plk*ptt13));


// Lagrangian bias terms (need to adapt to beyond LCDM)

      double sker = 2./3.-(1.-pow2((x-r)/(kargs[0]/k)));
      double dker = -2./7.*(1.-x2);
      double bs2 = -4./7.*(bias[0]-1.);
      double b3n = 32./315.*(bias[0]-1.);
      double G2ker = G2eds(k*r,kargs[0],(x-r)/(kargs[0]/k)); //
      double F2ker = F2eds(k*r,kargs[0],(x-r)/(kargs[0]/k)); //

      double lagdd =  prefac*pow4(Dl_spt)*((plkr*plkmp*r2*
              (2.*bias[0]*bias[1]*F2ker // b2
              +2.*bias[0]*bs2*F2ker*sker //bs2
              +pow2(bias[1])/2. //b22 a
              +bias[1]*bs2*sker //b2s2 a
              +pow2(bs2)*pow2(sker)/2.) //bs22 a

            - r2*pow2(plkr)*(pow2(bias[1])/2. //b22 b
                            +2.*bias[1]*bs2/3. // b2s2 b
                            +pow2(bs2)*2./9.))// bs22 b

            +plk*b3n*bias[0]*plkr*2.*105./16.*r2*(dker*sker+8./63.)); // b3 terms

      double lagdt = prefac*pow4(Dl_spt)*(-fl_spt*plkr*plkmp*r2*(bias[1]*G2ker //b2t
                                                 +bs2*G2ker*sker) // bs2t
              -fl_spt*plk*b3n*plkr*105./16.*r2*(dker*sker+8./63.));


// ab terms
      double abc[10];
      // 1st order kernels F1/G1(k-p)
      abc[0] = F1kmp_nk;
      abc[1] = G1kmp_nk/bias[0];
      // 1st order kernels F1/G1(p)
      abc[2] = F1p_nk;
      abc[3] = G1p_nk/bias[0];
      //symmetrized 2nd order kernels for ps F2/G2(p,k-p)
      abc[4] = F2_nk;
      abc[5] = G2_nk/bias[0];
      //symmetrized 2nd order kernels F2/G2(-p,k)
      abc[6] = F2B_nk;
      abc[7] = G2B_nk/bias[0];
      //symmetrized 2nd order kernels F2/G2(-k,k-p)
      abc[8] = F2C_nk;
      abc[9] = G2C_nk/bias[0];



      // A terms
      //u^2
      abct =  pow3(bias[0])/d*(
                -u0[0]*(F1_nk*r*(-2.*abc[8]*abc[1]*r*(r*x-1.)+abc[9]*(2.*abc[0]*x*d+abc[1]*r*(1.-x2)))*plk*plkmp // 10% discrep with analytic
                + abc[4]*r*(-2.*abc[2]*abc[1]*r*(-1.+r*x)+abc[3]*(2.*abc[0]*x*d+abc[1]*r*(1.-x2)))*plkr*plkmp // < 4%
                + F1_nk*r*(-2.*abc[7]*abc[2]*r*(r*x-1.)+abc[3]*(2.*abc[6]*x*d+abc[7]*r*(1.-x2)))*plk*plkr) // < 3%

      //u^4
      - u0[1]*(r*(2.*abc[8]*G1_nk/bias[0]*abc[1]*r*(r*x-1.)+abc[9]*(G1_nk/bias[0]*(-2.*abc[0]*x*d+abc[1]*r*(x2-1.))+F1_nk*abc[1]*(-2.*x+r*(-1.+3*x2))))*plk*plkmp // < 10% discrep with analytic up to k=0.5
              +r*(abc[4]*abc[1]*abc[3]*(-2.*x+r*(-1.+3.*x2))+abc[5]*(2.*abc[2]*abc[1]*r*(r*x-1.)+abc[3]*(-2.*abc[0]*x*d+abc[1]*r*(x2-1.))))*plkr*plkmp  // < 4% discrep up to k=0.5
              +r*(abc[3]*abc[7]*F1_nk*(-2.*x-r+3.*x2*r)+2.*abc[2]*G1_nk/bias[0]*abc[7]*r*(r*x-1.)+G1_nk/bias[0]*abc[3]*(-2.*abc[6]*x*d+abc[7]*r*(x2-1.)))*plk*plkr) // < 4 % discrep

      //u^6
      - u0[2]*(G1_nk/bias[0]*abc[1]*abc[9]*r*(r+2.*x-3.*r*x2)*plk*plkmp
                +abc[3]*abc[1]*abc[5]*r*(r+2.*x-3.*r*x2)*plkr*plkmp
                +abc[3]*G1_nk/bias[0]*abc[7]*r*(r+2.*x-3.*r*x2)*plkr*plk))

      // B terms
      + pow4(bias[0])*plkr*plkmp/(16.*d*d)*
      //u^2
      (u0[0]*(abc[1]*abc[3]*r2*(x2-1.)*(2.*abc[0]*d*(4.*abc[2]+3.*abc[3]*(x2-1.))
                 +abc[1]*r2*(x2-1.)*(6.*abc[2]+5.*abc[3]*(x2-1.))))
      //u^4
      +u0[1]*(abc[1]*abc[3]*r*(-4.*abc[0]*d*(3.*abc[3]*(x2-1.)*(-2.*x+r*(5.*x2-1.))+abc[2]*(-4*x+r*(-2.+6.*x2)))
                     -3.*abc[1]*r*(x2-1.)*(4.*abc[2]*(2.-6.*r*x+r*r*(-1.+5.*x2))
                     +abc[3]*(x2-1)*(6.-30.*r*x+5.*r*r*(-1.+7.*x2)))))
      //u^6
      +u0[2]*(abc[1]*abc[3]*r*(3.*abc[3]*abc[1]*(x2-1.)*(-8.*x+12.*r*(5.*x2-1.)-20.*r2*x*(7.*x2-3.)+5.*r3*(1.-14.*x2+21.*x4))
                      +2.*abc[0]*abc[3]*(4.*x*(3.-5.*x2)+r3*(3.-30.*x2+35.*x4)+r*(3.-54.*x2+75.*x4)+r2*(6.*x+40.*x3-70.*x4*x))
                      +2.*abc[2]*abc[1]*(-8.*x+12.*r*(3.*x2-1.)+r2*(36.*x-60.*x3)+r3*(3.-30.*x2+35.*x4))))
      //u^8
      +u0[3]*(pow2(abc[1]*abc[3])*r*(8.*x*(5.*x2-3.)-6.*r*(3.-30.*x2+35.*x4)
                    +6.*r2*x*(15.-70.*x2+63.*x4)
                    +r3*(5.-105.*x2+315.*x4-231.*x3*x3))));



    return u0[4]*(pow2(bias[0])*pdd + lagdd) - 2.*u0[0]*(bias[0]*pdt + lagdt) + u0[1]*ptt + prefac*abct;  // in other code multiply abc by 2

  }


  //// Integrands for numerical RSD SPT 1-loop spectrum //////
  static double pspt( const PowerSpectrum& P_L, double u0[], double vars[], int model, double bk, double k, double r, double x){
    double kargs[4],kv[3],xv[3], abct, pdd, pdd22, pdd13, pdt, pdt22, pdt13, ptt, ptt22, ptt13, prefac, plk, plkmp, plkr;
    double d,d2,r2,r3,x2,x3,x4,r4,x6,x5;
    prefac = k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt);
    d = 1.+ r*r - 2.*r*x;
    // tolerance for ode solver (see SpecialFunctions.cpp, initn2). This encounters a singularity if k'.-k' = -1 exactly.
    double tol = 1e-8;

    if(d < 1e-5){
        return 0;
      }

      // The integrated |k'| is parametrised as k*r
    IOW iow;
          kv[0] = k;
          kv[1] = k*r;
          kv[2] = kv[1];
          xv[0] = -1. + tol ;
          xv[1] = x;
          xv[2] = -x;
          kargs[0] = sqrt(kv[1]*kv[1]+kv[0]*kv[0]-2.*kv[1]*kv[0]*xv[1]);
          kargs[2] = sqrt(kv[1]*kv[1]+2.*kv[1]*kv[2]*xv[0]+kv[2]*kv[2]);
          kargs[1] = sqrt(kv[2]*kv[2]+2.*kv[2]*kv[0]*xv[2]+kv[0]*kv[0]);
          kargs[3] = sqrt(kv[0]*kv[0]+2.*kv[0]*kv[1]*xv[1]+kv[1]*kv[1]);
          iow.initn_rsd(vars[0],kv,xv,kargs,vars[1],vars[2],vars[3],vars[4],model);

          pdd22 = pow2(F2_nk);
          pdd13 = F1_nk*F3_nk;
          pdt22 = G2_nk*F2_nk;
          pdt13 = 0.5*(G1_nk*F3_nk + F1_nk*G3_nk);
          ptt22 = pow2(G2_nk);
          ptt13 = G1_nk * G3_nk;

          plk = P_L(k);
          plkmp = P_L(kargs[0]);
          plkr = P_L(k*r);

          d2 = pow2(d);
          r2 = pow2(r);
          r3 = r2*r;
          r4 = r2*r2;
          x2 = pow2(x);
          x3 = x2*x;
          x4 = x3*x;
          x5 = x4*x;
          x6 = x5*x;

  // 1-loop spectra
           pdd = bk*bk*u0[4]*(prefac*r2*2.*plkr*(plkmp*pdd22+ 3.*plk*pdd13));
           pdt = bk*u0[0]*(prefac*r2*2.*plkr*(plkmp*pdt22 + 3.*plk*pdt13));
           ptt = u0[1]*(prefac*r2*2.*plkr*(plkmp*ptt22 + 3.*plk*ptt13));


  // ab terms
          double abc[10];
          // 1st order kernels F1/G1(k-p)
          abc[0] = F1kmp_nk;
          abc[1] = G1kmp_nk/bk;
          // 1st order kernels F1/G1(p)
          abc[2] = F1p_nk;
          abc[3] = G1p_nk/bk;
          //symmetrized 2nd order kernels for ps F2/G2(p,k-p)
          abc[4] = F2_nk;
          abc[5] = G2_nk/bk;
          //symmetrized 2nd order kernels F2/G2(-p,k)
          abc[6] = F2B_nk;
          abc[7] = G2B_nk/bk;
          //symmetrized 2nd order kernels F2/G2(-k,k-p)
          abc[8] = F2C_nk;
          abc[9] = G2C_nk/bk;


          // A terms
          //u^2
          abct = pow3(bk)/d*(
                    -u0[0]*(F1_nk*r*(-2.*abc[8]*abc[1]*r*(r*x-1.)+abc[9]*(2.*abc[0]*x*d+abc[1]*r*(1.-x2)))*plk*plkmp // 10% discrep with analytic
                    + abc[4]*r*(-2.*abc[2]*abc[1]*r*(-1.+r*x)+abc[3]*(2.*abc[0]*x*d+abc[1]*r*(1.-x2)))*plkr*plkmp // < 4%
                    + F1_nk*r*(-2.*abc[7]*abc[2]*r*(r*x-1.)+abc[3]*(2.*abc[6]*x*d+abc[7]*r*(1.-x2)))*plk*plkr) // < 3%

          //u^4
          - u0[1]*(r*(2.*abc[8]*G1_nk/bk*abc[1]*r*(r*x-1.)+abc[9]*(G1_nk/bk*(-2.*abc[0]*x*d+abc[1]*r*(x2-1.))+F1_nk*abc[1]*(-2.*x+r*(-1.+3*x2))))*plk*plkmp // < 10% discrep with analytic up to k=0.5
                  +r*(abc[4]*abc[1]*abc[3]*(-2.*x+r*(-1.+3.*x2))+abc[5]*(2.*abc[2]*abc[1]*r*(r*x-1.)+abc[3]*(-2.*abc[0]*x*d+abc[1]*r*(x2-1.))))*plkr*plkmp  // < 4% discrep up to k=0.5
                  +r*(abc[3]*abc[7]*F1_nk*(-2.*x-r+3.*x2*r)+2.*abc[2]*G1_nk/bk*abc[7]*r*(r*x-1.)+G1_nk/bk*abc[3]*(-2.*abc[6]*x*d+abc[7]*r*(x2-1.)))*plk*plkr) // < 4 % discrep

          //u^6
          - u0[2]*(G1_nk/bk*abc[1]*abc[9]*r*(r+2.*x-3.*r*x2)*plk*plkmp
                    +abc[3]*abc[1]*abc[5]*r*(r+2.*x-3.*r*x2)*plkr*plkmp
                    +abc[3]*G1_nk/bk*abc[7]*r*(r+2.*x-3.*r*x2)*plkr*plk))

          // B terms
          + pow4(bk)*plkr*plkmp/(16.*d2)*
          //u^2
          (u0[0]*(abc[1]*abc[3]*r2*(x2-1.)*(2.*abc[0]*d*(4.*abc[2]+3.*abc[3]*(x2-1.))
                     +abc[1]*r2*(x2-1.)*(6.*abc[2]+5.*abc[3]*(x2-1.))))
          //u^4
          +u0[1]*(abc[1]*abc[3]*r*(-4.*abc[0]*d*(3.*abc[3]*(x2-1.)*(-2.*x+r*(5.*x2-1.))+abc[2]*(-4*x+r*(-2.+6.*x2)))
                         -3.*abc[1]*r*(x2-1.)*(4.*abc[2]*(2.-6.*r*x+r*r*(-1.+5.*x2))
                         +abc[3]*(x2-1)*(6.-30.*r*x+5.*r*r*(-1.+7.*x2)))))
          //u^6
          +u0[2]*(abc[1]*abc[3]*r*(3.*abc[3]*abc[1]*(x2-1.)*(-8.*x+12.*r*(5.*x2-1.)-20.*r2*x*(7.*x2-3.)+5.*r3*(1.-14.*x2+21.*x4))
                          +2.*abc[0]*abc[3]*(4.*x*(3.-5.*x2)+r3*(3.-30.*x2+35.*x4)+r*(3.-54.*x2+75.*x4)+r2*(6.*x+40.*x3-70.*x5))
                          +2.*abc[2]*abc[1]*(-8.*x+12.*r*(3.*x2-1.)+r2*(36.*x-60.*x3)+r3*(3.-30.*x2+35.*x4))))
          //u^8
          +u0[3]*(pow2(abc[1]*abc[3])*r*(8.*x*(5.*x2-3.)-6.*r*(3.-30.*x2+35.*x4)
                        +6.*r2*x*(15.-70.*x2+63.*x4)
                        +r3*(5.-105.*x2+315.*x4-231.*x6))))

          // C terms
            + pow4(bk)*plkr*plkmp/(16.*d2)*pow2(abc[3])*
            // u^2
            (u0[0]*(8.*pow2(abc[0])*d2*(1.-x2) - 12.*abc[0]*abc[1]*d*(r2-2.*r2*x2+r2*x4)
                    + pow2(abc[1])*(5.*r4 - 15.*r4*x2 + 15.*r4*x4 - 5.*r4*x6))
            // u^4
            +u0[1]*(8.*pow2(abc[0])*d2*(3.*x2-1.) - 4.*abc[0]*abc[1]*d*(4.-6.*r2-24.*r*x-4*x2+36.*r2*x2+24.*r*x3-30.*r2*x4)
                    + pow2(abc[1])*(36.*r2 - 15.*r4 - 120.*r3*x - 72.*r2*x2 + 135.*r4*x2 +
                      240.*r3*x3 + 36.*r2*x4 - 225.*r4*x4 - 120.*r3*x5 + 105.*r4*x6))
            // u^6
            +u0[2]*(-4.*abc[0]*abc[1]*d*(-4. + 3.*r2 + 24.*r*x + 12.*x2 - 30.*r2*x2 - 40.*r*x3 + 35.*r2*x4)
                    + pow2(abc[1])*(8. - 72.*r2 + 15*r4 - 96.*r*x + 240.*r3*x - 8.*x2 + 432.*r2*x2 -
                                      225.*r4*x2 + 96.*r*x3 - 800.*r3*x3 - 360.*r2*x4 + 525.*r4*x4 +
                                      560.*r3*x5 - 315.*r4*x6))
            // u^8
            +u0[3]*pow2(abc[1])*(-8. + 36.*r2 - 5.*r4 + 96.*r*x - 120.*r3*x + 24.*x2 - 360.*r2*x2 +
                                        105.*r4*x2 - 160.*r*x3 + 560.*r3*x3 + 420.*r2*x4 - 315.*r4*x4 -
                                        504.*r3*x5 + 231.*r4*x6));

        return  pdd - 2.*pdt + ptt + prefac*abct;
      }



// Velocity dispersion integrand
static real vel_disp_lin(const PowerSpectrum& P_L, double vars[], int model, double q){
  IOW iow;
  iow.initn_lin(vars[0], q, vars[1],vars[2], vars[3],vars[4],model);
  double growth = pow2(G1_nk/dnorm_spt);
  return  P_L(q)*growth/(6.*pow2(M_PI));
}


// sigma_v^2
real SPT::sigmav_init(double vars[], int model) const{
    return Integrate(bind(vel_disp_lin, cref(P_L), vars, model, std::placeholders::_1), QMINp, QMAXp, epsrel);
  }

/////// REDSHIFT SPACE MODIFIED GRAVITY SPECTRUM ////////
// a {0,..,3}  : 0=Kaiser, 1 = TNS q-bias [1507.01592], 2 = TNS lag bias (MG) - incomplete!, 3 = 1-loop SPT (MG) [see Eq.23 of 1006.0699 for example]
// b {1,2,3} :  1 = monopole, 2 = quadrupole, 3 = hexdecapole
// bias[] :  0 = linear bias, 1,2 = qbias param or lagrangian bias params (b_2, N)
// vars: 0 =  scale factor, 1= omega_m(z=0), 2 = mg param , 3 = mg param, 4 = mg param,
// pars[0] -  sigma_v free parameter for TNS and Kaiser or velocity dispersion
// err -  absolute error in differential equation solver

double SPT::PRSD_mg(int a, int b, double bias[], double vars[], int model, double pars[], double k, double err) const{
IOW iow;
double linear, nonlinear, u0x[8], bk, bl, stoch, kaiser_term, myd2, myf, k2, b2;
real KMAX = QMAXp/k;
real KMIN = QMINp/k;
double c[2] = {KMIN,-0.99999999};
double d[2] = {KMAX, 0.99999999};
switch (a) {
  case 0:
      bl = bias[0];
      iow.initn_lin(vars[0], k, vars[1],vars[2], vars[3], vars[4],model);
      linear = pow2(F1_nk*bl/dnorm_spt)*(factL(k, pars[0], 1., 1., 0, b, 7)*P_L(k) - 2.*(G1_nk/F1_nk/bl)*factL(k, pars[0], 1., 1., 1, b, 7)*P_L(k) + pow2(G1_nk/F1_nk/bl)*factL(k, pars[0], 1., 1., 2, b, 7)* P_L(k));

    return linear;

  case 1:
    bk =  bias[0]*sqrt((1.+bias[1]*pow2(k))/(1.+bias[2]*k));
// multipole prefactors or powers of mu
    u0x[0]=factL(k, pars[0], 1., 1., 1, b, 7);
    u0x[1]=factL(k, pars[0], 1., 1., 2, b, 7);
    u0x[2]=factL(k, pars[0], 1., 1., 3, b, 7);
    u0x[3]=factL(k, pars[0], 1., 1., 4, b, 7);
    u0x[4]=factL(k, pars[0], 1., 1., 0, b, 7);
    nonlinear = Integrate<2>(bind(ptns_qb,cref(P_L), u0x, vars, model, bk, k, std::placeholders::_1,std::placeholders::_2), c, d, err);
    linear = pow2(F1_nk*bias[0]/dnorm_spt)*(u0x[4]*P_L(k) - 2.*(G1_nk/F1_nk/bias[0])*u0x[0]*P_L(k) + pow2(G1_nk/F1_nk/bias[0])*u0x[1]*P_L(k));

    return linear + nonlinear;

  case 2:
// multipole prefactors or powers of mu
    u0x[0]=factL(k, pars[0], 1., 1., 1, b, 7);
    u0x[1]=factL(k, pars[0], 1., 1., 2, b, 7);
    u0x[2]=factL(k, pars[0], 1., 1., 3, b, 7);
    u0x[3]=factL(k, pars[0], 1., 1., 4, b, 7);
    u0x[4]=factL(k, pars[0], 1., 1., 0, b, 7);
    nonlinear = Integrate<2>(bind(ptns_lagb,cref(P_L), u0x, vars, model, bias, k,std::placeholders::_1,std::placeholders::_2), c, d, err);
    linear = pow2(F1_nk*bias[0]/dnorm_spt)*(u0x[4]*P_L(k) - 2.*(G1_nk/F1_nk/bias[0])*u0x[0]*P_L(k) + pow2(G1_nk/F1_nk/bias[0])*u0x[1]*P_L(k));
    stoch =  u0x[4]*bias[2];

    return linear + nonlinear + stoch;

  case 3:
      k2 = pow2(k);
      b2 = pow2(bias[0]);
  // multipole prefactors or powers of mu
      u0x[0] = factL(k,  0., 1., 1., 1 , b,  1); // u2
      u0x[1] = factL(k,  0., 1., 1., 2 , b,  1); // u4
      u0x[2] = factL(k,  0., 1., 1., 3 , b,  1); // u6
      u0x[3] = factL(k,  0., 1., 1., 4 , b,  1); // u8
      u0x[4] = factL(k,  0., 1., 1., 0 , b,  1); // u0
      u0x[5] = factL(k, pars[0]*k2, 1., 1., 0 , b,  8); // u0i : (k sigmav mu)^2  [pars[0] = sigmav]
      u0x[6] = factL(k, pars[0]*k2, 1., 1., 1 , b,  8); // u2i : (k sigmav mu)^2 u^2
      u0x[7] = factL(k, pars[0]*k2, 1., 1., 2 , b,  8); // u4i : (k sigmav mu)^2 u^4

      nonlinear = Integrate<2>(bind(pspt,cref(P_L), u0x, vars, model, bias[0], k, std::placeholders::_1,std::placeholders::_2), c, d, err);
      myd2 = pow2(F1_nk/dnorm_spt);
      myf = -G1_nk/F1_nk;
      linear = (u0x[4]*myd2*pow2(bias[0])*P_L(k) + 2.*u0x[0]*myd2*myf*bias[0]*P_L(k) + u0x[1]*myd2*pow2(myf)*b2*P_L(k));
      kaiser_term = -(u0x[5]*myd2*pow2(bias[0])*P_L(k) + 2.*u0x[6]*myd2*myf*bias[0]*P_L(k) + u0x[7]*myd2*pow2(myf)*b2*P_L(k));


      return linear + nonlinear + kaiser_term;

  default:
  warning("SPT: invalid indices, a = %d \n", a);
      return 0;
    }
  }
