#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/bind.hpp>
using boost::cref;


#include <cfloat>
#include <cmath>

#include "Common.h"
#include "Quadrature.h"
#include "SpecialFunctions.h"

//#include<omp.h>

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_bessel.h>

#include <math.h>       /* pow */
#include <iostream>
#include <stdlib.h>

/*NOTES TO TRAVELERS*/

/* This file is a bit of a mess but if you're looking for the modified gravity functions search for gamma2 or gamma3 (they were last sighted on lines 302-390, beware there be dragons ... )*/
/* the function that solves the evolution equations for these specified gammas is initn2 (or the older initn which I haven't removed yet) - you can search for that too but he's hiding somewhere deep in this jungle */


////////////////////////////////////////////////////////////
// Gauss-Legendre quadrature

const int integral::gl_num = 128;
const double integral::x[ 128 ] =
{-0.9998248879471319,-0.9990774599773759,-0.997733248625514,
 -0.9957927585349812,-0.9932571129002129,-0.9901278184917344,
 -0.9864067427245862,-0.9820961084357185,-0.9771984914639074,
 -0.9717168187471366,-0.9656543664319653,-0.9590147578536999,
 -0.9518019613412644,-0.9440202878302202,-0.9356743882779164,
 -0.9267692508789478,-0.9173101980809605,-0.9073028834017568,
 -0.8967532880491582,-0.8856677173453972,-0.8740527969580318,
 -0.8619154689395485,-0.849262987577969,-0.8361029150609068,
 -0.8224431169556438,-0.8082917575079137,-0.7936572947621933,
 -0.778548475506412,-0.7629743300440947,-0.746944166797062,
 -0.7304675667419088,-0.7135543776835874,-0.6962147083695143,
 -0.6784589224477193,-0.6602976322726461,-0.6417416925623076,
 -0.6228021939105849,-0.6034904561585486,-0.5838180216287631,
 -0.5637966482266181,-0.5434383024128104,-0.5227551520511755,
 -0.5017595591361445,-0.480464072404172,-0.4588814198335522,
 -0.4370245010371042,-0.414906379552275,-0.3925402750332674,
 -0.369939555349859,-0.3471177285976355,-0.3240884350244134,
 -0.3008654388776772,-0.2774626201779044,-0.2538939664226943,
 -0.23017356422666,-0.2063155909020792,-0.1823343059853372,
 -0.1582440427142249,-0.1340591994611878,-0.1097942311276437,
 -0.0854636405045155,-0.06108196960413957,-0.03666379096873349,
 -0.01222369896061576,0.01222369896061576,0.03666379096873349,
 0.06108196960413957,0.0854636405045155,0.1097942311276437,
 0.1340591994611878,0.1582440427142249,0.1823343059853372,
 0.2063155909020792,0.23017356422666,0.2538939664226943,
 0.2774626201779044,0.3008654388776772,0.3240884350244134,
 0.3471177285976355,0.369939555349859,0.3925402750332674,
 0.414906379552275,0.4370245010371042,0.4588814198335522,
 0.480464072404172,0.5017595591361445,0.5227551520511755,
 0.5434383024128104,0.5637966482266181,0.5838180216287631,
 0.6034904561585486,0.6228021939105849,0.6417416925623076,
 0.6602976322726461,0.6784589224477193,0.6962147083695143,
 0.7135543776835874,0.7304675667419088,0.746944166797062,
 0.7629743300440947,0.778548475506412,0.7936572947621933,
 0.8082917575079137,0.8224431169556438,0.8361029150609068,
 0.849262987577969,0.8619154689395485,0.8740527969580318,
 0.8856677173453972,0.8967532880491582,0.9073028834017568,
 0.9173101980809605,0.9267692508789478,0.9356743882779164,
 0.9440202878302202,0.9518019613412644,0.9590147578536999,
 0.9656543664319653,0.9717168187471366,0.9771984914639074,
 0.9820961084357185,0.9864067427245862,0.9901278184917344,
 0.9932571129002129,0.9957927585349812,0.997733248625514,
 0.9990774599773759,0.9998248879471319};
const double integral::w[ 128 ]=
{0.0004493809602920904,0.001045812679340349,
 0.00164250301866903,0.002238288430962619,0.002832751471457991,
 0.003425526040910216,0.004016254983738642,0.004604584256702955,
 0.00519016183267633,0.005772637542865699,0.006351663161707189,
 0.006926892566898814,0.007497981925634729,0.008064589890486058,
 0.00862637779861675,0.009183009871660874,0.009734153415006806,
 0.01027947901583216,0.01081866073950308,0.01135137632408042,
 0.01187730737274028,0.01239613954395092,0.01290756273926735,
 0.01341127128861633,0.01390696413295199,0.01439434500416685,
 0.01487312260214731,0.01534301076886514,0.01580372865939935,
 0.01625500090978519,0.0166965578015892,0.01712813542311138,
 0.0175494758271177,0.01796032718500869,0.01836044393733134,
 0.01874958694054471,0.01912752360995095,0.0194940280587066,
 0.01984888123283086,0.02019187104213004,0.02052279248696007,
 0.02084144778075115,0.02114764646822135,0.02144120553920846,
 0.02172194953805208,0.02198971066846049,0.02224432889379977,
 0.02248565203274497,0.02271353585023646,0.02292784414368685,
 0.02312844882438703,0.02331522999406276,0.02348807601653591,
 0.02364688358444762,0.0237915577810034,0.02392201213670346,
 0.02403816868102405,0.02413995798901928,0.02422731922281525,
 0.02430020016797187,0.02435855726469063,0.02440235563384958,
 0.02443156909785005,0.02444618019626252,0.02444618019626252,
 0.02443156909785005,0.02440235563384958,0.02435855726469063,
 0.02430020016797187,0.02422731922281525,0.02413995798901928,
 0.02403816868102405,0.02392201213670346,0.0237915577810034,
 0.02364688358444762,0.02348807601653591,0.02331522999406276,
 0.02312844882438703,0.02292784414368685,0.02271353585023646,
 0.02248565203274497,0.02224432889379977,0.02198971066846049,
 0.02172194953805208,0.02144120553920846,0.02114764646822135,
 0.02084144778075115,0.02052279248696007,0.02019187104213004,
 0.01984888123283086,0.0194940280587066,0.01912752360995095,
 0.01874958694054471,0.01836044393733134,0.01796032718500869,
 0.0175494758271177,0.01712813542311138,0.0166965578015892,
 0.01625500090978519,0.01580372865939935,0.01534301076886514,
 0.01487312260214731,0.01439434500416685,0.01390696413295199,
 0.01341127128861633,0.01290756273926735,0.01239613954395092,
 0.01187730737274028,0.01135137632408042,0.01081866073950308,
 0.01027947901583216,0.009734153415006806,0.009183009871660874,
 0.00862637779861675,0.008064589890486058,0.007497981925634729,
 0.006926892566898814,0.006351663161707189,0.005772637542865699,
 0.00519016183267633,0.004604584256702955,0.004016254983738642,
 0.003425526040910216,0.002832751471457991,0.002238288430962619,
 0.00164250301866903,0.001045812679340349,0.0004493809602920904};


 /* n = 512 */


//F1[k],G1[k]

double F1_nk;
double G1_nk;

double F1q_nk;
double G1q_nk;

//F1[k-p],G1[k-p]   (or F1/G2[-k-p] for bispectrum)
double F1kmp_nk[(n1+1)*n2];
double G1kmp_nk[(n1+1)*n2];

double F1kmpq_nk[(n1+1)*n2];
double G1kmpq_nk[(n1+1)*n2];

//F1[p], G1[p]
double F1p_nk[(n1+1)*n2];
double G1p_nk[(n1+1)*n2];

double F1pq_nk[(n1+1)*n2];
double G1pq_nk[(n1+1)*n2];


//2nd order used in Bispectrum term : F2/G2[p,k]
double F2A_nk[(n1+1)*n2];
double G2A_nk[(n1+1)*n2];

double F2Aq_nk[(n1+1)*n2];
double G2Aq_nk[(n1+1)*n2];

//2nd order used in Bispectrum term : F2/G2[-p,k]
double F2B_nk[(n1+1)*n2];
double G2B_nk[(n1+1)*n2];


//2nd order used in Bispectrum term : F2/G2[-k,k-p]
double F2C_nk[(n1+1)*n2];
double G2C_nk[(n1+1)*n2];

//2nd order used in S3 term : F2[-p,p] (or used in CDE for F2[p,p] and G2[p,p])
double F2D_nk[(n1+1)*n2];
double G2D_nk[(n1+1)*n2];


double F2Dq_nk[(n1+1)*n2];
double G2Dq_nk[(n1+1)*n2];


//2nd order used in P22 : F2/G2[k-p, p] (or F2/G2[-k-p,p] for bispectrum)
double F2_nk[(n1+1)*n2];
double G2_nk[(n1+1)*n2];

double F2q_nk[(n1+1)*n2];
double G2q_nk[(n1+1)*n2];

//3rd order used in P13 : F3/G3[k,p,-p]
double F3_nk[(n1+1)*n2];
double G3_nk[(n1+1)*n2];

double F3q_nk[(n1+1)*n2];
double G3q_nk[(n1+1)*n2];

  ///// EVOLUTION FACTORS ////////

  // LCDM
double Dl_spt;
  // DGP
double D_spt; //
double F_spt; //
double Cx_spt;
double I_spt;
double J_spt;
double K_spt;
double KL_spt;
double F3_spt;


//derivatives of evolution factors
  //LCDM
double fl_spt;
  //DGP
double fdgp_spt;
double Dd_spt;
double Fd_spt;
double Cdx_spt;
double Id_spt;
double Jd_spt;
double Kd_spt;
double KLd_spt;
double F3d_spt;

  //Hubble
double H_spt;

// Normalisation of linear power spectrum (Dl @ a=1)
double dnorm_spt;

double gnorm_spt;


//#include <gsl/gsl_sf_gamma.h>


//////// SPT EVOLUTION FACTORS AND NUMERICAL KERNEL CONSTRUCTION //////////
/* Power Spectrum angular integral limits require the following functions */

// Minimum and Maximum selection functions
double Min(real a, real b) {
real p;
if(a<b)
p=a;
else
p=b;
return p;
}

double Max(real a, real b) {
real p;
if(a<b)
p=b;
else
p=a;
return p;
}

// Normalized Hubble parameter H and a*H*dH/da for the background (LCDM)

double HA(double a, double omega0){
	double omegaL= 1.-omega0;
	return  sqrt(omega0/pow(a,3)+omegaL);}

double HA1(double a, double omega0){
	return -3.*omega0/(2.*pow(a,3));}

/*STANDARD KERNEL FUNCTIONS */
// Here u1 is the angle between vectors of magnitude k1 and k2.

//1-mu^2
static double ker1(double u1){
	return 1.-u1*u1;
}

//alpha(k1,k2)
double alpha(double k1, double k2, double u1){
	return 1.+k2*u1/k1;
}

//Symmetrized alpha : [alpha(k1,k2)+alpha(k2,k1)]/2
double alphas(double k1, double k2, double u1){
	return (alpha(k1,k2,u1)+alpha(k2,k1,u1))/2.;
}

//beta(k1,k2)
double beta1(double k1, double k2, double u1){
	return u1*(k1*k1+k2*k2+2.*k1*k2*u1)/(2.*k1*k2);
}



// normalized Hubble parameter and it's scale factor derivative for DGP background - unused at the moment

/*
static double HAD(double a, double omega0, double omegarc){
	double omegaL= 1.-omega0;
	return  sqrt(omega0/pow(a,3)+omegaL+omegarc) - sqrt(omegarc);}

static double HAD1(double a, double omega0, double omegarc){
	double omegaL= 1.-omega0;
	return -3*omega0/(2.*(1+ sqrt(omegarc)/HA(a,omega0))*pow(a,3));}
*/



                                        ////
                                    /////////////
                            //////////////////////////////
                      /////////////////////////////////////////
                ///////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////


//beta function for nDGP
inline double beta(double a, double omega0, double omegarc){
	return 1.+HA(a,omega0)/sqrt(omegarc)*(1.+HA1(a,omega0)/(3.*HA(a,omega0)*HA(a,omega0)));}


//MODIFIED GRAVITY MODEL FUNCTIONS:
//p1,p2,p3 are theory parameters
//k1,k2,k3 are the mode magnitudes
//u1 angle between k2 and k3; u2 between k1 and k2; u3 between k1 and k3


// default is 3 model params. Add in parameters p4,p5,... as required.
// Note on theories with more than 3 parameters:
// 1) Add in extra parameters to param_type3 structure below
// 2) Add in extra parameters to params of funcn1 system of ODEs
// 3) Add in extra parameters to initn function and in params passed by function to solver (this is right after initial conditions are passed)

double mu(double a, double k0, double omega0, double p1, double p2, double p3 ){
	double h0 = 1./2997.92458;
//    	return 1.; // GR
	  return 1. + pow2(k0/a)/(3.*(pow2(k0/a)+pow3(omega0/pow3(a)-4.*(omega0-1.))/(2.*p1/pow2(h0)*pow2(4.-3.*omega0)))); //f(R) Hu- Sawicki
// 	return 1.+1./(3.*beta(a,omega0,p1)); //nDGP
}

double gamma2(double a, double omega0, double k0, double k1, double k2, double u1, double p1, double p2, double p3 ){
	double h0 = 1./2997.92458;
// 	 return 0. ; // GR
double var1 = pow3(a);
double var2 = pow3(var1);
double var3 = pow2(h0);
double var4 = pow2(3.*omega0-4.);
double var5 = pow3(omega0-4.*var1*(omega0-1.))/(2.*var2*p1/var3*var4);

	 return -(9.*pow2(k0/a)*pow2(omega0/var1)*pow(omega0-4.*var1*(-1.+omega0),5))/
			    (48.*pow(var1,5)*pow2(p1/var3)*pow2(HA(a,omega0))*pow2(var4)
			   *(pow2(k0/a)+var5)
			   *(pow2(k1/a)+var5)
	 		   *(pow2(k2/a)+var5)); //f(R) Hu- Sawicki

//  	return -1./(HA(a,omega0)*HA(a,omega0)*24.*pow(beta(a,omega0,p1),3)*p1)*pow(omega0/(a*a*a),2)*ker1(u1); //nDGP
}

// partially symmetric gamma3 (in last 2 arguments) for DGP
double gamma3(double a, double omega0, double k0, double k1, double k2, double k3, double u1,double u2, double u3, double p1, double p2, double p3){
	double h0 = 1./2997.92458;
// 	return 0. ; // GR

 double var1 = pow3(a);
 double var2 = pow3(var1);
 double var3 = pow2(h0);
 double var4 = pow2(3.*omega0-4.);
 double var5 = pow3(omega0-4.*var1*(omega0-1.))/(2.*var2*p1/var3*var4);
 double var6 = pow(var1,5);
 double var7 = var6*pow2(var1);
 double k23 = sqrt(k2*k2+k3*k3+2.*k2*k3*u1);

	 return 1./3.*(pow2(k0/(HA(a,omega0)*a))*pow3(omega0/var1)/(36.
											*(pow2(k0/a)+var5)
											*(pow2(k1/a)+var5)
											*(pow2(k2/a)+var5)
											*(pow2(k3/a)+var5)
					 *(pow2(k23/a)+var5))
			 *(-45./8.*(pow2(k23/a)+var5)/(var7*pow3(p1/var3))*pow(omega0-4.*var1*(omega0-1.),7)/pow3(var4)
			   +pow2(9./(4.*var6*pow2(p1/var3))*pow(omega0-4.*var1*(omega0-1.),5)/pow2(var4))));

//  return 1./(HA(a,omega0)*HA(a,omega0)*144.*pow(beta(a,omega0,p1),5)*pow(p1,2))*pow(omega0/(a*a*a),3)*ker1(u1)*ker1((k2*u2+k3*u3)/sqrt(k2*k2+2.*k2*k3*u1+k3*k3)); //nDGP
}

////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////
                  /////////////////////////////////////////////////
                         //////////////////////////////////
                                  //////////////////
                                        /////////





/*EDS Kernels */

double F2eds(double k1, double k2, double u1){
	return 5./7. * alphas(k1,k2,u1) + 2./7.*beta1(k1,k2,u1);
}

double G2eds(double k1, double k2, double u1){
	return 3./7. * alphas(k1,k2,u1) + 4./7.*beta1(k1,k2,u1);
}

double F3eds(double k1, double k2, double k3, double x1, double x2, double x3){
double k23 = sqrt(k2*k2+k3*k3+2.*k2*k3*x1);
double k13 = sqrt(k1*k1+k3*k3+2.*k1*k3*x3);
double k12 = sqrt(k1*k1+k2*k2+2.*k1*k2*x2);

  return 1./3.*(2./63.*2.*beta1(k1,k23,(k2*x2+k3*x3)/k23)*(2.*beta1(k2,k3,x1)+6./4.*alphas(k2,k3,x1))
                +1./18.*alpha(k1,k23,(k2*x2+k3*x3)/k23)*(2.*beta1(k2,k3,x1)+10./2.*alphas(k2,k3,x1))
                  +1./9.*alpha(k23,k1,(k2*x2+k3*x3)/k23)*(2.*beta1(k2,k3,x1)+6./4.*alphas(k2,k3,x1))

                 +2./63.*2.*beta1(k3,k12,(k1*x3+k2*x1)/k12)*(2.*beta1(k1,k2,x2)+6./4.*alphas(k1,k2,x2))
                 +1./18.*alpha(k3,k12,(k1*x3+k2*x1)/k12)*(2.*beta1(k1,k2,x2)+10./2.*alphas(k1,k2,x2))
                 +1./9.*alpha(k12,k3,(k1*x3+k2*x1)/k12)*(2.*beta1(k1,k2,x2)+6./4.*alphas(k1,k2,x2))

                 +2./63.*2.*beta1(k2,k13,(k1*x2+k3*x1)/k13)*(2.*beta1(k1,k3,x3)+6./4.*alphas(k1,k3,x3))
                 +1./18.*alpha(k2,k13,(k1*x2+k3*x1)/k13)*(2.*beta1(k1,k3,x3)+10./2.*alphas(k1,k3,x3))
                 +1./9.*alpha(k13,k2,(k1*x2+k3*x1)/k13)*(2.*beta1(k1,k3,x3)+6./4.*alphas(k1,k3,x3)));
  }


double G3eds(double k1, double k2, double k3, double x1, double x2, double x3){
  double k23 = sqrt(k2*k2+k3*k3+2.*k2*k3*x1);
  double k13 = sqrt(k1*k1+k3*k3+2.*k1*k3*x3);
  double k12 = sqrt(k1*k1+k2*k2+2.*k1*k2*x2);

  return 1./3.*(2./21.*2.*beta1(k1,k23,(k2*x2+k3*x3)/k23)*(2.*beta1(k2,k3,x1)+6./4.*alphas(k2,k3,x1))
               +1./42.*alpha(k1,k23,(k2*x2+k3*x3)/k23)*(2.*beta1(k2,k3,x1)+10./2.*alphas(k2,k3,x1))
               +1./21.*alpha(k23,k1,(k2*x2+k3*x3)/k23)*(2.*beta1(k2,k3,x1)+6./4.*alphas(k2,k3,x1))

               +2./21.*2.*beta1(k3,k12,(k1*x3+k2*x1)/k12)*(2.*beta1(k1,k2,x2)+6./4.*alphas(k1,k2,x2))
               +1./42.*alpha(k3,k12,(k1*x3+k2*x1)/k12)*(2.*beta1(k1,k2,x2)+10./2.*alphas(k1,k2,x2))
               +1./21.*alpha(k12,k3,(k1*x3+k2*x1)/k12)*(2.*beta1(k1,k2,x2)+6./4.*alphas(k1,k2,x2))

               +2./21.*2.*beta1(k2,k13,(k1*x2+k3*x1)/k13)*(2.*beta1(k1,k3,x3)+6./4.*alphas(k1,k3,x3))
               +1./42.*alpha(k2,k13,(k1*x2+k3*x1)/k13)*(2.*beta1(k1,k3,x3)+10./2.*alphas(k1,k3,x3))
               +1./21.*alpha(k13,k2,(k1*x2+k3*x1)/k13)*(2.*beta1(k1,k3,x3)+6./4.*alphas(k1,k3,x3)));

}



double F3edsb(double k1, double k2, double k3, double k4, double k5, double k6, double x1, double x2, double x3){
double k23 = k4;//sqrt(k2*k2+k3*k3+2.*k2*k3*x1);
double k12 = k5;//sqrt(k1*k1+k2*k2+2.*k1*k2*x2);
double k13 = k6;//sqrt(k1*k1+k3*k3+2.*k1*k3*x3);

  return 1./3.*(2./63.*2.*beta1(k1,k23,(k2*x2+k3*x3)/k23)*(2.*beta1(k2,k3,x1)+6./4.*alphas(k2,k3,x1))
                +1./18.*alpha(k1,k23,(k2*x2+k3*x3)/k23)*(2.*beta1(k2,k3,x1)+10./2.*alphas(k2,k3,x1))
                  +1./9.*alpha(k23,k1,(k2*x2+k3*x3)/k23)*(2.*beta1(k2,k3,x1)+6./4.*alphas(k2,k3,x1))

                 +2./63.*2.*beta1(k3,k12,(k1*x3+k2*x1)/k12)*(2.*beta1(k1,k2,x2)+6./4.*alphas(k1,k2,x2))
                 +1./18.*alpha(k3,k12,(k1*x3+k2*x1)/k12)*(2.*beta1(k1,k2,x2)+10./2.*alphas(k1,k2,x2))
                 +1./9.*alpha(k12,k3,(k1*x3+k2*x1)/k12)*(2.*beta1(k1,k2,x2)+6./4.*alphas(k1,k2,x2))

                 +2./63.*2.*beta1(k2,k13,(k1*x2+k3*x1)/k13)*(2.*beta1(k1,k3,x3)+6./4.*alphas(k1,k3,x3))
                 +1./18.*alpha(k2,k13,(k1*x2+k3*x1)/k13)*(2.*beta1(k1,k3,x3)+10./2.*alphas(k1,k3,x3))
                 +1./9.*alpha(k13,k2,(k1*x2+k3*x1)/k13)*(2.*beta1(k1,k3,x3)+6./4.*alphas(k1,k3,x3)));
  }


double G3edsb(double k1, double k2, double k3, double k4, double k5, double k6, double x1, double x2, double x3){
  double k23 = k4;//sqrt(k2*k2+k3*k3+2.*k2*k3*x1);
  double k13 = k6;//sqrt(k1*k1+k3*k3+2.*k1*k3*x3);
  double k12 = k5;//sqrt(k1*k1+k2*k2+2.*k1*k2*x2);

  return 1./3.*(2./21.*2.*beta1(k1,k23,(k2*x2+k3*x3)/k23)*(2.*beta1(k2,k3,x1)+6./4.*alphas(k2,k3,x1))
               +1./42.*alpha(k1,k23,(k2*x2+k3*x3)/k23)*(2.*beta1(k2,k3,x1)+10./2.*alphas(k2,k3,x1))
               +1./21.*alpha(k23,k1,(k2*x2+k3*x3)/k23)*(2.*beta1(k2,k3,x1)+6./4.*alphas(k2,k3,x1))

               +2./21.*2.*beta1(k3,k12,(k1*x3+k2*x1)/k12)*(2.*beta1(k1,k2,x2)+6./4.*alphas(k1,k2,x2))
               +1./42.*alpha(k3,k12,(k1*x3+k2*x1)/k12)*(2.*beta1(k1,k2,x2)+10./2.*alphas(k1,k2,x2))
               +1./21.*alpha(k12,k3,(k1*x3+k2*x1)/k12)*(2.*beta1(k1,k2,x2)+6./4.*alphas(k1,k2,x2))

               +2./21.*2.*beta1(k2,k13,(k1*x2+k3*x1)/k13)*(2.*beta1(k1,k3,x3)+6./4.*alphas(k1,k3,x3))
               +1./42.*alpha(k2,k13,(k1*x2+k3*x1)/k13)*(2.*beta1(k1,k3,x3)+10./2.*alphas(k1,k3,x3))
               +1./21.*alpha(k13,k2,(k1*x2+k3*x1)/k13)*(2.*beta1(k1,k3,x3)+6./4.*alphas(k1,k3,x3)));

}



double F4eds(double k[], double x[]){
  // mags^2
  double k0s = pow2(k[0]);
  double k1s = pow2(k[1]);
  double k2s = pow2(k[2]);
  double k3s = pow2(k[3]);
  double k1234 = k0s+k1s+k2s+k3s
              + 2.*k[0]*k[1]*x[1] + 2.*k[0]*k[2]*x[2] + 2.*k[1]*k[2]*x[3]
              + 2.*k[0]*k[3]*x[5] + 2.*k[1]*k[3]*x[0] + 2.*k[2]*k[3]*x[4];
  double k123 = k0s+k1s+k2s + 2.*k[0]*k[1]*x[1] + 2.*k[0]*k[2]*x[2] + 2.*k[1]*k[2]*x[3];
  double k234 = k1s+k2s+k3s + 2.*k[1]*k[2]*x[3] + 2.*k[1]*k[3]*x[0] + 2.*k[2]*k[3]*x[4];
  double k341 = k2s+k3s+k0s + 2.*k[2]*k[3]*x[4] + 2.*k[2]*k[0]*x[2] + 2.*k[3]*k[0]*x[5];
  double k412 = k3s+k0s+k1s + 2.*k[3]*k[0]*x[5] + 2.*k[3]*k[1]*x[0] + 2.*k[0]*k[1]*x[1];
  double k12 = k0s+k1s+2.*k[0]*k[1]*x[1];
  double k23 = k1s+k2s+2.*k[1]*k[2]*x[3];
  double k34 = k2s+k3s+2.*k[2]*k[3]*x[4];
  double k41 = k3s+k0s+2.*k[3]*k[0]*x[5];
  double k13 = k0s+k2s+2.*k[0]*k[2]*x[2];
  double k24 = k1s+k3s+2.*k[1]*k[3]*x[0];

  double k12rt = sqrt(k12);
  double k23rt = sqrt(k23);
  double k34rt = sqrt(k34);
  double k41rt = sqrt(k41);
  double k13rt = sqrt(k13);
  double k24rt = sqrt(k24);
  // x index index
  // 01 = 1 | 02 = 2 | 03 = 5 | 12 = 3 | 13 = 0 | 23 = 4


  return 1./396.*(
        27.*(k0s+k[0]*k[1]*x[1]+k[0]*k[2]*x[2]+k[0]*k[3]*x[5])/k0s*F3edsb(k[1],k[2],k[3],k34rt,k23rt,k24rt,x[4],x[3],x[0])
      +(27.*(k234 + k[0]*(k[1]*x[1] + k[2]*x[2] + k[3]*x[5]))/k234 + 6.*k1234*k[0]*(k[1]*x[1]+k[2]*x[2]+k[3]*x[5])/k0s/k234)*G3edsb(k[1],k[2],k[3],k34rt,k23rt,k24rt,x[4],x[3],x[0])

    +27.*(k1s+k[1]*k[2]*x[3]+k[1]*k[3]*x[0]+k[1]*k[0]*x[1])/k1s*F3edsb(k[2],k[3],k[0],k41rt,k34rt,k13rt,x[5],x[4],x[2])
    +(27.*(k341 + k[1]*(k[2]*x[3] + k[3]*x[0] + k[0]*x[1]))/k341 + 6.*k1234*k[1]*(k[2]*x[3]+k[3]*x[0]+k[0]*x[1])/k1s/k341)*G3edsb(k[2],k[3],k[0],k41rt,k34rt,k13rt,x[5],x[4],x[2])

    +27.*(k2s+k[2]*k[3]*x[4]+k[2]*k[0]*x[2]+k[2]*k[1]*x[3])/k2s*F3edsb(k[3],k[0],k[1],k12rt,k41rt,k24rt,x[1],x[5],x[0])
    +(27.*(k412 + k[2]*(k[3]*x[4] + k[0]*x[2] + k[1]*x[3]))/k412 + 6.*k1234*k[2]*(k[3]*x[4]+k[0]*x[2]+k[1]*x[3])/k2s/k412)*G3edsb(k[3],k[0],k[1],k12rt,k41rt,k24rt,x[1],x[5],x[0])

    +27.*(k3s+k[3]*k[0]*x[5]+k[3]*k[1]*x[0]+k[3]*k[2]*x[4])/k3s*F3edsb(k[0],k[1],k[2],k23rt,k12rt,k13rt,x[3],x[1],x[2])
    +(27.*(k123 + k[3]*(k[0]*x[5] + k[1]*x[0] + k[2]*x[4]))/k123 + 6.*k1234*k[3]*(k[0]*x[5]+k[1]*x[0]+k[2]*x[4])/k3s/k123)*G3edsb(k[0],k[1],k[2],k23rt,k12rt,k13rt,x[3],x[1],x[2])

    +18.*(k12+k[2]*k[0]*x[2]+k[2]*k[1]*x[3]+k[3]*k[0]*x[5]+k[3]*k[1]*x[0])/k12*G2eds(k[0],k[1],x[1])*F2eds(k[2],k[3],x[4])
    +18.*(k23+k[3]*k[1]*x[0]+k[3]*k[2]*x[4]+k[0]*k[1]*x[1]+k[0]*k[2]*x[2])/k23*G2eds(k[1],k[2],x[3])*F2eds(k[3],k[0],x[5])
    +18.*(k34+k[0]*k[2]*x[2]+k[0]*k[3]*x[5]+k[1]*k[2]*x[3]+k[1]*k[3]*x[0])/k34*G2eds(k[2],k[3],x[4])*F2eds(k[0],k[1],x[1])
    +18.*(k13+k[0]*k[1]*x[1]+k[0]*k[3]*x[5]+k[1]*k[2]*x[3]+k[2]*k[3]*x[4])/k13*G2eds(k[0],k[2],x[2])*F2eds(k[1],k[3],x[0])
    +18.*(k24+k[1]*k[2]*x[3]+k[1]*k[0]*x[1]+k[3]*k[0]*x[5]+k[2]*k[3]*x[4])/k24*G2eds(k[1],k[3],x[0])*F2eds(k[0],k[2],x[2])
    +18.*(k41+k[0]*k[1]*x[1]+k[0]*k[2]*x[2]+k[3]*k[1]*x[0]+k[2]*k[3]*x[4])/k41*G2eds(k[0],k[3],x[5])*F2eds(k[1],k[2],x[3])

    +4.*k1234*(k[0]*k[2]*x[2]+k[0]*k[3]*x[5]+k[1]*k[2]*x[3]+k[1]*k[3]*x[0])/k12/k34*G2eds(k[0],k[1],x[1])*G2eds(k[2],k[3],x[4])
    +4.*k1234*(k[1]*k[3]*x[0]+k[1]*k[0]*x[1]+k[2]*k[3]*x[4]+k[2]*k[0]*x[2])/k23/k41*G2eds(k[1],k[2],x[3])*G2eds(k[3],k[0],x[5])
    +4.*k1234*(k[1]*k[0]*x[1]+k[1]*k[2]*x[3]+k[3]*k[0]*x[5]+k[3]*k[2]*x[4])/k24/k13*G2eds(k[1],k[3],x[0])*G2eds(k[0],k[2],x[2]));

}

double G4eds(double k[], double x[]){
  // mags^2
  double k0s = pow2(k[0]);
  double k1s = pow2(k[1]);
  double k2s = pow2(k[2]);
  double k3s = pow2(k[3]);

  double k1234 = k0s+k1s+k2s + k3s
              + 2.*k[0]*k[1]*x[1] + 2.*k[0]*k[2]*x[2] + 2.*k[1]*k[2]*x[3]
              + 2.*k[0]*k[3]*x[5] + 2.*k[1]*k[3]*x[0] + 2.*k[2]*k[3]*x[4];
  double k123 = k0s+k1s+k2s + 2.*k[0]*k[1]*x[1] + 2.*k[0]*k[2]*x[2] + 2.*k[1]*k[2]*x[3];
  double k234 = k1s+k2s+k3s + 2.*k[1]*k[2]*x[3] + 2.*k[1]*k[3]*x[0] + 2.*k[2]*k[3]*x[4];
  double k341 = k2s+k3s+k0s + 2.*k[2]*k[3]*x[4] + 2.*k[2]*k[0]*x[2] + 2.*k[3]*k[0]*x[5];
  double k412 = k3s+k0s+k1s + 2.*k[3]*k[0]*x[5] + 2.*k[3]*k[1]*x[0] + 2.*k[0]*k[1]*x[1];
  double k12 = k0s+k1s+2.*k[0]*k[1]*x[1];
  double k23 = k1s+k2s+2.*k[1]*k[2]*x[3];
  double k34 = k2s+k3s+2.*k[2]*k[3]*x[4];
  double k41 = k3s+k0s+2.*k[3]*k[0]*x[5];
  double k13 = k0s+k2s+2.*k[0]*k[2]*x[2];
  double k24 = k1s+k3s+2.*k[1]*k[3]*x[0];


  double k12rt = sqrt(k12);
  double k23rt = sqrt(k23);
  double k34rt = sqrt(k34);
  double k41rt = sqrt(k41);
  double k13rt = sqrt(k13);
  double k24rt = sqrt(k24);

  // x index index
  // 01 = 1 | 02 = 2 | 03 = 5 | 12 = 3 | 13 = 0 | 23 = 4

  return 1./396.*(
        9.*(k0s+k[0]*k[1]*x[1]+k[0]*k[2]*x[2]+k[0]*k[3]*x[5])/k0s*F3edsb(k[1],k[2],k[3],k34rt,k23rt,k24rt,x[4],x[3],x[0])
      +(9.*(k234 + k[0]*(k[1]*x[1] + k[2]*x[2] + k[3]*x[5]))/k234 + 24.*k1234*k[0]*(k[1]*x[1]+k[2]*x[2]+k[3]*x[5])/k0s/k234)*G3edsb(k[1],k[2],k[3],k34rt,k23rt,k24rt,x[4],x[3],x[0])

    +9.*(k1s+k[1]*k[2]*x[3]+k[1]*k[3]*x[0]+k[1]*k[0]*x[1])/k1s*F3edsb(k[2],k[3],k[0],k41rt,k34rt,k13rt,x[5],x[4],x[2])
    +(9.*(k341 + k[1]*(k[2]*x[3] + k[3]*x[0] + k[0]*x[1]))/k341 + 24.*k1234*k[1]*(k[2]*x[3]+k[3]*x[0]+k[0]*x[1])/k1s/k341)*G3edsb(k[2],k[3],k[0],k41rt,k34rt,k13rt,x[5],x[4],x[2])

    +9.*(k2s+k[2]*k[3]*x[4]+k[2]*k[0]*x[2]+k[2]*k[1]*x[3])/k2s*F3edsb(k[3],k[0],k[1],k12rt,k41rt,k24rt,x[1],x[5],x[0])
    +(9.*(k412 + k[2]*(k[3]*x[4] + k[0]*x[2] + k[1]*x[3]))/k412 + 24.*k1234*k[2]*(k[3]*x[4]+k[0]*x[2]+k[1]*x[3])/k2s/k412)*G3edsb(k[3],k[0],k[1],k12rt,k41rt,k24rt,x[1],x[5],x[0])

    +9.*(k3s+k[3]*k[0]*x[5]+k[3]*k[1]*x[0]+k[3]*k[2]*x[4])/k3s*F3edsb(k[0],k[1],k[2],k23rt,k12rt,k13rt,x[3],x[1],x[2])
    +(9.*(k123 + k[3]*(k[0]*x[5] + k[1]*x[0] + k[2]*x[4]))/k123 + 24.*k1234*k[3]*(k[0]*x[5]+k[1]*x[0]+k[2]*x[4])/k3s/k123)*G3edsb(k[0],k[1],k23rt,k12rt,k13rt,k[2],x[3],x[1],x[2])

    +6.*(k12+k[2]*k[0]*x[2]+k[2]*k[1]*x[3]+k[3]*k[0]*x[5]+k[3]*k[1]*x[0])/k12*G2eds(k[0],k[1],x[1])*F2eds(k[2],k[3],x[4])
    +6.*(k23+k[3]*k[1]*x[0]+k[3]*k[2]*x[4]+k[0]*k[1]*x[1]+k[0]*k[2]*x[2])/k23*G2eds(k[1],k[2],x[3])*F2eds(k[3],k[0],x[5])
    +6.*(k34+k[0]*k[2]*x[2]+k[0]*k[3]*x[5]+k[1]*k[2]*x[3]+k[1]*k[3]*x[0])/k34*G2eds(k[2],k[3],x[4])*F2eds(k[0],k[1],x[1])
    +6.*(k13+k[0]*k[1]*x[1]+k[0]*k[3]*x[5]+k[1]*k[2]*x[3]+k[2]*k[3]*x[4])/k13*G2eds(k[0],k[2],x[2])*F2eds(k[1],k[3],x[0])
    +6.*(k24+k[1]*k[2]*x[3]+k[1]*k[0]*x[1]+k[3]*k[0]*x[5]+k[2]*k[3]*x[4])/k24*G2eds(k[1],k[3],x[0])*F2eds(k[0],k[2],x[2])
    +6.*(k41+k[0]*k[1]*x[1]+k[0]*k[2]*x[2]+k[3]*k[1]*x[0]+k[2]*k[3]*x[4])/k41*G2eds(k[0],k[3],x[5])*F2eds(k[1],k[2],x[3])

    +16.*k1234*(k[0]*k[2]*x[2]+k[0]*k[3]*x[5]+k[1]*k[2]*x[3]+k[1]*k[3]*x[0])/k12/k34*G2eds(k[0],k[1],x[1])*G2eds(k[2],k[3],x[4])
    +16.*k1234*(k[1]*k[3]*x[0]+k[1]*k[0]*x[1]+k[2]*k[3]*x[4]+k[2]*k[0]*x[2])/k23/k41*G2eds(k[1],k[2],x[3])*G2eds(k[3],k[0],x[5])
    +16.*k1234*(k[1]*k[0]*x[1]+k[1]*k[2]*x[3]+k[3]*k[0]*x[5]+k[3]*k[2]*x[4])/k24/k13*G2eds(k[1],k[3],x[0])*G2eds(k[0],k[2],x[2]));
}



// x index index
// 12 = 3 | 13 = 4 | 23 = 1 | 14 = 5 | 24 = 2 | 34 = -1.
 double F4edsb(double k1, double k2, double k3, double k4, double x1, double x2,double x3, double x4, double x5){
// mags^2
double x6 = XMIN;
double k1s = pow2(k1);
double k2s = pow2(k2);
double k3s = pow2(k3);
double k4s = pow2(k4);
double k1234 = k1s+k2s+k3s + k4s
           + 2.*k1*k2*x3 + 2.*k1*k3*x4 + 2.*k2*k3*x1
           + 2.*k1*k4*x5 + 2.*k2*k4*x2 + 2.*k3*k4*x6;
double k123 = k1s+k2s+k3s + 2.*k1*k2*x3 + 2.*k1*k3*x4 + 2.*k2*k3*x1;
double k234 = k2s+k3s+k4s + 2.*k2*k3*x1 + 2.*k2*k4*x2 + 2.*k3*k4*x6;
double k341 = k3s+k4s+k1s + 2.*k3*k4*x6 + 2.*k3*k1*x4 + 2.*k4*k1*x5;
double k412 = k4s+k1s+k2s + 2.*k4*k1*x5 + 2.*k4*k2*x2 + 2.*k1*k2*x3;
double k12 = k1s+k2s+2.*k1*k2*x3;
double k23 = k2s+k3s+2.*k2*k3*x1;
double k34 = k3s+k4s+2.*k3*k4*x6;
double k41 = k4s+k1s+2.*k4*k1*x5;
double k13 = k1s+k3s+2.*k1*k3*x4;
double k24 = k2s+k4s+2.*k2*k4*x2;


double k12rt = sqrt(k12);
double k23rt = sqrt(k23);
double k34rt = sqrt(k34);
double k41rt = sqrt(k41);
double k13rt = sqrt(k13);
double k24rt = sqrt(k24);

// x index index
// 01 = 1 | 02 = 2 | 03 = 5 | 12 = 3 | 13 = 0 | 23 = 4

return 1./396.*(
     27.*(k1s+k1*k2*x3+k1*k3*x4+k1*k4*x5)/k1s*F3edsb(k2,k3,k4,k34rt,k23rt,k24rt,x6,x1,x2)
   +(27.*(k234 + k1*(k2*x3 + k3*x4 + k4*x5))/k234 + 6.*k1234*k1*(k2*x3+k3*x4+k4*x5)/k1s/k234)*G3edsb(k2,k3,k4,k34rt,k23rt,k24rt,x6,x1,x2)

 +27.*(k2s+k2*k3*x1+k2*k4*x2+k2*k1*x3)/k2s*F3edsb(k3,k4,k1,k41rt,k34rt,k13rt,x5,x6,x4)
 +(27.*(k341 + k2*(k3*x1 + k4*x2 + k1*x3))/k341 + 6.*k1234*k2*(k3*x1+k4*x2+k1*x3)/k2s/k341)*G3edsb(k3,k4,k1,k41rt,k34rt,k13rt,x5,x6,x4)

 +27.*(k3s+k3*k4*x6+k3*k1*x4+k3*k2*x1)/k3s*F3edsb(k4,k1,k2,k12rt,k41rt,k24rt,x3,x5,x2)
 +(27.*(k412 + k3*(k4*x6 + k1*x4 + k2*x1))/k412 + 6.*k1234*k3*(k4*x6+k1*x4+k2*x1)/k3s/k412)*G3edsb(k4,k1,k2,k12rt,k41rt,k24rt,x3,x5,x2)

 +27.*(k4s+k4*k1*x5+k4*k2*x2+k4*k3*x6)/k4s*F3edsb(k1,k2,k3,k23rt,k12rt,k13rt,x1,x3,x4)
 +(27.*(k123 + k4*(k1*x5 + k2*x2 + k3*x6))/k123 + 6.*k1234*k4*(k1*x5+k2*x2+k3*x6)/k4s/k123)*G3edsb(k1,k2,k3,k23rt,k12rt,k13rt,x1,x3,x4)

 +18.*(k12+k3*k1*x4+k3*k2*x1+k4*k1*x5+k4*k2*x2)/k12*G2eds(k1,k2,x3)*F2eds(k3,k4,x6)
 +18.*(k23+k4*k2*x2+k4*k3*x6+k1*k2*x3+k1*k3*x4)/k23*G2eds(k2,k3,x1)*F2eds(k4,k1,x5)
 +18.*(k34+k1*k3*x4+k1*k4*x5+k2*k3*x1+k2*k4*x2)/k34*G2eds(k3,k4,x6)*F2eds(k1,k2,x3)
 +18.*(k13+k1*k2*x3+k1*k4*x5+k2*k3*x1+k3*k4*x6)/k13*G2eds(k1,k3,x4)*F2eds(k2,k4,x2)
 +18.*(k24+k2*k3*x1+k2*k1*x3+k4*k1*x5+k3*k4*x6)/k24*G2eds(k2,k4,x2)*F2eds(k1,k3,x4)
 +18.*(k41+k1*k2*x3+k1*k3*x4+k4*k2*x2+k3*k4*x6)/k41*G2eds(k1,k4,x5)*F2eds(k2,k3,x1)

 +4.*k1234*(k1*k3*x4+k1*k4*x5+k2*k3*x1+k2*k4*x2)/k12/k34*G2eds(k1,k2,x3)*G2eds(k3,k4,x6)
 +4.*k1234*(k2*k4*x2+k2*k1*x3+k3*k4*x6+k3*k1*x4)/k23/k41*G2eds(k2,k3,x1)*G2eds(k4,k1,x5)
 +4.*k1234*(k2*k1*x3+k2*k3*x1+k4*k1*x5+k4*k3*x6)/k24/k13*G2eds(k2,k4,x2)*G2eds(k1,k3,x4));

}



double G4edsb(double k1, double k2, double k3, double k4, double x1, double x2,double x3, double x4, double x5){
  // mags^2
  double x6 = XMIN;
  double k1s = pow2(k1);
  double k2s = pow2(k2);
  double k3s = pow2(k3);
  double k4s = pow2(k4);
  double k1234 = k1s+k2s+k3s + k4s
             + 2.*k1*k2*x3 + 2.*k1*k3*x4 + 2.*k2*k3*x1
             + 2.*k1*k4*x5 + 2.*k2*k4*x2 + 2.*k3*k4*x6;
  double k123 = k1s+k2s+k3s + 2.*k1*k2*x3 + 2.*k1*k3*x4 + 2.*k2*k3*x1;
  double k234 = k2s+k3s+k4s + 2.*k2*k3*x1 + 2.*k2*k4*x2 + 2.*k3*k4*x6;
  double k341 = k3s+k4s+k1s + 2.*k3*k4*x6 + 2.*k3*k1*x4 + 2.*k4*k1*x5;
  double k412 = k4s+k1s+k2s + 2.*k4*k1*x5 + 2.*k4*k2*x2 + 2.*k1*k2*x3;
  double k12 = k1s+k2s+2.*k1*k2*x3;
  double k23 = k2s+k3s+2.*k2*k3*x1;
  double k34 = k3s+k4s+2.*k3*k4*x6;
  double k41 = k4s+k1s+2.*k4*k1*x5;
  double k13 = k1s+k3s+2.*k1*k3*x4;
  double k24 = k2s+k4s+2.*k2*k4*x2;


  double k12rt = sqrt(k12);
  double k23rt = sqrt(k23);
  double k34rt = sqrt(k34);
  double k41rt = sqrt(k41);
  double k13rt = sqrt(k13);
  double k24rt = sqrt(k24);

  // x index index
  // 01 = 1 | 02 = 2 | 03 = 5 | 12 = 3 | 13 = 0 | 23 = 4
    return 1./396.*(
         9.*(k1s+k1*k2*x3+k1*k3*x4+k1*k4*x5)/k1s*F3edsb(k2,k3,k4,k34rt,k23rt,k24rt,x6,x1,x2)
       +(9.*(k234 + k1*(k2*x3 + k3*x4 + k4*x5))/k234 + 24.*k1234*k1*(k2*x3+k3*x4+k4*x5)/k1s/k234)*G3edsb(k2,k3,k4,k34rt,k23rt,k24rt,x6,x1,x2)

     +9.*(k2s+k2*k3*x1+k2*k4*x2+k2*k1*x3)/k2s*F3edsb(k3,k4,k1,k41rt,k34rt,k13rt,x5,x6,x4)
     +(9.*(k341 + k2*(k3*x1 + k4*x2 + k1*x3))/k341 + 24.*k1234*k2*(k3*x1+k4*x2+k1*x3)/k2s/k341)*G3edsb(k3,k4,k1,k41rt,k34rt,k13rt,x5,x6,x4)

     +9.*(k3s+k3*k4*x6+k3*k1*x4+k3*k2*x1)/k3s*F3edsb(k4,k1,k2,k12rt,k41rt,k24rt,x3,x5,x2)
     +(9.*(k412 + k3*(k4*x6 + k1*x4 + k2*x1))/k412 + 24.*k1234*k3*(k4*x6+k1*x4+k2*x1)/k3s/k412)*G3edsb(k4,k1,k2,k12rt,k41rt,k24rt,x3,x5,x2)

     +9.*(k4s+k4*k1*x5+k4*k2*x2+k4*k3*x6)/k4s*F3edsb(k1,k2,k3,k23rt,k12rt,k13rt,x1,x3,x4)
     +(9.*(k123 + k4*(k1*x5 + k2*x2 + k3*x6))/k123 + 24.*k1234*k4*(k1*x5+k2*x2+k3*x6)/k4s/k123)*G3edsb(k1,k2,k3,k23rt,k12rt,k13rt,x1,x3,x4)

     +6.*(k12+k3*k1*x4+k3*k2*x1+k4*k1*x5+k4*k2*x2)/k12*G2eds(k1,k2,x3)*F2eds(k3,k4,x6)
     +6.*(k23+k4*k2*x2+k4*k3*x6+k1*k2*x3+k1*k3*x4)/k23*G2eds(k2,k3,x1)*F2eds(k4,k1,x5)
     +6.*(k34+k1*k3*x4+k1*k4*x5+k2*k3*x1+k2*k4*x2)/k34*G2eds(k3,k4,x6)*F2eds(k1,k2,x3)
     +6.*(k13+k1*k2*x3+k1*k4*x5+k2*k3*x1+k3*k4*x6)/k13*G2eds(k1,k3,x4)*F2eds(k2,k4,x2)
     +6.*(k24+k2*k3*x1+k2*k1*x3+k4*k1*x5+k3*k4*x6)/k24*G2eds(k2,k4,x2)*F2eds(k1,k3,x4)
     +6.*(k41+k1*k2*x3+k1*k3*x4+k4*k2*x2+k3*k4*x6)/k41*G2eds(k1,k4,x5)*F2eds(k2,k3,x1)

     +16.*k1234*(k1*k3*x4+k1*k4*x5+k2*k3*x1+k2*k4*x2)/k12/k34*G2eds(k1,k2,x3)*G2eds(k3,k4,x6)
     +16.*k1234*(k2*k4*x2+k2*k1*x3+k3*k4*x6+k3*k1*x4)/k23/k41*G2eds(k2,k3,x1)*G2eds(k4,k1,x5)
     +16.*k1234*(k2*k1*x3+k2*k3*x1+k4*k1*x5+k4*k3*x6)/k24/k13*G2eds(k2,k4,x2)*G2eds(k1,k3,x4));
}



// nDGP analytic kernels  (See 0902.0618v3 for expressions)
double C3SYM(double k1, double k2, double k3, double x1, double x2, double x3){

  return 2./3.*(beta1(k1,sqrt(k2*k2+k3*k3+2*k2*k3*x1),(k2*x2+k3*x3)/sqrt(k2*k2+k3*k3+2.*k2*k3*x1))*ker1(x1)
               +beta1(k3,sqrt(k1*k1+k2*k2+2*k1*k2*x2),(k1*x3+k2*x1)/sqrt(k1*k1+k2*k2+2.*k1*k2*x2))*ker1(x2)
               +beta1(k2,sqrt(k1*k1+k3*k3+2*k1*k3*x3),(k1*x2+k3*x1)/sqrt(k1*k1+k3*k3+2.*k1*k3*x3))*ker1(x3));
  }

double C3SYMb(double k[],double x[]){
//  k[3]=sqrt(k2*k2+k3*k3+2*k2*k3*x1)
//  k[4]=sqrt(k1*k1+k2*k2+2*k1*k2*x2)
//  k[5]=sqrt(k1*k1+k3*k3+2*k1*k3*x3)
// x[0] = k2.k3
// x[1] = k1.k2
// x[2] = k1.k3
return 2./3.*(beta1(k[0],k[3],(k[1]*x[1]+k[2]*x[2])/k[3])*ker1(x[0])
             +beta1(k[2],k[4],(k[0]*x[2]+k[1]*x[0])/k[4])*ker1(x[1])
             +beta1(k[1],k[5],(k[0]*x[1]+k[2]*x[0])/k[5])*ker1(x[2]));
}

double F3SYM(double k1, double k2, double k3, double x1, double x2, double x3){
  return 1./3.*(alpha(k1,sqrt(k2*k2+k3*k3+2.*k2*k3*x1),(k2*x2+k3*x3)/sqrt(k2*k2+k3*k3+2.*k2*k3*x1))*ker1(x1)
               +alpha(k3,sqrt(k1*k1+k2*k2+2.*k1*k2*x2),(k1*x3+k2*x1)/sqrt(k1*k1+k2*k2+2.*k1*k2*x2))*ker1(x2)
               +alpha(k2,sqrt(k1*k1+k3*k3+2.*k1*k3*x3),(k1*x2+k3*x1)/sqrt(k1*k1+k3*k3+2.*k1*k3*x3))*ker1(x3));
    }

double F3SYMb(double k[],double x[]){
return 1./3.*(alpha(k[0],k[3],(k[1]*x[1]+k[2]*x[2])/k[3])*ker1(x[0])
             +alpha(k[2],k[4],(k[0]*x[2]+k[1]*x[0])/k[4])*ker1(x[1])
             +alpha(k[1],k[5],(k[0]*x[1]+k[2]*x[0])/k[5])*ker1(x[2]));
}


////k2,p,p,XMIN,u,-u

double I3SYM(double k1, double k2, double k3, double x1, double x2, double x3){
  return 1./3.*(alpha(sqrt(k2*k2+k3*k3+2.*k2*k3*x1),k1,(k2*x2+k3*x3)/sqrt(k2*k2+k3*k3+2.*k2*k3*x1))*ker1(x1)
               +alpha(sqrt(k1*k1+k2*k2+2.*k1*k2*x2),k3,(k1*x3+k2*x1)/sqrt(k1*k1+k2*k2+2.*k1*k2*x2))*ker1(x2)
               +alpha(sqrt(k1*k1+k3*k3+2.*k1*k3*x3),k2,(k1*x2+k3*x1)/sqrt(k1*k1+k3*k3+2.*k1*k3*x3))*ker1(x3));
    }

double I3SYMb(double k[],double x[]){
return 1./3.*(alpha(k[3],k[0],(k[1]*x[1]+k[2]*x[2])/k[3])*ker1(x[0])
             +alpha(k[4],k[2],(k[0]*x[2]+k[1]*x[0])/k[4])*ker1(x[1])
             +alpha(k[5],k[1],(k[0]*x[1]+k[2]*x[0])/k[5])*ker1(x[2]));
}

double J3SYM(double k1, double k2, double k3, double x1, double x2, double x3){
  return 2./3.*(beta1(k2,k3,x1)*ker1((k2*x2+k3*x3)/sqrt(k2*k2+k3*k3+2.*k2*k3*x1))
               +beta1(k3,k1,x3)*ker1((k1*x2+k3*x1)/sqrt(k1*k1+k3*k3+2.*k1*k3*x3))
               +beta1(k1,k2,x2)*ker1((k1*x3+k2*x1)/sqrt(k1*k1+k2*k2+2.*k1*k2*x2)));
  }

double J3SYMb(double k[],double x[]){
  return 2./3.*(beta1(k[1],k[2],x[0])*ker1((k[1]*x[1]+k[2]*x[2])/k[3])
               +beta1(k[2],k[0],x[2])*ker1((k[0]*x[1]+k[2]*x[0])/k[5])
               +beta1(k[0],k[1],x[1])*ker1((k[0]*x[2]+k[1]*x[0])/k[4]));
  }

double K3SYM(double k1, double k2, double k3, double x1, double x2, double x3){
    return 1./3.*(alphas(k2,k3,x1)*ker1((k2*x2+k3*x3)/sqrt(k2*k2+k3*k3+2.*k2*k3*x1))
                 +alphas(k3,k1,x3)*ker1((k1*x2+k3*x1)/sqrt(k1*k1+k3*k3+2.*k1*k3*x3))
                 +alphas(k1,k2,x2)*ker1((k1*x3+k2*x1)/sqrt(k1*k1+k2*k2+2.*k1*k2*x2)));
    }


double K3SYMb(double k[],double x[]){
  return 1./3.*(alphas(k[1],k[2],x[0])*ker1((k[1]*x[1]+k[2]*x[2])/k[3])
               +alphas(k[2],k[0],x[2])*ker1((k[0]*x[1]+k[2]*x[0])/k[5])
               +alphas(k[0],k[1],x[1])*ker1((k[0]*x[2]+k[1]*x[0])/k[4]));
  }


double L3SYM(double k1, double k2, double k3, double x1, double x2, double x3){
    return 1./3.*(ker1(x1)*ker1((k2*x2+k3*x3)/sqrt(k2*k2+k3*k3+2.*k2*k3*x1))
                 +ker1(x3)*ker1((k1*x2+k3*x1)/sqrt(k1*k1+k3*k3+2.*k1*k3*x3))
                 +ker1(x2)*ker1((k1*x3+k2*x1)/sqrt(k1*k1+k2*k2+2.*k1*k2*x2)));
        }

double L3SYMb(double k[],double x[]){
  return 1./3.*(ker1(x[0])*ker1((k[1]*x[1]+k[2]*x[2])/k[3])
               +ker1(x[2])*ker1((k[0]*x[1]+k[2]*x[0])/k[5])
               +ker1(x[1])*ker1((k[0]*x[2]+k[1]*x[0])/k[4]));
  }




double F2ndgp(double k1, double k2, double u1){
 return pow2(1./dnorm_spt)*(pow2(D_spt)*F2eds(k1,k2,u1) + F_spt*(1.-pow2(u1)));
}

double G2ndgp(double k1, double k2, double u1){
 return pow2(1./dnorm_spt)*(-fdgp_spt*pow2(D_spt)*G2eds(k1,k2,u1) - Fd_spt/H_spt*(1-pow2(u1)));
}


inline double F3edsd(double k[], double x[]){
    return F3edsb(k[0],k[1],k[2],k[3],k[4],k[5],x[0],x[1],x[2]);
  }


inline double G3edsd(double k[], double x[]){
  return G3edsb(k[0],k[1],k[2],k[3],k[4],k[5],x[0],x[1],x[2]);
}



double F3ndgp(double k1, double k2, double k3, double x1, double x2, double x3){
  return pow3(1./dnorm_spt)*(pow3(D_spt)*F3eds(k1,k2,k3,x1,x2,x3)+Cx_spt*C3SYM(k1,k2,k3,x1,x2,x3)
          +F3_spt*F3SYM(k1,k2,k3,x1,x2,x3)+I_spt*I3SYM(k1,k2,k3,x1,x2,x3) + J_spt*J3SYM(k1,k2,k3,x1,x2,x3)
          +K_spt*K3SYM(k1,k2,k3,x1,x2,x3)+(KL_spt-K_spt)*L3SYM(k1,k2,k3,x1,x2,x3));
}

double G3ndgp(double k1, double k2, double k3, double x1, double x2, double x3){
  return pow3(1./dnorm_spt)*(-fdgp_spt*pow3(D_spt)*G3eds(k1,k2,k3,x1,x2,x3)-Cdx_spt/H_spt*C3SYM(k1,k2,k3,x1,x2,x3)
          -(F3d_spt/H_spt-D_spt*fdgp_spt*F_spt)*F3SYM(k1,k2,k3,x1,x2,x3)-(Id_spt/H_spt-D_spt*Fd_spt/H_spt)*I3SYM(k1,k2,k3,x1,x2,x3)
          - Jd_spt/H_spt*J3SYM(k1,k2,k3,x1,x2,x3) - Kd_spt/H_spt*K3SYM(k1,k2,k3,x1,x2,x3) - (KLd_spt-Kd_spt)/H_spt*L3SYM(k1,k2,k3,x1,x2,x3));
}


double F3bdgp(double k[],double x[]){
  double norm=pow3(1./dnorm_spt);
  return norm*(pow3(D_spt)*F3edsd(k,x) + Cx_spt*C3SYMb(k,x) + F3_spt*F3SYMb(k,x) + I_spt*I3SYMb(k,x) + J_spt*J3SYMb(k,x) + K_spt*K3SYMb(k,x) + (KL_spt-K_spt)*L3SYMb(k,x));
}

double G3bdgp(double k[],double x[]){
  double norm=pow3(1./dnorm_spt);
  return norm*(-fdgp_spt*pow3(D_spt)*G3edsd(k,x)-Cdx_spt/H_spt*C3SYMb(k,x)
          -(F3d_spt/H_spt-D_spt*fdgp_spt*F_spt)*F3SYMb(k,x)-(Id_spt/H_spt-D_spt*Fd_spt/H_spt)*I3SYMb(k,x)
          - Jd_spt/H_spt*J3SYMb(k,x) - Kd_spt/H_spt*K3SYMb(k,x) -  (KLd_spt-Kd_spt)/H_spt*L3SYMb(k,x));
}


/// 4TH ORDER OHHH MA GAWD /////

double evol4[34];
double devol4[34];

// Arguments are k[0-3] = k1,k2,k3,k4
// x index index
// 01 = 1 | 02 = 2 | 03 = 5 | 12 = 3 | 13 = 0 | 23 = 4

double F4dgp(double k[], double x[]){

  // mags^2
  double k0s = pow2(k[0]);
  double k1s = pow2(k[1]);
  double k2s = pow2(k[2]);
  double k3s = pow2(k[3]);

  double k1234 = k0s+k1s+k2s + k3s
              + 2.*k[0]*k[1]*x[1] + 2.*k[0]*k[2]*x[2] + 2.*k[1]*k[2]*x[3]
              + 2.*k[0]*k[3]*x[5] + 2.*k[1]*k[3]*x[0] + 2.*k[2]*k[3]*x[4];
  double k123 = k0s+k1s+k2s + 2.*k[0]*k[1]*x[1] + 2.*k[0]*k[2]*x[2] + 2.*k[1]*k[2]*x[3];
  double k234 = k1s+k2s+k3s + 2.*k[1]*k[2]*x[3] + 2.*k[1]*k[3]*x[0] + 2.*k[2]*k[3]*x[4];
  double k341 = k2s+k3s+k0s + 2.*k[2]*k[3]*x[4] + 2.*k[2]*k[0]*x[2] + 2.*k[3]*k[0]*x[5];
  double k412 = k3s+k0s+k1s + 2.*k[3]*k[0]*x[5] + 2.*k[3]*k[1]*x[0] + 2.*k[0]*k[1]*x[1];
  double k12 = k0s+k1s+2.*k[0]*k[1]*x[1];
  double k23 = k1s+k2s+2.*k[1]*k[2]*x[3];
  double k34 = k2s+k3s+2.*k[2]*k[3]*x[4];
  double k41 = k3s+k0s+2.*k[3]*k[0]*x[5];
  double k13 = k0s+k2s+2.*k[0]*k[2]*x[2];
  double k24 = k1s+k3s+2.*k[1]*k[3]*x[0];

  // x index legend
  // 01 = 1 | 02 = 2 | 03 = 5 | 12 = 3 | 13 = 0 | 23 = 4

// proper magnitudes and angles
  double krt2[6];
  double xba[4];
  double krt[4];

// doublets
     krt2[0] = sqrt(k12);
     krt2[1] = sqrt(k23);
     krt2[2] = sqrt(k34);
     krt2[3] = sqrt(k41);
     krt2[4] = sqrt(k13);
     krt2[5] = sqrt(k24);

//triplets
   krt[0] = sqrt(k234);
   krt[1] = sqrt(k341);
   krt[2] = sqrt(k412);
   krt[3] = sqrt(k123);

// x index legend
   // 01 = 1 | 02 = 2 | 03 = 5 | 12 = 3 | 13 = 0 | 23 = 4

// angles between triplets and singlets
   xba[0] = (k[1]*x[1]+k[2]*x[2]+k[3]*x[5])/krt[0]; // k234.k1
   xba[1] = (k[2]*x[3]+k[3]*x[0]+k[0]*x[1])/krt[1]; // k341.k2
   xba[2] = (k[3]*x[4]+k[0]*x[2]+k[1]*x[3])/krt[2]; // k412.k3
   xba[3] = (k[0]*x[5]+k[1]*x[0]+k[2]*x[4])/krt[3]; // k123.k4


// basic kernel functions
  double abt[4][4];

// factor of 2 comes from beta(k_i,k_jkl)*(theta_3 theta_1 + theta_1theta_3) => cyclic perms are the same
  abt[0][0]= 2.*beta1(k[0],krt[0],xba[0]);
  abt[0][1]= 2.*beta1(k[1],krt[1],xba[1]);
  abt[0][2]= 2.*beta1(k[2],krt[2],xba[2]);
  abt[0][3]= 2.*beta1(k[3],krt[3],xba[3]);

// alpha(k_i,k_jkl)
  abt[1][0]= alpha(k[0],krt[0],xba[0]);
  abt[1][1]= alpha(k[1],krt[1],xba[1]);
  abt[1][2]= alpha(k[2],krt[2],xba[2]);
  abt[1][3]= alpha(k[3],krt[3],xba[3]);


// alpha(k_jkl,k_i)
  abt[2][0]= alpha(krt[0],k[0],xba[0]);
  abt[2][1]= alpha(krt[1],k[1],xba[1]);
  abt[2][2]= alpha(krt[2],k[2],xba[2]);
  abt[2][3]= alpha(krt[3],k[3],xba[3]);

// (1-u^2_i,jkl)
  abt[3][0]= 2.*ker1(xba[0]);
  abt[3][1]= 2.*ker1(xba[1]);
  abt[3][2]= 2.*ker1(xba[2]);
  abt[3][3]= 2.*ker1(xba[3]);

// Individual kernel arguments for feeding into C3,F3,I3,J3,K3,L3
  double karg[4][6];
  double xarg[4][3];

// x index legend
// 01 = 1 | 02 = 2 | 03 = 5 | 12 = 3 | 13 = 0 | 23 = 4

// k2,k3,k4, k34,k23,k24
  karg[0][0] = k[1];
  karg[0][1] = k[2];
  karg[0][2] = k[3];
  karg[0][3] = krt2[2];
  karg[0][4] = krt2[1];
  karg[0][5] = krt2[5];
  xarg[0][0] = x[4];
  xarg[0][1] = x[3];
  xarg[0][2] = x[0];

// k1,k3,k4, k34,k13,k14
  karg[1][0] = k[0];
  karg[1][1] = k[2];
  karg[1][2] = k[3];
  karg[1][3] = krt2[2];
  karg[1][4] = krt2[4];
  karg[1][5] = krt2[3];
  xarg[1][0] = x[4];
  xarg[1][1] = x[2];
  xarg[1][2] = x[5];

//k1,k2,k4,k24,k12,k14
  karg[2][0] = k[0];
  karg[2][1] = k[1];
  karg[2][2] = k[3];
  karg[2][3] = krt2[5];
  karg[2][4] = krt2[0];
  karg[2][5] = krt2[3];
  xarg[2][0] = x[0];
  xarg[2][1] = x[1];
  xarg[2][2] = x[5];

//k1,k2,k3, k23,k12,k13
  karg[3][0] = k[0];
  karg[3][1] = k[1];
  karg[3][2] = k[2];
  karg[3][3] = krt2[1];
  karg[3][4] = krt2[0];
  karg[3][5] = krt2[4];
  xarg[3][0] = x[3];
  xarg[3][1] = x[1];
  xarg[3][2] = x[2];


double terms[16];

for(int i=0;i<16;i++){
  terms[i]=0.;
}

// COMPUTE A,D,E,I

//A: beta(k123,k4) theta_dgp(k1,k2,k3) term
//D: alpha(k1,k234) delta_dgp(k2,k3,k4) term
//E: alpha(k234,k1) theta_dgp(k2,k3,k4) term
//I: gam2(k1,k234) delta_dgp(k2,k3,k4) term
double kernels[6];
double kargi[6],xargi[3];

for(int n = 0; n<4; n++){ // index to add all permutations
          for(int i=0;i<6;i++){
            kargi[i]=karg[n][i];
          }
          for(int j=0; j<3; j++){
            xargi[j]=xarg[n][j];
          }
          kernels[0] = C3SYMb(kargi,xargi);
          kernels[1] = F3SYMb(kargi,xargi);
          kernels[2] = I3SYMb(kargi,xargi);
          kernels[3] = J3SYMb(kargi,xargi);
          kernels[4] = L3SYMb(kargi,xargi);
          kernels[5] = K3SYMb(kargi,xargi);

              for(int m = 0; m<4; m++){ // index to choose A,D,E, or I
                  for(int o=0;o<6;o++){ // index to choose scale dep kernel (and evol factor )
                      terms[m] += evol4[6*m+o+5]*1./4.*abt[m][n]*kernels[o];
                    }
                  }
            }



//COMPUTE L

for(int i=0; i<4; i++){
  for(int j=0;j<6;j++){
    kargi[j]=karg[i][j];
  }
  for(int j=0; j<3; j++){
    xargi[j]=xarg[i][j];
    }
        terms[4] += evol4[31]*1./4.*abt[3][i]*F3edsd(kargi,xargi); // gam2(k1,k234) delta_gr(k2,k3,k4) term
    }

// 22 terms

double xba2[3];
double abt2[4][6];


   // x index
   // 01 = 1 | 02 = 2 | 03 = 5 | 12 = 3 | 13 = 0 | 23 = 4

// angle between k_ij and k_kl
   xba2[0] = (k[0]*k[2]*x[2]+k[0]*k[3]*x[5]+k[1]*k[2]*x[3]+k[1]*k[3]*x[0])/krt2[0]/krt2[2]; // k12.k34
   xba2[1] = (k[0]*k[1]*x[1]+k[0]*k[3]*x[5]+k[1]*k[2]*x[3]+k[2]*k[3]*x[4])/krt2[4]/krt2[5]; // k13.k24
   xba2[2] = (k[0]*k[1]*x[1]+k[1]*k[3]*x[0]+k[0]*k[2]*x[2]+k[2]*k[3]*x[4])/krt2[1]/krt2[3]; // k23.k14

  // permutations of basic kernels
  abt2[0][0]= beta1(krt2[0],krt2[2],xba2[0]);
  abt2[0][1]= beta1(krt2[4],krt2[5],xba2[1]);
  abt2[0][2]= beta1(krt2[1],krt2[3],xba2[2]);
  abt2[0][3]= abt2[0][1];
  abt2[0][4]= abt2[0][2];
  abt2[0][5]= abt2[0][0];


  // doublets

  abt2[1][0]= alpha(krt2[0],krt2[2],xba2[0]);
  abt2[1][1]= alpha(krt2[4],krt2[5],xba2[1]);
  abt2[1][2]= alpha(krt2[1],krt2[3],xba2[2]);
  abt2[1][3]= alpha(krt2[5],krt2[4],xba2[1]);
  abt2[1][4]= alpha(krt2[3],krt2[1],xba2[2]);
  abt2[1][5]= alpha(krt2[2],krt2[0],xba2[0]);

  abt2[2][0]= ker1(xba2[0]);
  abt2[2][1]= ker1(xba2[1]);
  abt2[2][2]= ker1(xba2[2]);
  abt2[2][3]= abt2[2][1];
  abt2[2][4]= abt2[2][2];
  abt2[2][5]= abt2[2][0];

  // Individual kernel arguments
  // 01 = 1 | 02 = 2 | 03 = 5 | 12 = 3 | 13 = 0 | 23 = 4

    double kernels2[6][6];

    kernels2[0][0] = F2eds(k[0],k[1],x[1]);
    kernels2[0][1] = F2eds(k[0],k[2],x[2]);
    kernels2[0][2] = F2eds(k[1],k[2],x[3]);
    kernels2[0][3] = F2eds(k[1],k[3],x[0]);
    kernels2[0][4] = F2eds(k[0],k[3],x[5]);
    kernels2[0][5] = F2eds(k[2],k[3],x[4]);

    kernels2[1][0] = G2eds(k[0],k[1],x[1]); //
    kernels2[1][1] = G2eds(k[0],k[2],x[2]); //
    kernels2[1][2] = G2eds(k[1],k[2],x[3]); //
    kernels2[1][3] = G2eds(k[1],k[3],x[0]); //
    kernels2[1][4] = G2eds(k[0],k[3],x[5]); //
    kernels2[1][5] = G2eds(k[2],k[3],x[4]); //

    kernels2[2][0] = ker1(x[1]);
    kernels2[2][1] = ker1(x[2]);
    kernels2[2][2] = ker1(x[3]);
    kernels2[2][3] = ker1(x[0]);
    kernels2[2][4] = ker1(x[5]);
    kernels2[2][5] = ker1(x[4]);

// INVERT POSITION

    kernels2[3][0] = kernels2[0][5];
    kernels2[3][1] = kernels2[0][3];
    kernels2[3][2] = kernels2[0][4];
    kernels2[3][3] = kernels2[0][1];
    kernels2[3][4] = kernels2[0][2];
    kernels2[3][5] = kernels2[0][0];

    kernels2[4][0] = kernels2[1][5];
    kernels2[4][1] = kernels2[1][3];
    kernels2[4][2] = kernels2[1][4];
    kernels2[4][3] = kernels2[1][1];
    kernels2[4][4] = kernels2[1][2];
    kernels2[4][5] = kernels2[1][0];

    kernels2[5][0] = kernels2[2][5]; //
    kernels2[5][1] = kernels2[2][3]; //
    kernels2[5][2] = kernels2[2][4]; //
    kernels2[5][3] = kernels2[2][1]; //
    kernels2[5][4] = kernels2[2][2]; //
    kernels2[5][5] = kernels2[2][0]; //


// COMPUTE B,C,F,G,H,J,K,M (22 terms)
//correction:
//I THINK THERE SHOULD BE A FACTOR OF 2 IN TERMS[5] HERE! NUMERICAL NEEDS TO BE CORRECTED TO ADJUST FOR THIS

// add all permutations
for(int m = 0; m<6; m++){
          terms[5]+=evol4[0]*1./6.*(abt2[0][m]*kernels2[1][m]*kernels2[5][m]); // B checked (theta_dgp(k1,k2) theta_gr(k3,k4) beta(k12,k34) terms)
          terms[6]+=evol4[1]*1./6.*(abt2[0][m]*kernels2[2][m]*kernels2[5][m]); // C checked (theta_dgp(k1,k2) theta_dgp(k3,k4) beta(k12,k34) terms)
          terms[7]+=evol4[2]*1./6.*(abt2[1][m]*kernels2[1][m]*kernels2[5][m]); // F checked (theta_gr(k1,k2) delta_dgp(k3,k4) alpha(k12,k34) terms)
          terms[8]+=evol4[3]*1./6.*(abt2[1][m]*kernels2[2][m]*kernels2[3][m]); // G checked (theta_dgp(k1,k2) delta_gr(k3,k4) alpha(k12,k34) terms)
          terms[9]+=evol4[4]*1./6.*(abt2[1][m]*kernels2[2][m]*kernels2[5][m]); // H checked (theta_dgp(k1,k2) delta_dgp(k3,k4) alpha(k12,k34) terms)
          terms[10]+=evol4[29]*1./6.*(abt2[2][m]*kernels2[0][m]*kernels2[5][m]); // J checked (delta_gr(k1,k2) delta_dgp(k3,k4) gam2(k12,k34)  + gam3(k12,k3,k4) delta_gr(k1,k2) terms)
          terms[11]+=evol4[30]*1./6.*(abt2[2][m]*kernels2[2][m]*kernels2[5][m]); // K checked (delta_dgp(k1,k2) delta_dgp(k3,k4) gam2(k12,k34) +gam3(k12,k3,k4) delta_dgp(k1,k2) + gam4b(k1,k2,k3,k4) terms)
          terms[12]+=evol4[31]*1./6.*(abt2[2][m]*kernels2[0][m]*kernels2[3][m]); // M checked (delta_gr(k1,k2) delta_gr(k3,k4) gam2(k12,k34)  terms)
   }


// REMAINING TERMS: N,O

double kernelsn[12];
// reminder of x indices
// 01 = 1 | 02 = 2 | 03 = 5 | 12 = 3 | 13 = 0 | 23 = 4

// gamma3 and gamma4 kernel components with arguments commented on left
kernelsn[0]=ker1((k[2]*x[3]+k[3]*x[0])/krt2[2]); // k2.k34
kernelsn[1]=ker1((k[1]*x[3]+k[3]*x[4])/krt2[5]); // k3.k24
kernelsn[2]=ker1((k[1]*x[0]+k[2]*x[4])/krt2[1]); // k4.k23
kernelsn[3]=ker1((k[2]*x[2]+k[3]*x[5])/krt2[2]); // k1.k34
kernelsn[4]=ker1((k[0]*x[2]+k[3]*x[4])/krt2[3]); // k3.k41
kernelsn[5]=ker1((k[0]*x[5]+k[2]*x[4])/krt2[4]); // k4.k13
kernelsn[6]=ker1((k[1]*x[1]+k[3]*x[5])/krt2[5]); // k1.k24
kernelsn[7]=ker1((k[0]*x[1]+k[3]*x[0])/krt2[3]); // k2.k41
kernelsn[8]=ker1((k[0]*x[5]+k[1]*x[0])/krt2[0]); // k4.k12
kernelsn[9]=ker1((k[1]*x[1]+k[2]*x[2])/krt2[1]); // k1.k23
kernelsn[10]=ker1((k[0]*x[1]+k[2]*x[3])/krt2[4]); // k2.k13
kernelsn[11]=ker1((k[0]*x[2]+k[1]*x[3])/krt2[0]); // k3.k12

// 01 = 1 | 02 = 2 | 03 = 5 | 12 = 3 | 13 = 0 | 23 = 4

double common =abt[3][0]/2.*(kernels2[2][5]*kernelsn[0] + kernels2[2][3]*kernelsn[1] + kernels2[2][2]*kernelsn[2])
              +abt[3][1]/2.*(kernels2[2][5]*kernelsn[3] + kernels2[2][4]*kernelsn[4] + kernels2[2][1]*kernelsn[5])
              +abt[3][2]/2.*(kernels2[2][3]*kernelsn[6] + kernels2[2][4]*kernelsn[7] + kernels2[2][0]*kernelsn[8])
              +abt[3][3]/2.*(kernels2[2][2]*kernelsn[9] + kernels2[2][1]*kernelsn[10] + kernels2[2][0]*kernelsn[11]);

// N
terms[13] = evol4[32]*1./12.*common; // (gam3(k1,k2,k34) delta_dgp(k3,k4) + gam3(k1,k23,k4) delta_dgp(k2,k3) + gam4a(k1,k2,k3,k4) terms)

// O
terms[14] = evol4[33]*1./12.*(abt[3][0]/2.*(kernels2[0][5]*kernelsn[0] + kernels2[0][3]*kernelsn[1] + kernels2[0][2]*kernelsn[2])
                            +abt[3][1]/2.*(kernels2[0][5]*kernelsn[3] + kernels2[0][4]*kernelsn[4] + kernels2[0][1]*kernelsn[5])
                            +abt[3][2]/2.*(kernels2[0][3]*kernelsn[6] + kernels2[0][4]*kernelsn[7] + kernels2[0][0]*kernelsn[8])
                            +abt[3][3]/2.*(kernels2[0][2]*kernelsn[9] + kernels2[0][1]*kernelsn[10] + kernels2[0][0]*kernelsn[11])); // (gam3(k1,k2,k34) delta_gr(k3,k4)  + gam3(k1,k23,k4) delta_gr(k2,k3)  terms)



terms[15]= pow4(D_spt)/396.*(
     27.*(k0s+k[0]*k[1]*x[1]+k[0]*k[2]*x[2]+k[0]*k[3]*x[5])/k0s*F3edsb(k[1],k[2],k[3],krt2[2],krt2[1],krt2[5],x[4],x[3],x[0])
    +(27.*(k234 + k[0]*(k[1]*x[1] + k[2]*x[2] + k[3]*x[5]))/k234 + 6.*k1234*k[0]*(k[1]*x[1]+k[2]*x[2]+k[3]*x[5])/k0s/k234)*G3edsb(k[1],k[2],k[3],krt2[2],krt2[1],krt2[5],x[4],x[3],x[0])

  +27.*(k1s+k[1]*k[2]*x[3]+k[1]*k[3]*x[0]+k[1]*k[0]*x[1])/k1s*F3edsb(k[2],k[3],k[0],krt2[3],krt2[2],krt2[4],x[5],x[4],x[2])
  +(27.*(k341 + k[1]*(k[2]*x[3] + k[3]*x[0] + k[0]*x[1]))/k341 + 6.*k1234*k[1]*(k[2]*x[3]+k[3]*x[0]+k[0]*x[1])/k1s/k341)*G3edsb(k[2],k[3],k[0],krt2[3],krt2[2],krt2[4],x[5],x[4],x[2])

  +27.*(k2s+k[2]*k[3]*x[4]+k[2]*k[0]*x[2]+k[2]*k[1]*x[3])/k2s*F3edsb(k[3],k[0],k[1],krt2[0],krt2[3],krt2[5],x[1],x[5],x[0])
  +(27.*(k412 + k[2]*(k[3]*x[4] + k[0]*x[2] + k[1]*x[3]))/k412 + 6.*k1234*k[2]*(k[3]*x[4]+k[0]*x[2]+k[1]*x[3])/k2s/k412)*G3edsb(k[3],k[0],k[1],krt2[0],krt2[3],krt2[5],x[1],x[5],x[0])

  +27.*(k3s+k[3]*k[0]*x[5]+k[3]*k[1]*x[0]+k[3]*k[2]*x[4])/k3s*F3edsb(k[0],k[1],k[2],krt2[1],krt2[0],krt2[4],x[3],x[1],x[2])
  +(27.*(k123 + k[3]*(k[0]*x[5] + k[1]*x[0] + k[2]*x[4]))/k123 + 6.*k1234*k[3]*(k[0]*x[5]+k[1]*x[0]+k[2]*x[4])/k3s/k123)*G3edsb(k[0],k[1],k[2],krt2[1],krt2[0],krt2[4],x[3],x[1],x[2])

  +18.*(k12+k[2]*k[0]*x[2]+k[2]*k[1]*x[3]+k[3]*k[0]*x[5]+k[3]*k[1]*x[0])/k12*kernels2[1][0]*kernels2[0][5]
  +18.*(k23+k[3]*k[1]*x[0]+k[3]*k[2]*x[4]+k[0]*k[1]*x[1]+k[0]*k[2]*x[2])/k23*kernels2[1][2]*kernels2[0][4]
  +18.*(k34+k[0]*k[2]*x[2]+k[0]*k[3]*x[5]+k[1]*k[2]*x[3]+k[1]*k[3]*x[0])/k34*kernels2[1][5]*kernels2[0][0]
  +18.*(k13+k[0]*k[1]*x[1]+k[0]*k[3]*x[5]+k[1]*k[2]*x[3]+k[2]*k[3]*x[4])/k13*kernels2[1][1]*kernels2[0][3]
  +18.*(k24+k[1]*k[2]*x[3]+k[1]*k[0]*x[1]+k[3]*k[0]*x[5]+k[2]*k[3]*x[4])/k24*kernels2[1][3]*kernels2[0][1]
  +18.*(k41+k[0]*k[1]*x[1]+k[0]*k[2]*x[2]+k[3]*k[1]*x[0]+k[2]*k[3]*x[4])/k41*kernels2[1][4]*kernels2[0][2]

  +4.*k1234*(k[0]*k[2]*x[2]+k[0]*k[3]*x[5]+k[1]*k[2]*x[3]+k[1]*k[3]*x[0])/k12/k34*kernels2[1][0]*kernels2[1][5]
  +4.*k1234*(k[1]*k[3]*x[0]+k[1]*k[0]*x[1]+k[2]*k[3]*x[4]+k[2]*k[0]*x[2])/k23/k41*kernels2[1][2]*kernels2[1][4]
  +4.*k1234*(k[1]*k[0]*x[1]+k[1]*k[2]*x[3]+k[3]*k[0]*x[5]+k[3]*k[2]*x[4])/k24/k13*kernels2[1][3]*kernels2[1][1]); // checked

double total = 0.;

for(int i=0; i<16; i++){
  total += terms[i];
}
  return pow4(1./dnorm_spt)*(total);

}



/// Velocity kernel
double G4dgp(double k[], double x[]){

  // mags^2
  double k0s = pow2(k[0]);
  double k1s = pow2(k[1]);
  double k2s = pow2(k[2]);
  double k3s = pow2(k[3]);

  double k1234 = k0s+k1s+k2s + k3s
              + 2.*k[0]*k[1]*x[1] + 2.*k[0]*k[2]*x[2] + 2.*k[1]*k[2]*x[3]
              + 2.*k[0]*k[3]*x[5] + 2.*k[1]*k[3]*x[0] + 2.*k[2]*k[3]*x[4];
  double k123 = k0s+k1s+k2s + 2.*k[0]*k[1]*x[1] + 2.*k[0]*k[2]*x[2] + 2.*k[1]*k[2]*x[3];
  double k234 = k1s+k2s+k3s + 2.*k[1]*k[2]*x[3] + 2.*k[1]*k[3]*x[0] + 2.*k[2]*k[3]*x[4];
  double k341 = k2s+k3s+k0s + 2.*k[2]*k[3]*x[4] + 2.*k[2]*k[0]*x[2] + 2.*k[3]*k[0]*x[5];
  double k412 = k3s+k0s+k1s + 2.*k[3]*k[0]*x[5] + 2.*k[3]*k[1]*x[0] + 2.*k[0]*k[1]*x[1];
  double k12 = k0s+k1s+2.*k[0]*k[1]*x[1];
  double k23 = k1s+k2s+2.*k[1]*k[2]*x[3];
  double k34 = k2s+k3s+2.*k[2]*k[3]*x[4];
  double k41 = k3s+k0s+2.*k[3]*k[0]*x[5];
  double k13 = k0s+k2s+2.*k[0]*k[2]*x[2];
  double k24 = k1s+k3s+2.*k[1]*k[3]*x[0];

  // x index legend
  // 01 = 1 | 02 = 2 | 03 = 5 | 12 = 3 | 13 = 0 | 23 = 4

// proper magnitudes and angles
  double krt2[6];
  double xba[4];
  double krt[4];

// doublets
     krt2[0] = sqrt(k12);
     krt2[1] = sqrt(k23);
     krt2[2] = sqrt(k34);
     krt2[3] = sqrt(k41);
     krt2[4] = sqrt(k13);
     krt2[5] = sqrt(k24);

//triplets
   krt[0] = sqrt(k234);
   krt[1] = sqrt(k341);
   krt[2] = sqrt(k412);
   krt[3] = sqrt(k123);

// x index legend
   // 01 = 1 | 02 = 2 | 03 = 5 | 12 = 3 | 13 = 0 | 23 = 4

// angles between triplets and singlets
   xba[0] = (k[1]*x[1]+k[2]*x[2]+k[3]*x[5])/krt[0]; // k234.k1
   xba[1] = (k[2]*x[3]+k[3]*x[0]+k[0]*x[1])/krt[1]; // k341.k2
   xba[2] = (k[3]*x[4]+k[0]*x[2]+k[1]*x[3])/krt[2]; // k412.k3
   xba[3] = (k[0]*x[5]+k[1]*x[0]+k[2]*x[4])/krt[3]; // k123.k4


// basic kernel functions
  double abt[4][4];

// factor of 2 comes from beta(k_i,k_jkl)*(theta_3 theta_1 + theta_1theta_3) => cyclic perms are the same
  abt[0][0]= 2.*beta1(k[0],krt[0],xba[0]);
  abt[0][1]= 2.*beta1(k[1],krt[1],xba[1]);
  abt[0][2]= 2.*beta1(k[2],krt[2],xba[2]);
  abt[0][3]= 2.*beta1(k[3],krt[3],xba[3]);

// alpha(k_i,k_jkl)
  abt[1][0]= alpha(k[0],krt[0],xba[0]);
  abt[1][1]= alpha(k[1],krt[1],xba[1]);
  abt[1][2]= alpha(k[2],krt[2],xba[2]);
  abt[1][3]= alpha(k[3],krt[3],xba[3]);


// alpha(k_jkl,k_i)
  abt[2][0]= alpha(krt[0],k[0],xba[0]);
  abt[2][1]= alpha(krt[1],k[1],xba[1]);
  abt[2][2]= alpha(krt[2],k[2],xba[2]);
  abt[2][3]= alpha(krt[3],k[3],xba[3]);

// (1-u^2_i,jkl)
  abt[3][0]= 2.*ker1(xba[0]);
  abt[3][1]= 2.*ker1(xba[1]);
  abt[3][2]= 2.*ker1(xba[2]);
  abt[3][3]= 2.*ker1(xba[3]);

// Individual kernel arguments for feeding into C3,F3,I3,J3,K3,L3
  double karg[4][6];
  double xarg[4][3];

// x index legend
// 01 = 1 | 02 = 2 | 03 = 5 | 12 = 3 | 13 = 0 | 23 = 4

// k2,k3,k4, k34,k23,k24
  karg[0][0] = k[1];
  karg[0][1] = k[2];
  karg[0][2] = k[3];
  karg[0][3] = krt2[2];
  karg[0][4] = krt2[1];
  karg[0][5] = krt2[5];
  xarg[0][0] = x[4];
  xarg[0][1] = x[3];
  xarg[0][2] = x[0];

// k1,k3,k4, k34,k13,k14
  karg[1][0] = k[0];
  karg[1][1] = k[2];
  karg[1][2] = k[3];
  karg[1][3] = krt2[2];
  karg[1][4] = krt2[4];
  karg[1][5] = krt2[3];
  xarg[1][0] = x[4];
  xarg[1][1] = x[2];
  xarg[1][2] = x[5];

//k1,k2,k4,k24,k12,k14
  karg[2][0] = k[0];
  karg[2][1] = k[1];
  karg[2][2] = k[3];
  karg[2][3] = krt2[5];
  karg[2][4] = krt2[0];
  karg[2][5] = krt2[3];
  xarg[2][0] = x[0];
  xarg[2][1] = x[1];
  xarg[2][2] = x[5];

//k1,k2,k3, k23,k12,k13
  karg[3][0] = k[0];
  karg[3][1] = k[1];
  karg[3][2] = k[2];
  karg[3][3] = krt2[1];
  karg[3][4] = krt2[0];
  karg[3][5] = krt2[4];
  xarg[3][0] = x[3];
  xarg[3][1] = x[1];
  xarg[3][2] = x[2];


double terms[16];

for(int i=0;i<16;i++){
  terms[i]=0.;
}

// COMPUTE A,D,E,I

//A: beta(k123,k4) theta_dgp(k1,k2,k3) term
//D: alpha(k1,k234) delta_dgp(k2,k3,k4) term
//E: alpha(k234,k1) theta_dgp(k2,k3,k4) term
//I: gam2(k1,k234) delta_dgp(k2,k3,k4) term
double kernels[6];
double kargi[6],xargi[3];

for(int n = 0; n<4; n++){ // index to add all permutations
          for(int i=0;i<6;i++){
            kargi[i]=karg[n][i];
          }
          for(int j=0; j<3; j++){
            xargi[j]=xarg[n][j];
          }
          kernels[0] = C3SYMb(kargi,xargi);
          kernels[1] = F3SYMb(kargi,xargi);
          kernels[2] = I3SYMb(kargi,xargi);
          kernels[3] = J3SYMb(kargi,xargi);
          kernels[4] = L3SYMb(kargi,xargi);
          kernels[5] = K3SYMb(kargi,xargi);

              for(int m = 0; m<4; m++){ // index to choose A,D,E, or I
                  for(int o=0;o<6;o++){ // index to choose scale dep kernel (and evol factor )
                      terms[m] += devol4[6*m+o+5]*1./4.*abt[m][n]*kernels[o];
                    }
                  }
            }



//COMPUTE L

for(int i=0; i<4; i++){
  for(int j=0;j<6;j++){
    kargi[j]=karg[i][j];
  }
  for(int j=0; j<3; j++){
    xargi[j]=xarg[i][j];
    }
        terms[4] += devol4[31]*1./4.*abt[3][i]*F3edsd(kargi,xargi); // gam2(k1,k234) delta_gr(k2,k3,k4) term
    }

// 22 terms

double xba2[3];
double abt2[4][6];


   // x index
   // 01 = 1 | 02 = 2 | 03 = 5 | 12 = 3 | 13 = 0 | 23 = 4

// angle between k_ij and k_kl
   xba2[0] = (k[0]*k[2]*x[2]+k[0]*k[3]*x[5]+k[1]*k[2]*x[3]+k[1]*k[3]*x[0])/krt2[0]/krt2[2]; // k12.k34
   xba2[1] = (k[0]*k[1]*x[1]+k[0]*k[3]*x[5]+k[1]*k[2]*x[3]+k[2]*k[3]*x[4])/krt2[4]/krt2[5]; // k13.k24
   xba2[2] = (k[0]*k[1]*x[1]+k[1]*k[3]*x[0]+k[0]*k[2]*x[2]+k[2]*k[3]*x[4])/krt2[1]/krt2[3]; // k23.k14

  // permutations of basic kernels
  abt2[0][0]= beta1(krt2[0],krt2[2],xba2[0]);
  abt2[0][1]= beta1(krt2[4],krt2[5],xba2[1]);
  abt2[0][2]= beta1(krt2[1],krt2[3],xba2[2]);
  abt2[0][3]= abt2[0][1];
  abt2[0][4]= abt2[0][2];
  abt2[0][5]= abt2[0][0];


  // doublets

  abt2[1][0]= alpha(krt2[0],krt2[2],xba2[0]);
  abt2[1][1]= alpha(krt2[4],krt2[5],xba2[1]);
  abt2[1][2]= alpha(krt2[1],krt2[3],xba2[2]);
  abt2[1][3]= alpha(krt2[5],krt2[4],xba2[1]);
  abt2[1][4]= alpha(krt2[3],krt2[1],xba2[2]);
  abt2[1][5]= alpha(krt2[2],krt2[0],xba2[0]);

  abt2[2][0]= ker1(xba2[0]);
  abt2[2][1]= ker1(xba2[1]);
  abt2[2][2]= ker1(xba2[2]);
  abt2[2][3]= abt2[2][1];
  abt2[2][4]= abt2[2][2];
  abt2[2][5]= abt2[2][0];

  // Individual kernel arguments
  // 01 = 1 | 02 = 2 | 03 = 5 | 12 = 3 | 13 = 0 | 23 = 4

    double kernels2[6][6];

    kernels2[0][0] = F2eds(k[0],k[1],x[1]);
    kernels2[0][1] = F2eds(k[0],k[2],x[2]);
    kernels2[0][2] = F2eds(k[1],k[2],x[3]);
    kernels2[0][3] = F2eds(k[1],k[3],x[0]);
    kernels2[0][4] = F2eds(k[0],k[3],x[5]);
    kernels2[0][5] = F2eds(k[2],k[3],x[4]);

    kernels2[1][0] = G2eds(k[0],k[1],x[1]); //
    kernels2[1][1] = G2eds(k[0],k[2],x[2]); //
    kernels2[1][2] = G2eds(k[1],k[2],x[3]); //
    kernels2[1][3] = G2eds(k[1],k[3],x[0]); //
    kernels2[1][4] = G2eds(k[0],k[3],x[5]); //
    kernels2[1][5] = G2eds(k[2],k[3],x[4]); //

    kernels2[2][0] = ker1(x[1]);
    kernels2[2][1] = ker1(x[2]);
    kernels2[2][2] = ker1(x[3]);
    kernels2[2][3] = ker1(x[0]);
    kernels2[2][4] = ker1(x[5]);
    kernels2[2][5] = ker1(x[4]);

// INVERT POSITION

    kernels2[3][0] = kernels2[0][5];
    kernels2[3][1] = kernels2[0][3];
    kernels2[3][2] = kernels2[0][4];
    kernels2[3][3] = kernels2[0][1];
    kernels2[3][4] = kernels2[0][2];
    kernels2[3][5] = kernels2[0][0];

    kernels2[4][0] = kernels2[1][5];
    kernels2[4][1] = kernels2[1][3];
    kernels2[4][2] = kernels2[1][4];
    kernels2[4][3] = kernels2[1][1];
    kernels2[4][4] = kernels2[1][2];
    kernels2[4][5] = kernels2[1][0];

    kernels2[5][0] = kernels2[2][5]; //
    kernels2[5][1] = kernels2[2][3]; //
    kernels2[5][2] = kernels2[2][4]; //
    kernels2[5][3] = kernels2[2][1]; //
    kernels2[5][4] = kernels2[2][2]; //
    kernels2[5][5] = kernels2[2][0]; //


// COMPUTE B,C,F,G,H,J,K,M (22 terms)
//correction:
//I THINK THERE SHOULD BE A FACTOR OF 2 IN TERMS[5] HERE! NUMERICAL NEEDS TO BE CORRECTED TO ADJUST FOR THIS

// add all permutations
for(int m = 0; m<6; m++){
          terms[5]+=devol4[0]*1./6.*(abt2[0][m]*kernels2[1][m]*kernels2[5][m]); // B checked (theta_dgp(k1,k2) theta_gr(k3,k4) beta(k12,k34) terms)
          terms[6]+=devol4[1]*1./6.*(abt2[0][m]*kernels2[2][m]*kernels2[5][m]); // C checked (theta_dgp(k1,k2) theta_dgp(k3,k4) beta(k12,k34) terms)
          terms[7]+=devol4[2]*1./6.*(abt2[1][m]*kernels2[1][m]*kernels2[5][m]); // F checked (theta_gr(k1,k2) delta_dgp(k3,k4) alpha(k12,k34) terms)
          terms[8]+=devol4[3]*1./6.*(abt2[1][m]*kernels2[2][m]*kernels2[3][m]); // G checked (theta_dgp(k1,k2) delta_gr(k3,k4) alpha(k12,k34) terms)
          terms[9]+=devol4[4]*1./6.*(abt2[1][m]*kernels2[2][m]*kernels2[5][m]); // H checked (theta_dgp(k1,k2) delta_dgp(k3,k4) alpha(k12,k34) terms)
          terms[10]+=devol4[29]*1./6.*(abt2[2][m]*kernels2[0][m]*kernels2[5][m]); // J checked (delta_gr(k1,k2) delta_dgp(k3,k4) gam2(k12,k34)  + gam3(k12,k3,k4) delta_gr(k1,k2) terms)
          terms[11]+=devol4[30]*1./6.*(abt2[2][m]*kernels2[2][m]*kernels2[5][m]); // K checked (delta_dgp(k1,k2) delta_dgp(k3,k4) gam2(k12,k34) +gam3(k12,k3,k4) delta_dgp(k1,k2) + gam4b(k1,k2,k3,k4) terms)
          terms[12]+=devol4[31]*1./6.*(abt2[2][m]*kernels2[0][m]*kernels2[3][m]); // M checked (delta_gr(k1,k2) delta_gr(k3,k4) gam2(k12,k34)  terms)
   }


// REMAINING TERMS: N,O

double kernelsn[12];
// reminder of x indices
// 01 = 1 | 02 = 2 | 03 = 5 | 12 = 3 | 13 = 0 | 23 = 4

// gamma3 and gamma4 kernel components with arguments commented on left
kernelsn[0]=ker1((k[2]*x[3]+k[3]*x[0])/krt2[2]); // k2.k34
kernelsn[1]=ker1((k[1]*x[3]+k[3]*x[4])/krt2[5]); // k3.k24
kernelsn[2]=ker1((k[1]*x[0]+k[2]*x[4])/krt2[1]); // k4.k23
kernelsn[3]=ker1((k[2]*x[2]+k[3]*x[5])/krt2[2]); // k1.k34
kernelsn[4]=ker1((k[0]*x[2]+k[3]*x[4])/krt2[3]); // k3.k41
kernelsn[5]=ker1((k[0]*x[5]+k[2]*x[4])/krt2[4]); // k4.k13
kernelsn[6]=ker1((k[1]*x[1]+k[3]*x[5])/krt2[5]); // k1.k24
kernelsn[7]=ker1((k[0]*x[1]+k[3]*x[0])/krt2[3]); // k2.k41
kernelsn[8]=ker1((k[0]*x[5]+k[1]*x[0])/krt2[0]); // k4.k12
kernelsn[9]=ker1((k[1]*x[1]+k[2]*x[2])/krt2[1]); // k1.k23
kernelsn[10]=ker1((k[0]*x[1]+k[2]*x[3])/krt2[4]); // k2.k13
kernelsn[11]=ker1((k[0]*x[2]+k[1]*x[3])/krt2[0]); // k3.k12

// 01 = 1 | 02 = 2 | 03 = 5 | 12 = 3 | 13 = 0 | 23 = 4

double common =abt[3][0]/2.*(kernels2[2][5]*kernelsn[0] + kernels2[2][3]*kernelsn[1] + kernels2[2][2]*kernelsn[2])
              +abt[3][1]/2.*(kernels2[2][5]*kernelsn[3] + kernels2[2][4]*kernelsn[4] + kernels2[2][1]*kernelsn[5])
              +abt[3][2]/2.*(kernels2[2][3]*kernelsn[6] + kernels2[2][4]*kernelsn[7] + kernels2[2][0]*kernelsn[8])
              +abt[3][3]/2.*(kernels2[2][2]*kernelsn[9] + kernels2[2][1]*kernelsn[10] + kernels2[2][0]*kernelsn[11]);

// N
terms[13] = devol4[32]*1./12.*common; // (gam3(k1,k2,k34) delta_dgp(k3,k4) + gam3(k1,k23,k4) delta_dgp(k2,k3) + gam4a(k1,k2,k3,k4) terms)

// O
terms[14] = devol4[33]*1./12.*(abt[3][0]/2.*(kernels2[0][5]*kernelsn[0] + kernels2[0][3]*kernelsn[1] + kernels2[0][2]*kernelsn[2])
                            +abt[3][1]/2.*(kernels2[0][5]*kernelsn[3] + kernels2[0][4]*kernelsn[4] + kernels2[0][1]*kernelsn[5])
                            +abt[3][2]/2.*(kernels2[0][3]*kernelsn[6] + kernels2[0][4]*kernelsn[7] + kernels2[0][0]*kernelsn[8])
                            +abt[3][3]/2.*(kernels2[0][2]*kernelsn[9] + kernels2[0][1]*kernelsn[10] + kernels2[0][0]*kernelsn[11])); // (gam3(k1,k2,k34) delta_gr(k3,k4)  + gam3(k1,k23,k4) delta_gr(k2,k3)  terms)



terms[15]= fdgp_spt*pow4(D_spt)/396.*(
     9.*(k0s+k[0]*k[1]*x[1]+k[0]*k[2]*x[2]+k[0]*k[3]*x[5])/k0s*F3edsb(k[1],k[2],k[3],krt2[2],krt2[1],krt2[5],x[4],x[3],x[0])
    +(9.*(k234 + k[0]*(k[1]*x[1] + k[2]*x[2] + k[3]*x[5]))/k234 + 24.*k1234*k[0]*(k[1]*x[1]+k[2]*x[2]+k[3]*x[5])/k0s/k234)*G3edsb(k[1],k[2],k[3],krt2[2],krt2[1],krt2[5],x[4],x[3],x[0])

  +9.*(k1s+k[1]*k[2]*x[3]+k[1]*k[3]*x[0]+k[1]*k[0]*x[1])/k1s*F3edsb(k[2],k[3],k[0],krt2[3],krt2[2],krt2[4],x[5],x[4],x[2])
  +(9.*(k341 + k[1]*(k[2]*x[3] + k[3]*x[0] + k[0]*x[1]))/k341 + 24.*k1234*k[1]*(k[2]*x[3]+k[3]*x[0]+k[0]*x[1])/k1s/k341)*G3edsb(k[2],k[3],k[0],krt2[3],krt2[2],krt2[4],x[5],x[4],x[2])

  +9.*(k2s+k[2]*k[3]*x[4]+k[2]*k[0]*x[2]+k[2]*k[1]*x[3])/k2s*F3edsb(k[3],k[0],k[1],krt2[0],krt2[3],krt2[5],x[1],x[5],x[0])
  +(9.*(k412 + k[2]*(k[3]*x[4] + k[0]*x[2] + k[1]*x[3]))/k412 + 24.*k1234*k[2]*(k[3]*x[4]+k[0]*x[2]+k[1]*x[3])/k2s/k412)*G3edsb(k[3],k[0],k[1],krt2[0],krt2[3],krt2[5],x[1],x[5],x[0])

  +9.*(k3s+k[3]*k[0]*x[5]+k[3]*k[1]*x[0]+k[3]*k[2]*x[4])/k3s*F3edsb(k[0],k[1],k[2],krt2[1],krt2[0],krt2[4],x[3],x[1],x[2])
  +(9.*(k123 + k[3]*(k[0]*x[5] + k[1]*x[0] + k[2]*x[4]))/k123 + 24.*k1234*k[3]*(k[0]*x[5]+k[1]*x[0]+k[2]*x[4])/k3s/k123)*G3edsb(k[0],k[1],k[2],krt2[1],krt2[0],krt2[4],x[3],x[1],x[2])

  +6.*(k12+k[2]*k[0]*x[2]+k[2]*k[1]*x[3]+k[3]*k[0]*x[5]+k[3]*k[1]*x[0])/k12*kernels2[1][0]*kernels2[0][5]
  +6.*(k23+k[3]*k[1]*x[0]+k[3]*k[2]*x[4]+k[0]*k[1]*x[1]+k[0]*k[2]*x[2])/k23*kernels2[1][2]*kernels2[0][4]
  +6.*(k34+k[0]*k[2]*x[2]+k[0]*k[3]*x[5]+k[1]*k[2]*x[3]+k[1]*k[3]*x[0])/k34*kernels2[1][5]*kernels2[0][0]
  +6.*(k13+k[0]*k[1]*x[1]+k[0]*k[3]*x[5]+k[1]*k[2]*x[3]+k[2]*k[3]*x[4])/k13*kernels2[1][1]*kernels2[0][3]
  +6.*(k24+k[1]*k[2]*x[3]+k[1]*k[0]*x[1]+k[3]*k[0]*x[5]+k[2]*k[3]*x[4])/k24*kernels2[1][3]*kernels2[0][1]
  +6.*(k41+k[0]*k[1]*x[1]+k[0]*k[2]*x[2]+k[3]*k[1]*x[0]+k[2]*k[3]*x[4])/k41*kernels2[1][4]*kernels2[0][2]

  +16.*k1234*(k[0]*k[2]*x[2]+k[0]*k[3]*x[5]+k[1]*k[2]*x[3]+k[1]*k[3]*x[0])/k12/k34*kernels2[1][0]*kernels2[1][5]
  +16.*k1234*(k[1]*k[3]*x[0]+k[1]*k[0]*x[1]+k[2]*k[3]*x[4]+k[2]*k[0]*x[2])/k23/k41*kernels2[1][2]*kernels2[1][4]
  +16.*k1234*(k[1]*k[0]*x[1]+k[1]*k[2]*x[3]+k[3]*k[0]*x[5]+k[3]*k[2]*x[4])/k24/k13*kernels2[1][3]*kernels2[1][1]); // checked


double total = 0.;

for(int i=0; i<16; i++){
  total += terms[i];
}
  return -pow4(1./dnorm_spt)*(total);

}

/// Fitting f2 kernels of 1805.10567 calibrated to GR sims as in the same paper

static double QVAL(double x){
  return (4.-pow(2.,x))/(1.+pow(2.,x+1.));
}
double F2fit(double vars[], double k1, double k2, double x){
  double a1,a2,a3,a4,a5,a6,a7,a8,a9,neff1,neff2,s8,knl,kap,lam;
  s8 = vars[0];
  knl = vars[1];
  kap = vars[2];
  lam = vars[3];
  neff1 =vars[4];
  neff2 = vars[5];

  a1=0.484;
  a2=3.74;
  a3=-0.849;
  a4=0.392;
  a5=1.013;
  a6=-0.575;
  a7=0.128;
  a8=-0.722;
  a9=-0.926;

  double val1=pow(k1/knl*a1,neff1+a2);
  double val2a=pow(k1/knl*a7,neff1+3.+a8);
  double val2b=pow(k1/knl*a7,neff1+3.5+a8);
  double val3a=pow(k1/knl*a5,neff1+3.+a9);
  double val3b=pow(k1/knl*a5,neff1+3.5+a9);

  double vval1=pow(k2/knl*a1,neff2+a2);
  double vval2a=pow(k2/knl*a7,neff2+3.+a8);
  double vval2b=pow(k2/knl*a7,neff2+3.5+a8);
  double vval3a=pow(k2/knl*a5,neff2+3.+a9);
  double vval3b=pow(k2/knl*a5,neff2+3.5+a9);


    double aterm1 = (1.+ pow(s8,a6)*sqrt(0.7*QVAL(neff1))*val1) /(1.+val1);
    double bterm1 = (1.+ 0.2*a3*(neff1+3.)*val2a) /(1.+val2b);
    double cterm1 = (1.+ 4.5*a4/(1.5+pow4(neff1+3.)) * val3a) /(1.+val3b);

    double aterm2 = (1.+ pow(s8,a6)*sqrt(0.7*QVAL(neff2))*vval1) /(1.+vval1);
    double bterm2 = (1.+ 0.2*a3*(neff2+3.)*vval2a) /(1.+vval2b);
    double cterm2 = (1.+ 4.5*a4/(1.5+pow4(neff2+3.)) * vval3a) /(1.+vval3b);

  return (kap-2./7.*lam)*aterm1*aterm2 + kap*x*(k1*k1+k2*k2)/(2.*k1*k2)*bterm1*bterm2 + lam*2./7.*cterm1*cterm2*pow2(x);

}


//scoccimaro's version
double F2fitsc(double vars[], double k1, double k2, double x){
  double a1,a2,a3,a4,a5,a6,a7,a8,a9,neff1,neff2,s8,knl,kap,lam;
  s8 = vars[0];
  knl = vars[1];
  kap = vars[2];
  lam = vars[3];
  neff1 =vars[4];
  neff2 = vars[5];

  a1=0.25;
  a2=3.5;
  a3=2.;
  a4=1.;
  a5=2.;
  a6=-0.2;
  a7=1.;
  a8=0.;
  a9=0.;

  double val1=pow(k1/knl*a1,neff1+a2);
  double val2a=pow(k1/knl*a7,neff1+3.+a8);
  double val2b=pow(k1/knl*a7,neff1+3.5+a8);
  double val3a=pow(k1/knl*a5,neff1+3.+a9);
  double val3b=pow(k1/knl*a5,neff1+3.5+a9);

  double vval1=pow(k2/knl*a1,neff2+a2);
  double vval2a=pow(k2/knl*a7,neff2+3.+a8);
  double vval2b=pow(k2/knl*a7,neff2+3.5+a8);
  double vval3a=pow(k2/knl*a5,neff2+3.+a9);
  double vval3b=pow(k2/knl*a5,neff2+3.5+a9);


    double aterm1 = (1.+pow(s8,a6)*sqrt(0.7*QVAL(neff1))*val1)/(1.+val1);
    double bterm1 = (1.+0.2*a3*(neff1+3.)*val2a)/(1.+val2b);
    double cterm1 = (1.+ 4.5*a4/(1.5+pow4(neff1+3.)) * val3a)/(1.+val3b);

    double aterm2 = (1.+pow(s8,a6)*sqrt(0.7*QVAL(neff2))*vval1)/(1.+vval1);
    double bterm2 = (1.+0.2*a3*(neff2+3.)*vval2a)/(1.+vval2b);
    double cterm2 = (1.+ 4.5*a4/(1.5+pow4(neff2+3.)) * vval3a)/(1.+vval3b);

  return  (kap-2./7.*lam)*aterm1*aterm2 + kap*x*(k1*k1+k2*k2)/(2.*k1*k2)*bterm1*bterm2 + lam*2./7.*cterm1*cterm2*pow2(x);

}


/* Differential equations for separable and non-separable cases */

///SEPARABLE/////

/*cosmological parameters for evolution factors (applied to analytic solutions) */

struct param_type2 {
	real omega00;
  real par1;
	real par2;
	real par3;
  int par4;
};

/* Evolution factor equations for separable DGP and LCDM solutions */
int func(double a, const double G[], double F[], void *params){
	param_type2 p = *(param_type2 *)(params);
	real omega0 = p.omega00;
	real p1 = p.par1;
	real p2 = p.par2;
	real p3 = p.par3;

  double hub1 = pow2(HA(a,omega0));
  double ap5 = pow(a,5);

	// D
	F[0] = G[1];
	F[1] = -1./a*(3.+HA1(a,omega0)/(hub1))*G[1]+3./2.*omega0/(hub1*ap5)*mu(a,1.,omega0,p1,p2,p3)*G[0];

	//E
	F[2] = G[3];
	F[3] = -1./a*(3.+HA1(a,omega0)/(hub1))*G[3]+3./2.*omega0/(hub1*ap5)*mu(a,1.,omega0,p1,p2,p3)*G[2]+G[1]*G[1]+3./2.*omega0/(hub1*ap5)*mu(a,1.,omega0,p1,p2,p3)*G[0]*G[0];

	//F
	F[4] = G[5];
	F[5] = -1./a*(3.+HA1(a,omega0)/(hub1))*G[5]+3./2.*omega0/(hub1*ap5)*mu(a,1.,omega0,p1,p2,p3)*G[4]-1./(a*a*hub1*24.*pow(beta(a,omega0,p1),3)*p1) *pow(omega0/(a*a*a),2)*G[0]*G[0];

	//C
	F[6] =G[7];
	F[7] = -1./a*(3.+HA1(a,omega0)/(hub1))*G[7]+3./2.*omega0/(hub1*ap5)*mu(a,1.,omega0,p1,p2,p3)*G[6]+G[1]*(G[5]-G[1]);

	//I
	F[8] = G[9] ;
	F[9] = -1./a*(3.+HA1(a,omega0)/(hub1))*G[9]+3./2.*omega0/(hub1*ap5)*mu(a,1.,omega0,p1,p2,p3)*G[8]+3./2.*omega0/(hub1*ap5)*mu(a,1.,omega0,p1,p2,p3)*(G[4]-G[0])*G[0]-1./(a*a*hub1*24.*pow(beta(a,omega0,p1),3)*p1)*pow(omega0/(a*a*a),2)*G[0]*G[0]*G[0]+G[1]*(G[5]-G[1]);


	//KL
	F[10] = G[11] ;
	F[11] = -1./a*(3.+HA1(a,omega0)/(hub1))*G[11]+3./2.*omega0/(hub1*ap5)*mu(a,1.,omega0,p1,p2,p3)*G[10]
	-1./(a*a*hub1*12.*pow(beta(a,omega0,p1),3)*p1)*pow(omega0/(a*a*a),2)*G[0]*(G[2]+G[4]-2*G[0])
	+1./(a*a*hub1*144.*pow(beta(a,omega0,p1),5)*pow(p1,2))*pow(omega0/(a*a*a),3)*G[0]*G[0]*G[0];


	//	LCDM
	F[12]= G[13];
	F[13]= -1./a*(3.+HA1(a,omega0)/(hub1))*G[13]+3./2.*omega0/(hub1*ap5)*G[12];


 // J and L for 3rd order ndgp kernels

 //J
 F[14] = G[15];
 F[15] = -1./a*(3.+HA1(a,omega0)/(hub1))*G[15]+3./2.*omega0/(hub1*ap5)*mu(a,1.,omega0,p1,p2,p3)*G[14]-1./(a*a*hub1*24.*pow(beta(a,omega0,p1),3)*p1) *pow(omega0/(a*a*a),2)*G[0]*G[0]*G[0]/7.;

 //K
 F[16] = G[17];
 F[17] = -1./a*(3.+HA1(a,omega0)/(hub1))*G[17]+3./2.*omega0/(hub1*ap5)*mu(a,1.,omega0,p1,p2,p3)*G[16]-1./(a*a*hub1*24.*pow(beta(a,omega0,p1),3)*p1) *pow(omega0/(a*a*a),2)*G[0]*G[0]*G[0]*5./7.;

 //F3
 F[18] = G[19] ;
 F[19] = -1./a*(3.+HA1(a,omega0)/(hub1))*G[19]+3./2.*omega0/(hub1*ap5)*mu(a,1.,omega0,p1,p2,p3)*G[18]+3./2.*omega0/(hub1*ap5)*mu(a,1.,omega0,p1,p2,p3)*(G[4]-G[0])*G[0]+G[1]*(G[5]-G[1]);


	return GSL_SUCCESS;
}

//Jacobian of the system required when calling the system evolver, the below is neither needed, nor correct for solving
int jac (double a, const double G[], double *dfdy, double dfdt[], void *params)
{
	return GSL_SUCCESS;
}

//loops for evol factors, A is the scale factor
//Normalized to D[a] = a at early times and a = 1 today.
void IOW::inite(double A, double omega0, double par1, double par2, double par3)
{

	double a = 0.00001;
	double G[20] = {a,1., a, 1., a,  1., a, 1.,  a, 1., a, 1., a, 1., a, 1., a, 1.,a,1.}; // initial conditions
	struct param_type2 params = {omega0,par1,par2,par3,1};

// Solutions of evolution factors @ a=A
	gsl_odeiv2_system sys = {func, jac, 20, &params};
	gsl_odeiv2_driver * d =
	gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
								  1e-6, 1e-6, 1e-6);

	int status1 = gsl_odeiv2_driver_apply (d, &a, A, G);


  //Hubble
  H_spt  = HA(A,omega0);
    double mod = 1.;
		D_spt  = G[0];
		F_spt  = mod*(G[4]-G[0]);
		Cx_spt = mod*(G[6]-G[0]);
		I_spt  = mod*(G[8]-G[0]);
    J_spt  = mod*(G[14]-G[0]);
    K_spt  = mod*(G[16]-G[0]);
		KL_spt = mod*(G[10]-G[0]);
    F3_spt = mod*(G[18]-G[0]);


		//Their time (not scale factor!) derivatives
		fdgp_spt =  G[1]/G[0]*A;
		Dd_spt   =  G[1]*H_spt*A;
		Fd_spt   = (G[5]-G[1])*H_spt*A;
		Cdx_spt  = (G[7]-G[1])*H_spt*A;
		Id_spt   = (G[9]-G[1])*H_spt*A;
    Jd_spt   = (G[15]-G[1])*H_spt*A;
    Kd_spt   = (G[17]-G[1])*H_spt*A;
		KLd_spt  = (G[11]-G[1])*H_spt*A;
    F3d_spt  = (G[19]-G[1])*H_spt*A;

		//LCDM
		Dl_spt = G[12];
		fl_spt = G[13]/G[12]*A;

		gsl_odeiv2_driver_free(d);

// Solution for LCDM linear growth @ a=1 for our normalisation of the linear power spectrum
	gsl_odeiv2_driver * d1 =
	gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
								   1e-6, 1e-6, 1e-6);

	int status2 = gsl_odeiv2_driver_apply (d1, &a, 1. , G);

	dnorm_spt = G[12];

	gsl_odeiv2_driver_free(d1);


}

// //loops for evol factors, A is the scale factor
// //Normalized to D[a] = a at early times and a = 1 today.
// void IOW::inite(double A, double omega0, double par1, double par2, double par3)
// {
//
// 	double a = 0.00001;
// 	double G[20] = {a,1., a, 1., a,  1., a, 1.,  a, 1., a, 1., a, 1., a, 1., a, 1.,a,1.}; // initial conditions
// 	struct param_type2 params = {omega0,par1,par2,par3,1};
//
// // Solutions of evolution factors @ a=A
// 	gsl_odeiv2_system sys = {func, jac, 20, &params};
// 	gsl_odeiv2_driver * d =
// 	gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
// 								  1e-6, 1e-6, 1e-6);
//
// 	int status1 = gsl_odeiv2_driver_apply (d, &a, A, G);
//
//
//   //Hubble
//   H_spt  = HA(A,omega0);
//     double mod = 1.;
// 		D_spt  = G[0];
// 		F_spt  = mod*(G[4]-G[0]);
// 		Cx_spt = mod*(G[6]-G[0]);
// 		I_spt  = mod*(G[8]-G[0]);
//     J_spt  = mod*(G[14]-G[0]);
//     K_spt  = mod*(G[16]-G[0]);
// 		KL_spt = mod*(G[10]-G[0]);
//     F3_spt = mod*(G[18]-G[0]);
//
//
// 		//Their time (not scale factor!) derivatives
// 		fdgp_spt =  G[1]/G[0]*A;
// 		Dd_spt   =  G[1]*H_spt*A;
// 		Fd_spt   = (G[5]-G[1])*H_spt*A;
// 		Cdx_spt  = (G[7]-G[1])*H_spt*A;
// 		Id_spt   = (G[9]-G[1])*H_spt*A;
//     Jd_spt   = (G[15]-G[1])*H_spt*A;
//     Kd_spt   = (G[17]-G[1])*H_spt*A;
// 		KLd_spt  = (G[11]-G[1])*H_spt*A;
//     F3d_spt  = (G[19]-G[1])*H_spt*A;
//
// 		//LCDM
// 		Dl_spt = G[12];
// 		fl_spt = G[13]/G[12]*A;
//
// 		gsl_odeiv2_driver_free(d);
//
// // Solution for LCDM linear growth @ a=1 for our normalisation of the linear power spectrum
// 	gsl_odeiv2_driver * d1 =
// 	gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
// 								   1e-6, 1e-6, 1e-6);
//
// 	int status2 = gsl_odeiv2_driver_apply (d1, &a, 1. , G);
//
// 	dnorm_spt = G[12];
//
// 	gsl_odeiv2_driver_free(d1);
//
//
// }

//// Larger set for F4 and G4 dgp kernel evolution factors

int funce2(double a, const double G[], double F[], void *params){
	param_type2 p = *(param_type2 *)(params);
	real omega0 = p.omega00;
	real p1 = p.par1;
	real p2 = p.par2;
	real p3 = p.par3;

  double hub1 = pow2(HA(a,omega0));
  double hub2 = HA1(a,omega0)/hub1;
  double ap5 = pow(a,5);
  double mua = mu(a,1.,omega0,p1,p2,p3);
  double betacub = pow(beta(a,omega0,p1),3);
  double omsqr = pow(omega0/(a*a*a),2);

	// D
	F[0] = G[1];
	F[1] = -1./a*(3.+hub2)*G[1]+3./2.*omega0/(hub1*ap5)*mua*G[0];

	//E
	F[2] = G[3];
	F[3] = -1./a*(3.+hub2)*G[3]+3./2.*omega0/(hub1*ap5)*mua*G[2]+G[1]*G[1]+3./2.*omega0/(hub1*ap5)*mua*G[0]*G[0];

	//F
	F[4] = G[5];
	F[5] = -1./a*(3.+hub2)*G[5]+3./2.*omega0/(hub1*ap5)*mua*G[4]-1./(a*a*hub1*24.*betacub*p1) *omsqr*G[0]*G[0];

	//C
	F[6] =G[7];
	F[7] = -1./a*(3.+hub2)*G[7]+3./2.*omega0/(hub1*ap5)*mua*G[6]+G[1]*(G[5]-G[1]);

	//I
	F[8] = G[9] ;
	F[9] = -1./a*(3.+hub2)*G[9]+3./2.*omega0/(hub1*ap5)*mua*G[8] + 3./2.*omega0/(hub1*ap5)*mua*(G[4]-G[0])*G[0] - 1./(a*a*hub1*24.*betacub*p1)*omsqr*G[0]*G[0]*G[0] + G[1]*(G[5]-G[1]);

	//KL
	F[10] = G[11] ;
	F[11] = -1./a*(3.+hub2)*G[11]+3./2.*omega0/(hub1*ap5)*mua*G[10]
	-1./(a*a*hub1*12.*betacub*p1)*omsqr*G[0]*(G[2]+G[4]-2.*G[0])
	+1./(a*a*hub1*144.*pow(beta(a,omega0,p1),5)*pow(p1,2))*pow(omega0/(a*a*a),3)*G[0]*G[0]*G[0];


	//	LCDM
	F[12]= G[13];
	F[13]= -1./a*(3.+hub2)*G[13]+3./2.*omega0/(hub1*ap5)*G[12];


 // J,F3 K for 3rd order ndgp kernels

 //J
 F[14] = G[15];
 F[15] = -1./a*(3.+hub2)*G[15]+3./2.*omega0/(hub1*ap5)*mua*G[14] - 1./(a*a*hub1*24.*betacub*p1) *omsqr*G[0]*G[0]*G[0]/7.;

 //K
 F[16] = G[17];
 F[17] = -1./a*(3.+hub2)*G[17]+3./2.*omega0/(hub1*ap5)*mua*G[16] - 1./(a*a*hub1*24.*betacub*p1) *omsqr*G[0]*G[0]*G[0]*5./7.;

 //F3
 F[18] = G[19] ;
 F[19] = -1./a*(3.+hub2)*G[19]+3./2.*omega0/(hub1*ap5)*mua*G[18] + 3./2.*omega0/(hub1*ap5)*mua*(G[4]-G[0])*G[0] + G[1]*(G[5]-G[1]);


// 4th order equations A-O (15)

// B

F[20] =G[21];
F[21] = -1./a*(3.+hub2)*G[21]+3./2.*omega0/(hub1*ap5)*mua*G[20]+2.*G[1]*G[0]*(G[5]-G[1]);

// C

F[22] =G[23];
F[23] = -1./a*(3.+hub2)*G[23]+3./2.*omega0/(hub1*ap5)*mua*G[22]+(G[5]-G[1])*(G[5]-G[1]);

// F
F[24] =G[25];
F[25] = -1./a*(3.+hub2)*G[25]+3./2.*omega0/(hub1*ap5)*mua*G[24]
              + 3./2.*omega0/(hub1*ap5)*mua*G[0]*G[0]*(G[4]-G[0]) + G[1]*G[1]*(G[4]-G[0]) + G[1]*G[0]*(G[5]-G[1]);

// G
F[26] =G[27];
F[27] = -1./a*(3.+hub2)*G[27]+3./2.*omega0/(hub1*ap5)*mua*G[26]
              + 3./2.*omega0/(hub1*ap5)*mua*G[0]*G[0]*(G[4]-G[0]) - 1./(a*a*hub1*24.*betacub*p1)*omsqr*G[0]*G[0]*G[0]*G[0] + 2.*G[1]*G[0]*(G[5]-G[1]);

// H
F[28] =G[29];
F[29] = -1./a*(3.+hub2)*G[29]+3./2.*omega0/(hub1*ap5)*mua*G[28]
              + 3./2.*omega0/(hub1*ap5)*mua*(G[4]-G[0])*(G[4]-G[0]) - 1./(a*a*hub1*24.*betacub*p1)*omsqr*G[0]*G[0]*(G[4]-G[0]) + (G[5]-G[1])*(G[5]-G[1]);


//A
//Ai
F[30] =G[31];
F[31] = -1./a*(3.+hub2)*G[31]+3./2.*omega0/(hub1*ap5)*mua*G[30]+G[1]*(G[7]-G[1]); //C'
//Aii
F[32] =G[33];
F[33] = -1./a*(3.+hub2)*G[33]+3./2.*omega0/(hub1*ap5)*mua*G[32]+G[1]*((G[19]-G[1]) -G[1]*(G[4]-G[0])); //F3'-D1'F2
//Aiii
F[34] =G[35];
F[35] = -1./a*(3.+hub2)*G[35]+3./2.*omega0/(hub1*ap5)*mua*G[34]+G[1]*((G[9]-G[1])-G[0]*(G[5]-G[1])); // I3'-DF2'
//Aiv
F[36] =G[37];
F[37] = -1./a*(3.+hub2)*G[37]+3./2.*omega0/(hub1*ap5)*mua*G[36]+G[1]*(G[15]-G[1]);// J3'
//Av
F[38] =G[39];
F[39] = -1./a*(3.+hub2)*G[39]+3./2.*omega0/(hub1*ap5)*mua*G[38]+G[1]*(G[11]-G[17]); //L3'
//Avi
F[40] =G[41];
F[41] = -1./a*(3.+hub2)*G[41]+3./2.*omega0/(hub1*ap5)*mua*G[40]+G[1]*(G[17]-G[1]); //K3'

//D

//Di
F[42] =G[43];
F[43] = -1./a*(3.+hub2)*G[43]+3./2.*omega0/(hub1*ap5)*mua*G[42]
              + 3./2.*omega0/(hub1*ap5)*mua*G[0]*(G[6]-G[0]) +G[1]*(G[7]-G[1]); //C3
//Dii
F[44] =G[45];
F[45] = -1./a*(3.+hub2)*G[45]+3./2.*omega0/(hub1*ap5)*mua*G[44]
              + 3./2.*omega0/(hub1*ap5)*mua*G[0]*(G[18]-G[0]) +G[1]*(G[19]-G[1]); //F3
//Diii
F[46] =G[47];
F[47] = -1./a*(3.+hub2)*G[47]+3./2.*omega0/(hub1*ap5)*mua*G[46]
              + 3./2.*omega0/(hub1*ap5)*mua*G[0]*(G[8]-G[0]) + G[1]*(G[9]-G[1]); // I3


//Div
F[48] =G[49];
F[49] = -1./a*(3.+hub2)*G[49]+3./2.*omega0/(hub1*ap5)*mua*G[48]
              + 3./2.*omega0/(hub1*ap5)*mua*G[0]*(G[14]-G[0]) +G[1]*(G[15]-G[1]); // J3
//Dv
F[50] =G[51];
F[51] = -1./a*(3.+hub2)*G[51]+3./2.*omega0/(hub1*ap5)*mua*G[50]
              + 3./2.*omega0/(hub1*ap5)*mua*G[0]*(G[10]-G[16]) +G[1]*(G[11]-G[17]); // L3
//Dvi
F[52] =G[53];
F[53] = -1./a*(3.+hub2)*G[53]+3./2.*omega0/(hub1*ap5)*mua*G[52]
              + 3./2.*omega0/(hub1*ap5)*mua*G[0]*(G[16]-G[0]) +G[1]*(G[17]-G[1]); //K3


// E
//Ei
F[54] =G[55];
F[55] = -1./a*(3.+hub2)*G[55]+3./2.*omega0/(hub1*ap5)*mua*G[54]
              + 3./2.*omega0/(hub1*ap5)*mua*G[0]*(G[6]-G[0])  + G[1]*(G[5]-G[1])*G[0] + G[1]*(G[7]-G[1]); //C3

//Eii
F[56] =G[57];
F[57] = -1./a*(3.+hub2)*G[57]+3./2.*omega0/(hub1*ap5)*mua*G[56]
              + 3./2.*omega0/(hub1*ap5)*mua*G[0]*(G[18]-G[0]) + 3./2.*omega0/(hub1*ap5)*mua*(G[4]-G[0])*G[0]*G[0] + G[1]*(G[5]-G[1])*G[0]
              - 3./2.*omega0/(hub1*ap5)*mua*G[0]*G[0]*(G[4]-G[0])
              + G[1]*(G[19]-G[1]) - G[1]*G[1]*(G[4]-G[0]) - G[1]*G[0]*(G[5]-G[1]); //F3
//Eiii
F[58] =G[59];
F[59] = -1./a*(3.+hub2)*G[59]+3./2.*omega0/(hub1*ap5)*mua*G[58]
              + 3./2.*omega0/(hub1*ap5)*mua*G[0]*(G[8]-G[0]) + 3./2.*omega0/(hub1*ap5)*mua*(G[4]-G[0])*G[0]*G[0] + G[1]*(G[5]-G[1])*G[0]
              -3./2.*omega0/(hub1*ap5)*mua*G[0]*G[0]*(G[4]-G[0])
              + G[1]*(G[9]-G[1]) - 2.*G[1]*G[0]*(G[5]-G[1]); // I3


//Eiv
F[60] =G[61];
F[61] = -1./a*(3.+hub2)*G[61]+3./2.*omega0/(hub1*ap5)*mua*G[60]
              + 3./2.*omega0/(hub1*ap5)*mua*G[0]*(G[14]-G[0])  - 1./(a*a*hub1*24.*betacub*p1)*omsqr*G[0]*G[0]*G[0]*G[0]/7. + G[1]*(G[15]-G[1]); // J3


//Ev (L term)
F[62] =G[63];
F[63] = -1./a*(3.+hub2)*G[63]+3./2.*omega0/(hub1*ap5)*mua*G[62]
            +3./2.*omega0/(hub1*ap5)*mua*G[0]*(G[10]-G[16])
            -1./(a*a*hub1*12.*betacub*p1)*omsqr*G[0]*G[0]*(G[2]+G[4]-2.*G[0])
          	+1./(a*a*hub1*144.*pow(beta(a,omega0,p1),5)*p1*p1)*pow(omega0/(a*a*a),3)*G[0]*G[0]*G[0]*G[0]
            + 1./(a*a*hub1*24.*betacub*p1) *omsqr*G[0]*G[0]*G[0]*G[0]*5./7.
            + G[1]*(G[11]-G[17]);  // L3

            //          + 3./2.*omega0/(hub1*ap5)*mua*G[0]*(G[10]-G[16])
            //          -1./(a*a*hub1*24.*betacub*p1)*omsqr*G[0]*G[0]*(G[4]-G[0])
            //          +1./(a*a*hub1*144.*pow(beta(a,omega0,p1),5)*pow(p1,2))*pow(omega0/(a*a*a),3)*G[0]*G[0]*G[0]*G[0]

//Evi (K term)
F[64] =G[65];
F[65] = -1./a*(3.+hub2)*G[65]+3./2.*omega0/(hub1*ap5)*mua*G[64]
              + 3./2.*omega0/(hub1*ap5)*mua*G[0]*(G[16]-G[0]) - 1./(a*a*hub1*24.*betacub*p1)*omsqr*G[0]*G[0]*G[0]*G[0]*5./7. + G[1]*(G[17]-G[1]); //K3




//Ii
F[66] =G[67];
F[67] = -1./a*(3.+hub2)*G[67]+3./2.*omega0/(hub1*ap5)*mua*G[66]
              -1./(a*a*hub1*24.*betacub*p1) *omsqr*G[0]*(G[6]-G[0]);

//Iii
F[68] =G[69];
F[69] = -1./a*(3.+hub2)*G[69]+3./2.*omega0/(hub1*ap5)*mua*G[68]
            -1./(a*a*hub1*24.*betacub*p1) *omsqr*G[0]*(G[18]-G[0]);

//Iiii
F[70] =G[71];
F[71] = -1./a*(3.+hub2)*G[71]+3./2.*omega0/(hub1*ap5)*mua*G[70]
            -1./(a*a*hub1*24.*betacub*p1) *omsqr*G[0]*(G[8]-G[0]);

//Iiv
F[72] =G[73];
F[73] = -1./a*(3.+hub2)*G[73]+3./2.*omega0/(hub1*ap5)*mua*G[72]
            -1./(a*a*hub1*24.*betacub*p1) *omsqr*G[0]*(G[14]-G[0]);

//Iv // L
F[74] =G[75];
F[75] = -1./a*(3.+hub2)*G[75]+3./2.*omega0/(hub1*ap5)*mua*G[74]
            -1./(a*a*hub1*24.*betacub*p1) *omsqr*G[0]*(G[10]-G[16]);

//Ivi // K
F[76] =G[77];
F[77] = -1./a*(3.+hub2)*G[77]+3./2.*omega0/(hub1*ap5)*mua*G[76]
            -1./(a*a*hub1*24.*betacub*p1) *omsqr*G[0]*(G[16]-G[0]);





// J
F[78] =G[79];
F[79] = -1./a*(3.+hub2)*G[79]+3./2.*omega0/(hub1*ap5)*mua*G[78]
            -2./(a*a*hub1*24.*betacub*p1) *omsqr*G[0]*G[0]*(G[4]-G[0])
            +1./(a*a*hub1*144.*pow(beta(a,omega0,p1),5)*pow(p1,2))*pow(omega0/(a*a*a),3)*G[0]*G[0]*G[0]*G[0];

// K
F[80] =G[81];
F[81] = -1./a*(3.+hub2)*G[81]+3./2.*omega0/(hub1*ap5)*mua*G[80]
            -1./(a*a*hub1*24.*betacub*p1) *omsqr*(G[4]-G[0])*(G[4]-G[0])
            +1./(a*a*hub1*144.*pow(beta(a,omega0,p1),5)*pow(p1,2))*pow(omega0/(a*a*a),3)*G[0]*G[0]*(G[4]-G[0])
            -2./(a*a*hub1*1728.*pow(beta(a,omega0,p1),7)*pow(p1,3))*pow2(omsqr)*G[0]*G[0]*G[0]*G[0];

//L AND M
F[82] =G[83];
F[83] = -1./a*(3.+hub2)*G[83]+3./2.*omega0/(hub1*ap5)*mua*G[82]
            -1./(a*a*hub1*24.*betacub*p1)*omsqr*G[0]*G[0]*G[0]*G[0];

// N
F[84] =G[85];
F[85] = -1./a*(3.+hub2)*G[85]+3./2.*omega0/(hub1*ap5)*mua*G[84]
            +2./(a*a*hub1*144.*pow(beta(a,omega0,p1),5)*pow(p1,2))*pow(omega0/(a*a*a),3)*G[0]*G[0]*(G[4]-G[0])
            -1./(a*a*hub1*2.*1728.*pow(beta(a,omega0,p1),7)*pow(p1,3))*pow(omega0/(a*a*a),4)*G[0]*G[0]*G[0]*G[0];

// O
F[86] =G[87];
F[87] = -1./a*(3.+hub2)*G[87]+3./2.*omega0/(hub1*ap5)*mua*G[86]
            +2./(a*a*hub1*144.*pow(beta(a,omega0,p1),5)*pow(p1,2))*pow(omega0/(a*a*a),3)*G[0]*G[0]*G[0]*G[0];

return GSL_SUCCESS;
}


//loops for evol factors, A is the scale factor
//Normalized to D[a] = a at early times and a = 1 today.
void IOW::inite2(double A, double omega0, double par1, double par2, double par3)
{

	double a = 0.00001;
	double G[88] = {a,1., a, 1., a,  1., a, 1.,  a, 1., a, 1., a, 1., a, 1., a, 1.,a,1.,
                  a,1., a, 1., a,  1., a, 1.,  a, 1., a, 1., a, 1., a, 1., a, 1.,a,1.,
                  a,1., a, 1., a,  1., a, 1.,  a, 1., a, 1., a, 1., a, 1., a, 1.,a,1.,
                  a,1., a, 1., a,  1., a, 1.,  a, 1., a, 1., a, 1., a, 1., a, 1.,a,1.,a, 1., a, 1., a, 1.,a,1.}; // initial conditions

  struct param_type2 params = {omega0,par1,par2,par3,1};

// Solutions of evolution factors @ a=A
	gsl_odeiv2_system sys = {funce2, jac, 88, &params};
	gsl_odeiv2_driver * d =
	gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
								  1e-3, 1e-3, 1e-1);

	int status1 = gsl_odeiv2_driver_apply (d, &a, A, G);

  //Hubble
  H_spt  = HA(A,omega0);

// 1st, 2nd and 3rd order growth functions

    double mod = 1.;
		D_spt  = G[0];
		F_spt  =mod*(G[4]-G[0]);
		Cx_spt =mod*(G[6]-G[0]);
		I_spt  =mod*(G[8]-G[0]);
    J_spt  =mod*(G[14]-G[0]);
    K_spt  =mod*(G[16]-G[0]);
		KL_spt =mod*(G[10]-G[0]);
    F3_spt =mod*(G[18]-G[0]);


		//Their time (not scale factor!) derivatives
		fdgp_spt =  G[1]/G[0]*A;
		Dd_spt   =  G[1]*H_spt*A;
		Fd_spt   = mod*(G[5]-G[1])*H_spt*A;
		Cdx_spt  = mod*(G[7]-G[1])*H_spt*A;
		Id_spt   = mod*(G[9]-G[1])*H_spt*A;
    Jd_spt   = mod*(G[15]-G[1])*H_spt*A;
    Kd_spt   = mod*(G[17]-G[1])*H_spt*A;
		KLd_spt  = mod*(G[11]-G[1])*H_spt*A;
    F3d_spt  = mod*(G[19]-G[1])*H_spt*A;

		//LCDM
		Dl_spt = G[0];//G[12];
		fl_spt = G[1]/G[0]*A;//G[13]/G[12]*A;

// 4th order evolution factors
    for(int myint=0; myint<34; myint++){
          evol4[myint] = mod*(G[20+2*myint]-G[0]);
        }
    for(int myint=0; myint<34; myint++){
          devol4[myint] = mod*(G[21+2*myint]-G[1])*H_spt*A;
        }

// add in other factors for G4

// D4
devol4[11] += -Dd_spt/H_spt*Cx_spt;
devol4[12] += -Dd_spt/H_spt*F3_spt;
devol4[13] += -Dd_spt/H_spt*I_spt;
devol4[14] += -Dd_spt/H_spt*J_spt;
devol4[15] += -Dd_spt/H_spt*(KL_spt-K_spt);
devol4[16] += -Dd_spt/H_spt*K_spt;
//E4
devol4[17] += -D_spt*Cdx_spt/H_spt;
devol4[18] += -D_spt*(F3d_spt/H_spt-D_spt*fdgp_spt*F_spt);
devol4[19] += -D_spt*(Id_spt/H_spt-D_spt*Fd_spt/H_spt);
devol4[20] += -D_spt*Jd_spt/H_spt;
devol4[21] += -D_spt*(KLd_spt-Kd_spt)/H_spt;
devol4[22] += -D_spt*Kd_spt/H_spt;
// F4
devol4[2] += -F_spt*Dd_spt*D_spt/H_spt;
//G4
devol4[3] += -Fd_spt*D_spt*D_spt/H_spt;
//G4
devol4[4] += -Fd_spt*F_spt/H_spt;

		gsl_odeiv2_driver_free(d);

// Solution for LCDM linear growth @ a=1 for our normalisation of the linear power spectrum
	gsl_odeiv2_driver * d1 =
	gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
								   1e-6, 1e-6, 1e-6);

	int status2 = gsl_odeiv2_driver_apply (d1, &a, 1. , G);

	dnorm_spt = G[12];

	gsl_odeiv2_driver_free(d1);

}




///////// NON-SEPARABLE //////////////

/* Euler and Continuity equations for numerical kernels */

//HA2 = -dH/dt/H^2
double HA2(double a, double omega0){
	return 3.*omega0/(2.*HA(a,omega0)*HA(a,omega0)*pow(a,3))
	;}

/* Parameters passed to system of Euler and continuity equations*/
// k (magnitude) and x (angular) values for the system of equations
// k1.k2=k1k2x2 , k1.k3 = k1k3x3, k2.k3=k2k3x1
// EDIT : Add in new gravity parameters to this list as required (par1 is the nDGP parameter omega_rc as default)
struct param_type3 {
  real kk1;
  real xx1;
	real kk2;
  real xx2;
	real kk3;
  real xx3;
	real omega00;
	real par1;
	real par2;
	real par3;
  real arg1;
  real arg2;
  real arg3;
  real arg4;

};


// G[i] are used to construct symmetrized kernels used in G3,F3 equations as well as RSD A-term

int funcn1(double a, const double G[], double F[], void *params)
{
	param_type3 p = *(param_type3 *)(params);
	real k1 = p.kk1;
	real x1= p.xx1;
	real k2 = p.kk2;
	real x2= p.xx2;
	real k3 = p.kk3;
	real x3= p.xx3;
	real omega0 = p.omega00;
	real p1 = p.par1;
	real p2 = p.par2;
	real p3 = p.par3;
  real karg1 = p.arg1;
  real karg2 = p.arg2;
  real karg3 = p.arg3;
  real karg4 = p.arg4;


  double a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14;
  double b1,b2,b3,b4,b5,b6,b7;
  double hub = HA2(a,omega0);

  a1 = alpha(k2,karg1,(k1*x2-k2)/karg1);
  a2 = alpha(karg1,k2,(k1*x2-k2)/karg1);
  a3 = alpha(k1,karg1,(k2*x2-k1)/karg1);
  a4 = alpha(karg1,k1,(k2*x2-k1)/karg1);
  a5 = alpha(k2,k1,x2);
  a6 = alpha(k1,k2,x2);
  a7 = alpha(k3,k1,x3);
  a8 = alpha(k1,k3,x3);

  a9 = alpha(k1,karg3,(k2*x2+k3*x3)/karg3);
  a10 = alpha(k3,karg4,(k1*x3+k2*x1)/karg4);
  a11 = alpha(k2,karg2,(k3*x1+k1*x2)/karg2);
  a12 = alpha(karg3,k1,(k2*x2+k3*x3)/karg3);
  a13 = alpha(karg4,k3,(k1*x3+k2*x1)/karg4);
  a14 = alpha(karg2,k2,(k3*x1+k1*x2)/karg2);

  b1 = beta1(k2,karg1,(k1*x2-k2)/karg1);
  b2 = beta1(k1,karg1,(k2*x2-k1)/karg1);
  b3 = beta1(k1,k2,x2);
  b4 = beta1(k1,k3,x3);
  b5 = beta1(k1,karg3,(k2*x2+k3*x3)/karg3);
  b6 = beta1(k3,karg4,(k1*x3+k2*x1)/karg4);
  b7 = beta1(k2,karg2,(k3*x1+k1*x2)/karg2);





	/* 1st order */

	//1. F1/G1(k)
	F[0] = -G[1]/a;
	F[1] =1./a*(-(2.-hub)*G[1]-hub*G[0]*mu(a,k1,omega0,p1,p2,p3));


	//2. F1/G1(k-p)
	F[2] = -G[3]/a;
	F[3] =1./a*(-(2.-hub)*G[3]-hub*G[2]*mu(a,karg1,omega0,p1,p2,p3));

	//3. F1/G1(p)
	F[4] = -G[5]/a;
	F[5] =1./a*(-(2.-hub)*G[5]-hub*G[4]*mu(a,k2,omega0,p1,p2,p3));

	/* 2nd order */

	//4. F2/G2(p,k-p) (P22)
	F[6] =1./a*(-(a1*G[5]*G[2]+a2*G[3]*G[4])/2.-G[7]);
	F[7] =1./a*(-(2.-hub)*G[7]-hub*G[6]*mu(a,k1,omega0,p1,p2,p3) - gamma2(a, omega0, k1, k2,karg1,(k1*x2-k2)/karg1,p1,p2,p3)*G[4]*G[2] - b1*G[5]*G[3]);


	//5. F2/G2(-k,k-p)
	F[8] =1./a*(-(a3*G[1]*G[2]+a4*G[3]*G[0])/2.-G[9]);
	F[9] =1./a*(-(2.-hub)*G[9]-hub*G[8]*mu(a,k2,omega0,p1,p2,p3) - gamma2(a, omega0, k2, k1, karg1,(k2*x2-k1)/karg1,p1,p2,p3)*G[2]*G[0] - b2*G[3]*G[1]);



	//6. F2/G2(p,-p)
	F[10] =1./a*(-G[11]) ; // set 1st term to 0 because x1 is not identically -1 (needed to keep evolution at 0 for GR)
	F[11] =1./a*(-(2.-hub)*G[11]-hub*G[10]*mu(a,0.,omega0,p1,p2,p3) - gamma2(a, omega0, 0., k2,k3,x1,p1,p2,p3)*G[4]*G[4]);


	//7. F2/G2(p,k)
	F[12] =1./a*(-(a5*G[5]*G[0]+a6*G[1]*G[4])/2.-G[13]) ;
	F[13] =1./a*(-(2.-hub)*G[13]-hub*G[12]*mu(a,karg4,omega0,p1,p2,p3) - gamma2(a, omega0, karg4, k2,k1,x2,p1,p2,p3)*G[4]*G[0]-b3*G[5]*G[1]);



	//8. F2/G2(-p,k)=F2/G2(p,-k)
	F[14] =1./a*(-(a7*G[5]*G[0]+a8*G[1]*G[4])/2.-G[15]) ;
	F[15] =1./a*(-(2.-hub)*G[15]-hub*G[14]*mu(a,karg1,omega0,p1,p2,p3) -gamma2(a, omega0, karg1, k3,k1,x3,p1,p2,p3)*G[4]*G[0]-b4*G[5]*G[1]);


	//9. 3rd order  ~ F3/G3(k,p,-p)
	F[16] = - 1./(3.*a)*(a9*G[10]*G[1]
						+ a10*G[12]*G[5]
						+ a11*G[14]*G[5]

						+ a13*G[13]*G[4]
						+ a14*G[15]*G[4]
						+ a12*G[11]*G[0]

						+3.*G[17]) ;


	F[17] =1./(3.*a)*(-3.*(2.-hub)*G[17]-3.*hub*G[16]*mu(a,k1,omega0,p1,p2,p3)

					 -2.*b5*G[1]*G[11]

					 -2.*b6*G[5]*G[13]

					 -2.*b7*G[5]*G[15]

					 -2.*gamma2(a, omega0, k1, k1,karg3,(k2*x2+k3*x3)/karg3,p1,p2,p3)*G[0]*G[10]
					 -2.*gamma2(a, omega0, k1, k3,karg4,(k1*x3+k2*x1)/karg4,p1,p2,p3)*G[4]*G[12]
					 -2.*gamma2(a, omega0, k1, k2,karg2,(k3*x1+k1*x2)/karg2,p1,p2,p3)*G[4]*G[14]

					   -(gamma3(a, omega0, k1, k1,k2,k3,x1,x2,x3,p1,p2,p3)
					   + gamma3(a, omega0, k1, k3,k1,k2,x2,x3,x1,p1,p2,p3)
					   + gamma3(a, omega0, k1, k2,k3,k1,x3,x1,x2,p1,p2,p3))*G[4]*G[4]*G[0]);




	return GSL_SUCCESS;
}

//Loops to construct numerical kernels
// Kmax and Kmin indicates the maximum and minimum values of K
// k loop dictates the angular parameter x =k1.k/(k^2*r)
// j loop dictates the magnitude of k1 (k1=rk)

// Default initialization at only scale factor A
// k is your k vector dependence (See mu(a,k))
// YMAX = QMAXp/kmin and YMIN=QMINp/kmax where kmin and kmax are the k-limits and QMINp and QMAXp are k1 limits (integrated over)
// set k = 1 and YMIN = QMINp/kmax, YMAX = QMAXp/kmin for no k -dependence (initialize once)
// set kmin = kmax = k for k dependence (initialize for each k)

// EDIT : Add in new gravity parameters to function's input as required (par1 is the nDGP parameter omega_rc as default)
// EDIT : if more than 1 gravity parameter is needed, one must also edit initn function in SPT.h file.

void IOW::initn(double A, double YMIN, double YMAX, double k, double omega0, double par1, double par2, double par3)
{
//#pragma omp parallel for schedule(dynamic)

for (int j = 0; j< n2; ++j)
{
  /* k - magnitude sampling */
    double R = (YMIN) * exp(j*log(YMAX/YMIN)/(n2*1.-1.)); //exponential sampling
  //  	double R = (j*1.)/(n2*1.-1)*(YMAX-YMIN) + YMIN; // linear sampling
  //	double R = pow2(j*1./(n2*1.-1))*YMAX+YMIN; //quadratic sampling

  /*Parameter selection for kernels */

  double k1 = k;
  double k2 = k*R;
  double  k3 = k2;

  double x1 = XMIN;
  double karg3 = sqrt(k2*k2+2.*k2*k3*x1+k3*k3);

	for ( int i = 0; i < n1; ++i)
	{

      double x2 = x128[i]; // x128 goes from - 0.99995 to 0.99995
			double x3 = -x2;

      double  karg1 = sqrt(k2*k2+k1*k1-2.*k2*k1*x2);
      double  karg2 = sqrt(k3*k3+2.*k3*k1*x3+k1*k1);
      double  karg4 = sqrt(k1*k1+2.*k1*k2*x2+k2*k2);


			// Initial scale factor for solving system of equations
			//Sets 1-3
				double a = 0.0001;

				// Einstein de Sitter initial condtions for 2nd order kernels

		//		4. F2/G2(p,k-p) (P22)
				// double EdsF2=a*a*(10./14.*alphas(k2,karg1,(k1*x2-k2)/karg1)+ 2./7.*beta1(k2,karg1,(k1*x2-k2)/karg1));
				// double EdsG2=-a*a*(6./14.*alphas(k2,karg1,(k1*x2-k2)/karg1)+ 4./7.*beta1(k2,karg1,(k1*x2-k2)/karg1));
        //
				// //5. F2/G2(-k,k-p)
				// double EdsF2C=a*a*(10./14.*alphas(k1,karg1,(k2*x2-k1)/karg1)+ 2./7.*beta1(k1,karg1,(k2*x2-k1)/karg1));
				// double EdsG2C=-a*a*(6./14.*alphas(k1,karg1,(k2*x2-k1)/karg1) + 4./7.*beta1(k1,karg1,(k2*x2-k1)/karg1));
        //
				// //6. F2/G2(p,-p)
				// double EdsF2A= 0.;
				// double EdsG2A= 0.; //Setting manually to 0 because xmin is not identically -1
        //
				// //7. F2/G2(p,k)
				// double EdsF2B= a*a*(10./14.*alphas(k2,k1,x2)+ 2./7.*beta1(k2,k1,x2));
				// double EdsG2B= -a*a*(6./14.*alphas(k2,k1,x2) + 4./7.*beta1(k2,k1,x2));
				// //8. F2/G2(-p,k)=F2/G2(p,-k)
				// double Edsf2d= a*a*(10./14.*alphas(k3,k1,x3)+ 2./7.*beta1(k3,k1,x3));
				// double Edsg2d= -a*a*(6./14.*alphas(k3,k1,x3)+ 4./7.*beta1(k3,k1,x3));
        //
				// // Einstein de Sitter initial condtions for 3rd order kernels
				// //9. 3rd order  ~ F3/G3(k,p,-p)
				// double Edsf3 = pow3(a)/3.*(2./63.*beta1(k1,karg3,(k2*x2+k3*x3)/karg3)*(beta1(k2,k3,x1)+6./4.*alphas(k2,k3,x1))
				// 							+1./18.*alpha(k1,karg3,(k2*x2+k3*x3)/karg3)*(beta1(k2,k3,x1)+10./2.*alphas(k2,k3,x1))
				// 						    +1./9.*alpha(karg3,k1,(k2*x2+k3*x3)/karg3)*(beta1(k2,k3,x1)+6./4.*alphas(k2,k3,x1))
        //
				// 						   +2./63.*beta1(k3,karg4,(k1*x3+k2*x1)/karg4)*(beta1(k1,k2,x2)+6./4.*alphas(k1,k2,x2))
				// 						   +1./18.*alpha(k3,karg4,(k1*x3+k2*x1)/karg4)*(beta1(k1,k2,x2)+10./2.*alphas(k1,k2,x2))
				// 						   +1./9.*alpha(karg4,k3,(k1*x3+k2*x1)/karg4)*(beta1(k1,k2,x2)+6./4.*alphas(k1,k2,x2))
        //
				// 						   +2./63.*beta1(k2,karg2,(k1*x2+k3*x1)/karg2)*(beta1(k1,k3,x3)+6./4.*alphas(k1,k3,x3))
				// 						   +1./18.*alpha(k2,karg2,(k1*x2+k3*x1)/karg2)*(beta1(k1,k3,x3)+10./2.*alphas(k1,k3,x3))
				// 						   +1./9.*alpha(karg2,k2,(k1*x2+k3*x1)/karg2)*(beta1(k1,k3,x3)+6./4.*alphas(k1,k3,x3)));
        //
				// double Edsg3 = -pow3(a)/3.*(2./21.*beta1(k1,karg3,(k2*x2+k3*x3)/karg3)*(beta1(k2,k3,x1)+6./4.*alphas(k2,k3,x1))
				// 					   +1./42.*alpha(k1,karg3,(k2*x2+k3*x3)/karg3)*(beta1(k2,k3,x1)+10./2.*alphas(k2,k3,x1))
				// 					   +1./21.*alpha(karg3,k1,(k2*x2+k3*x3)/karg3)*(beta1(k2,k3,x1)+6./4.*alphas(k2,k3,x1))
        //
				// 					   +2./21.*beta1(k3,karg4,(k1*x3+k2*x1)/karg4)*(beta1(k1,k2,x2)+6./4.*alphas(k1,k2,x2))
				// 					   +1./42.*alpha(k3,karg4,(k1*x3+k2*x1)/karg4)*(beta1(k1,k2,x2)+10./2.*alphas(k1,k2,x2))
				// 					   +1./21.*alpha(karg4,k3,(k1*x3+k2*x1)/karg4)*(beta1(k1,k2,x2)+6./4.*alphas(k1,k2,x2))
        //
				// 					   +2./21.*beta1(k2,karg2,(k1*x2+k3*x1)/karg2)*(beta1(k1,k3,x3)+6./4.*alphas(k1,k3,x3))
				// 					   +1./42.*alpha(k2,karg2,(k1*x2+k3*x1)/karg2)*(beta1(k1,k3,x3)+10./2.*alphas(k1,k3,x3))
				// 					   +1./21.*alpha(karg2,k2,(k1*x2+k3*x1)/karg2)*(beta1(k1,k3,x3)+6./4.*alphas(k1,k3,x3)));


		//	  double G[18] = { a,-a,a,-a,a,-a,  EdsF2 ,EdsG2,  EdsF2C, EdsG2C, EdsF2A, EdsG2A , EdsF2B,  EdsG2B, Edsf2d, Edsg2d, Edsf3 ,Edsg3}; // initial conditions
// Non-Eds ICs
			  double G[18] = { a,-a,a,-a,a,-a, 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

			/*Parameters passed to system of equations */
// EDIT : If more than one gravity parameter is used, add them after p1
				struct param_type3 my_params1 = {k1,x1,k2,x2,k3,x3,omega0, par1, par2, par3,karg1,karg2,karg3,karg4};

				gsl_odeiv2_system sys = {funcn1, jac, 18, &my_params1};

			//  this gives the system, the method, the initial interval step, the absolute error and the relative error.
				gsl_odeiv2_driver * d =
				gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
										   1e-3, 1e-3, 0.);

				int status1 = gsl_odeiv2_driver_apply (d, &a, A, G);


				/*Allocation of array values */

			//F1(k;a), G1(k;a)
			F1_nk = G[0];
			G1_nk = G[1];

			// F1(k-p;a), G1(k-p;a)
			F1kmp_nk[i*n2 + j] = G[2];
			G1kmp_nk[i*n2 + j] = G[3];

			//F1(p;a), G1(p;a)
			F1p_nk[i*n2 + j] = G[4];
			G1p_nk[i*n2 + j] = G[5];

			/*2nd order*/

			//F2/G2(p,k-p) (P22)
			F2_nk[i*n2 + j] =  G[6];
			G2_nk[i*n2 + j] =  G[7];

			//F2/G2(-k,k-p)
			F2C_nk[i*n2 + j] =  G[8];
			G2C_nk[i*n2 + j] =  G[9];

      //F2/G2(p,-p)
      F2D_nk[i*n2 + j] =  G[10];
      G2D_nk[i*n2 + j] =  G[11];

			//F2/G2(p,k)
			F2A_nk[i*n2 + j] =  G[12];
			G2A_nk[i*n2 + j] =  G[13];

			//F2/G2(p,-k)
			F2B_nk[i*n2 + j] =  G[14];
			G2B_nk[i*n2 + j] =  G[15];

			/* 3rd order */
			//F3/G3[k,p,-p]
			F3_nk[i*n2 + j]= G[16];
			G3_nk[i*n2 + j]= G[17];

			gsl_odeiv2_driver_free(d);
		}
	}

}



/* Modified friction term equations - see Baldi and Fergus : https://academic.oup.com/mnras/article-lookup/doi/10.1093/mnras/stw2702 */
// variable omega models
// CPL model

static double HAdecpl(double a, double omega0, double w0, double wa){
  double A = -3.*(1.+w0+wa);
  double omegaf = pow(a,A)*exp(3*(-1+a)*wa);
	double omegaL= (1.-omega0)*omegaf;
	return  sqrt(omega0/pow(a,3)+omegaL);
}
// aH dH/da
static double HA1decpl(double a, double omega0, double w0, double wa){
  double A = -3.*(1.+w0+wa);
  double omegaf = pow(a,A)*exp(3*(-1+a)*wa);
	double omegaL= (1.-omega0)*omegaf;
  return -3*omega0/(2.*pow(a,3)) + (A+3*a*wa)*omegaL/2.;
}
//HA2de = -dH/dt/H^2 = -a dH/da / H
static double HA2decpl(double a, double omega0, double w0, double wa){
  return -HA1decpl(a,omega0,w0,wa)/pow2(HAdecpl(a,omega0,w0,wa));
}
//3/(2H^2) * Omega_m
static double HA2de2cpl(double a, double omega0, double w0, double wa){
  return 3.*omega0/(2.*pow2(HAdecpl(a,omega0,w0,wa))*pow(a,3));
}

// fixed omega
// H/H0
static double HAde(double a, double omega0, double p){
  double A = -3.*(1.+p);
  double omegaf = pow(a,A);
	double omegaL= (1.-omega0)*omegaf;
	return  sqrt(omega0/pow(a,3)+omegaL);
}
// aH dH/da / H0^2
static double HA1de(double a, double omega0, double p){
  double A = -3.*(1.+p);
  double omegaf = pow(a,A);
  double omegaL= (1.-omega0)*omegaf;
	return -3*omega0/(2.*pow(a,3)) + A*omegaL/2.;
}
//HA2de = -dH/dt/H^2 = -a dH/da / H
static double HA2de(double a, double omega0, double p){
  return -HA1de(a,omega0,p)/pow2(HAde(a,omega0,p));
}
//3/(2H^2) * Omega_m
static double HA2de2(double a, double omega0, double p){
  return 3.*omega0/(2.*pow2(HAde(a,omega0,p))*pow(a,3));
}



// HYP model
static double omega_hyp(double w0, double wa,double av){
  double whyp = w0+wa/2.*(1+tanh(1./av-3.5));
  return 3.*(1.+whyp)/av;
}
static double omega_var_hyp(double w0, double wa, double a){
  return Integrate(bind(omega_hyp, w0, wa, _1), a, 1, 1e-4);
}
static double HAdehyp(double a, double omega0, double p){
  double omegaf = exp(p);
	double omegaL= (1.-omega0)*omegaf;
	return  sqrt(omega0/pow(a,3)+omegaL);
}
// aH dH/da
static double HA1dehyp(double a, double omega0, double p){
  double omegaf = exp(p);
  double omegaL= (1.-omega0)*omegaf;
	return -3*omega0/(2.*pow(a,3)) + p*a*omegaL/2.;
}
//HA2de = -dH/dt/H^2 = -a dH/da / H
static double HA2dehyp(double a, double omega0, double p){
  return -HA1dehyp(a,omega0,p)/pow2(HAdehyp(a,omega0,p));
}
//3/(2H^2) * Omega_m
static double HA2de2hyp(double a, double omega0, double p){
  return 3.*omega0/(2.*pow2(HAdehyp(a,omega0,p))*pow(a,3));
}


struct param_type4 {
  real kk1;
  real xx1;
	real kk2;
  real xx2;
	real kk3;
  real xx3;
  real arg1;
  real arg2;
  real arg3;
  real arg4;
	real omega00;
	real par1;
	real par2;
	real par3;
  int par4;

};


int funcn2(double a, const double G[], double F[], void *params)
{
	param_type4 p = *(param_type4 *)(params);
	real k1 = p.kk1;
	real x1= p.xx1;
	real k2 = p.kk2;
	real x2= p.xx2;
	real k3 = p.kk3;
	real x3= p.xx3;

  real karg1 = p.arg1;
  real karg2 = p.arg2;
  real karg3 = p.arg3;
  real karg4 = p.arg4;

  real omega0 = p.omega00;
	real p1 = p.par1;
	real p2 = p.par2;
	real p3 = p.par3;
  int p4 = p.par4;


  // Friction term : See Eq. 6-9 of Baldi and Fergus
  // Constant = 3/8pi * 1e28(m^2->bn) / 1.783*10^(27) (kg->Gev/c^2) /G * 1/3.086*10^19 (km/Mpc) * c (dimensionless) = 0.0974655 (so that p1 = h*xi as in Baldi)
  double omegaL,omegaf,A,fric,hade1,hade2;

  switch (p4) {
    case 1:
     A = -3.*(1.+p2);
     omegaf = pow(a,A);
     omegaL= (1.-omega0)*omegaf;
     fric = (1.+p2)*omegaL*p1/HAde(a,omega0,p2)*0.0974655;
     hade1 = HA2de(a,omega0,p2);
     hade2 = HA2de2(a,omega0,p2);
    break;
    case 2:
     A = -3*(1+p2+p3);
     omegaf = pow(a,A)*exp(3.*(-1.+a)*p3);
     omegaL= (1.-omega0)*omegaf;
     fric = (1.+(p2+(1.-a)*p3))*omegaL*p1/HAdecpl(a,omega0,p2,p3)*0.0974655;
     hade1 = HA2decpl(a,omega0,p2,p3);
     hade2 = HA2de2cpl(a,omega0,p2,p3);
    break;
    case 3:
     A = omega_var_hyp(p2,p3,a);
     omegaL= (1.-omega0)*exp(A);
     fric = (1.+(p2+p3/2.*(1.+tanh(1./a-3.5))))*omegaL*p1/HAdehyp(a,omega0,A)*0.0974655;
     hade1 = HA2dehyp(a,omega0,A);
     hade2 = HA2de2hyp(a,omega0,A);
    break;
  }
	/* 1st order */
	//1. F1/G1(k)

	F[0] = -G[1]/a;
	F[1] =1./a*(-(2.+fric-hade1)*G[1]-hade2*G[0]);

	//2. F1/G1(k-p)
	F[2] = -G[3]/a;
	F[3] =1./a*(-(2.+fric-hade1)*G[3]-hade2*G[2]);

	//3. F1/G1(p)
	F[4] = -G[5]/a;
	F[5] =1./a*(-(2.+fric-hade1)*G[5]-hade2*G[4]);

	/* 2nd order */

	//4. F2/G2(p,k-p) (P22)
	F[6] =1./a*(-(alpha(k2,karg4,(k1*x2-k2)/karg4)*G[5]*G[2]+alpha(karg4,k2,(k1*x2-k2)/karg4)*G[3]*G[4])/2.-G[7]);
	F[7] =1./a*(-(2.+fric-hade1)*G[7]-hade2*G[6] - beta1(k2,karg4,(k1*x2-k2)/karg4)*G[5]*G[3]);


	//5. F2/G2(-k,k-p)
	F[8] =1./a*(-(alpha(k1,karg4,(k2*x2-k1)/karg4)*G[1]*G[2]+alpha(karg4,k1,(k2*x2-k1)/karg4)*G[3]*G[0])/2.-G[9]);
	F[9] =1./a*(-(2.+fric-hade1)*G[9]-hade2*G[8] - beta1(k1,karg4,(k2*x2-k1)/karg4)*G[3]*G[1]);


	//6. F2/G2(p,-p)
	F[10] =1./a*(-alpha(k2,k3,x1)*G[5]*G[4]*0.-G[11]) ; // set 1st term to 0 because x1 is not identically -1 (needed to keep evolution at 0 for GR)
	F[11] =1./a*(-(2.+fric-hade1)*G[11]-hade2*G[10] -beta1(k2,k3,x1)*G[5]*G[5]*0.);


	//7. F2/G2(p,k)
	F[12] =1./a*(-(alpha(k2,k1,x2)*G[5]*G[0]+alpha(k1,k2,x2)*G[1]*G[4])/2.-G[13]) ;
	F[13] =1./a*(-(2.+fric-hade1)*G[13]-hade2*G[12] -beta1(k2,k1,x2)*G[5]*G[1]);



	//8. F2/G2(-p,k)=F2/G2(p,-k)
	F[14] =1./a*(-(alpha(k3,k1,x3)*G[5]*G[0]+alpha(k1,k3,x3)*G[1]*G[4])/2.-G[15]) ;
	F[15] =1./a*(-(2.+fric-hade1)*G[15]-hade2*G[14] -beta1(k3,k1,x3)*G[5]*G[1]);


	//9. 3rd order  ~ F3/G3(k,p,-p)
	F[16] = - 1./(3.*a)*(alpha(k1,karg1,(k2*x2+k3*x3)/karg1)*G[10]*G[1]
						+ alpha(k3,karg2,(k1*x3+k2*x1)/karg2)*G[12]*G[5]
						+ alpha(k2,karg3,(k3*x1+k1*x2)/karg3)*G[14]*G[5]

						+ alpha(karg2,k3,(k1*x3+k2*x1)/karg2)*G[13]*G[4]
						+ alpha(karg3,k2,(k3*x1+k1*x2)/karg3)*G[15]*G[4]
						+ alpha(karg1,k1,(k2*x2+k3*x3)/karg1)*G[11]*G[0]

						+3.*G[17]) ;


	F[17] =1./(3.*a)*(-3.*(2.+fric-hade1)*G[17]-3.*hade2*G[16]

					 -2.*beta1(k1,karg1,(k2*x2+k3*x3)/karg1)*G[1]*G[11]

					 -2.*beta1(k3,karg2,(k1*x3+k2*x1)/karg2)*G[5]*G[13]

					 -2.*beta1(k2,karg3,(k3*x1+k1*x2)/karg3)*G[5]*G[15]);


	return GSL_SUCCESS;
}


// For wCDM normalization
int funcn3(double a, const double G[], double F[], void *params)
{
	param_type2 p = *(param_type2 *)(params);
	real omega0 = p.omega00;
	real p1 = p.par1;
  real p2;
	int p3 = p.par4;
  if (p3>1.1) {
    p2 = -0.8;
  }
  else{
    p2 = p.par2;
  }
	/* 1st order */
	//1. F1/G1(k)
	F[0] = -G[1]/a;
	F[1] =1./a*(-(2.-HA2de(a,omega0,p2))*G[1]-HA2de2(a,omega0,p2)*G[0]);
	return GSL_SUCCESS;
}



void IOW::initn_de(double A, double YMIN, double YMAX, double k, double omega0, double par1, double par2, double par3, int par4)
{
//#pragma omp parallel for schedule(dynamic)
    real x1 = XMIN;
    real k1 = k;
		for (int j = 0; j< n2; ++j)
		{
    real k2,k3,karg1;
  /* k - magnitude sampling */
    double R = (YMIN) * exp(j*log(YMAX/YMIN)/(n2*1.-1.)); //exponential sampling
    k2 = k*R;
    k3 = k2;
    karg1 = sqrt(k2*k2+2.*k2*k3*x1+k3*k3);

      for ( int i = 0; i < n1; ++i)
      {

        real x2,x3,karg2,karg3,karg4;
				/*Parameter selection for kernels */

        x2 = x128[i]; // x128 goes from - 0.99995 to 0.99995
				x3 = -x2;

        karg2 = sqrt(k1*k1+2.*k1*k2*x2+k2*k2);
        karg3 = sqrt(k3*k3+2.*k3*k1*x3+k1*k1);
        karg4 = sqrt(k2*k2+k1*k1-2.*k2*k1*x2);


			// Initial scale factor for solving system of equations
			//Sets 1-3
				double a = 0.00001;

// Non-Eds ICs
			  double G[18] = { a,-a,a,-a,a,-a, 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

			/*Parameters passed to system of equations */
// EDIT : If more than one gravity parameter is used, add them after p1
				struct param_type4 my_params1 = {k1,x1,k2,x2,k3,x3,karg1,karg2,karg3,karg4, omega0, par1, par2, par3,par4};

				gsl_odeiv2_system sys = {funcn2, jac, 18, &my_params1};

			//  this gives the system, the method, the initial interval step, the absolute error and the relative error.
				gsl_odeiv2_driver * d =
				gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
										   1e-3, 1e-3, 0.0);

				int status1 = gsl_odeiv2_driver_apply (d, &a, A, G);


				/*Allocation of array values */

			//F1(k;a), G1(k;a)
			F1_nk = G[0];
			G1_nk = G[1];

			// F1(k-p;a), G1(k-p;a)
			F1kmp_nk[i*n2 + j] = G[2];
			G1kmp_nk[i*n2 + j] = G[3];

			//F1(p;a), G1(p;a)
			F1p_nk[i*n2 + j] = G[4];
			G1p_nk[i*n2 + j] = G[5];

			/*2nd order*/

			//F2/G2(p,k-p) (P22)
			F2_nk[i*n2 + j] =  G[6];
			G2_nk[i*n2 + j] =  G[7];

			//F2/G2(-k,k-p)
			F2C_nk[i*n2 + j] =  G[8];
			G2C_nk[i*n2 + j] =  G[9];

      //F2/G2(p,-p)
      F2D_nk[i*n2 + j] =  G[10];
      G2D_nk[i*n2 + j] =  G[11];

			//F2/G2(p,k)
			F2A_nk[i*n2 + j] =  G[12];
			G2A_nk[i*n2 + j] =  G[13];

			//F2/G2(p,-k)
			F2B_nk[i*n2 + j] =  G[14];
			G2B_nk[i*n2 + j] =  G[15];

			/* 3rd order */
			//F3/G3[k,p,-p]
			F3_nk[i*n2 + j]= G[16];
			G3_nk[i*n2 + j]= G[17];

			gsl_odeiv2_driver_free(d);
		}
	}
//  Solution for wCDM linear growth @ a=1 for our normalisation of the linear power spectrum
  	double a = 0.0001;
  	double G[2] = {a,-a}; // initial conditions
  	struct param_type2 params_my = {omega0,par1,par2,1.,par4};

  // Solutions of evolution factors @ a=A
  	gsl_odeiv2_system sys3 = {funcn3, jac, 2, &params_my};
  	gsl_odeiv2_driver * d1 =
  	gsl_odeiv2_driver_alloc_y_new (&sys3, gsl_odeiv2_step_rk8pd,
  								   1e-3, 1e-3, 0.0);

  	int status2 = gsl_odeiv2_driver_apply (d1, &a, 1. , G);

  	dnorm_spt = G[0];

  	gsl_odeiv2_driver_free(d1);

}


/// the new solver for numerical equations - optimised and solved in SPT.cpp's ploopn2 functions

void IOW::initn2(double A, double k[], double x[], double kargs[], double omega0, double par1, double par2, double par3)
{

			// Initial scale factor for solving system of equations
				double a = 0.00001;

        // Non-Eds ICs
			  double G[18] = { a,-a,a,-a,a,-a, 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

			/*Parameters passed to system of equations */
      // EDIT : If more than one gravity parameter is used, add them after p1
				struct param_type3 my_params1 = {k[0],x[0],k[1],x[1],k[2],x[2],omega0, par1, par2, par3,kargs[0],kargs[1],kargs[2],kargs[3]};

				gsl_odeiv2_system sys = {funcn1, jac, 18, &my_params1};

			//  this gives the system, the method, the initial interval step, the absolute error and the relative error.
				gsl_odeiv2_driver * d =
				gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
										   1e-5, 1e-4, 1e-2);

				int status1 = gsl_odeiv2_driver_apply (d, &a, A, G);


				/*Allocation of array values */

			//F1(k;a), G1(k;a)
			F1_nk = G[0];
			G1_nk = G[1];

			// F1(k-p;a), G1(k-p;a)
			F1kmp_nk[0] = G[2];
			G1kmp_nk[0] = G[3];

			//F1(p;a), G1(p;a)
			F1p_nk[0] = G[4];
			G1p_nk[0] = G[5];

			/*2nd order*/

			//F2/G2(p,k-p) (P22)
			F2_nk[0] =  G[6];
			G2_nk[0] =  G[7];

			//F2/G2(-k,k-p)
			F2C_nk[0] =  G[8];
			G2C_nk[0] =  G[9];

      //F2/G2(p,-p)
      F2D_nk[0] =  G[10];
      G2D_nk[0] =  G[11];

			//F2/G2(p,k)
			F2A_nk[0] =  G[12];
			G2A_nk[0] =  G[13];

			//F2/G2(p,-k)
			F2B_nk[0] =  G[14];
			G2B_nk[0] =  G[15];

			/* 3rd order */
			//F3/G3[k,p,-p]
			F3_nk[0]= G[16];
			G3_nk[0]= G[17];

			gsl_odeiv2_driver_free(d);
}


void IOW::initn2_ide(double A, double k[], double x[], double kargs[], double omega0, double par1, double par2, double par3, int par4)
{
			// Initial scale factor for solving system of equations
				double a = 0.00001;

        // Non-Eds ICs
			  double G[18] = { a,-a,a,-a,a,-a, 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

			/*Parameters passed to system of equations */
      // EDIT : If more than one gravity parameter is used, add them after p1

				struct param_type4 my_params1 = {k[0],x[0],k[1],x[1],k[2],x[2], kargs[2],kargs[3],kargs[1],kargs[0], omega0, par1, par2, par3,par4};

				gsl_odeiv2_system sys = {funcn2, jac, 18, &my_params1};


			//  this gives the system, the method, the initial interval step, the absolute error and the relative error.
				gsl_odeiv2_driver * d =
				gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
										   1e-3, 1e-3, 0.0);

				int status1 = gsl_odeiv2_driver_apply (d, &a, A, G);


				/*Allocation of array values */

			//F1(k;a), G1(k;a)
			F1_nk = G[0];
			G1_nk = G[1];

			// F1(k-p;a), G1(k-p;a)
			F1kmp_nk[0] = G[2];
			G1kmp_nk[0] = G[3];

			//F1(p;a), G1(p;a)
			F1p_nk[0] = G[4];
			G1p_nk[0] = G[5];

			/*2nd order*/

			//F2/G2(p,k-p) (P22)
			F2_nk[0] =  G[6];
			G2_nk[0] =  G[7];

			//F2/G2(-k,k-p)
			F2C_nk[0] =  G[8];
			G2C_nk[0] =  G[9];

      //F2/G2(p,-p)
      F2D_nk[0] =  G[10];
      G2D_nk[0] =  G[11];

			//F2/G2(p,k)
			F2A_nk[0] =  G[12];
			G2A_nk[0] =  G[13];

			//F2/G2(p,-k)
			F2B_nk[0] =  G[14];
			G2B_nk[0] =  G[15];

			/* 3rd order */
			//F3/G3[k,p,-p]
			F3_nk[0]= G[16];
			G3_nk[0]= G[17];

			gsl_odeiv2_driver_free(d);

}


// linear growth factor for w=fixed with interaction
int funcn4_ide(double a, const double G[], double F[], void *params)
{
	param_type2 p = *(param_type2 *)(params);
	real omega0 = p.omega00;
	real p1 = p.par1;
	real p2 = p.par2;
  real p3 = p.par3;
  int p4 = p.par4;

		double A, omegaf, omegaL, fric,hade1,hade2;
    switch (p4) {
      case 1:
       A = -3.*(1.+p2);
       omegaf = pow(a,A);
       omegaL= (1.-omega0)*omegaf;
       fric = (1.+p2)*omegaL*p1/HAde(a,omega0,p2)*0.0974655;
       hade1 = HA2de(a,omega0,p2);
       hade2 = HA2de2(a,omega0,p2);
      break;
      case 2:
       A = -3*(1+p2+p3);
       omegaf = pow(a,A)*exp(3.*(-1.+a)*p3);
       omegaL= (1.-omega0)*omegaf;
       fric = (1.+(p2+(1.-a)*p3))*omegaL*p1/HAdecpl(a,omega0,p2,p3)*0.0974655;
       hade1 = HA2decpl(a,omega0,p2,p3);
       hade2 = HA2de2cpl(a,omega0,p2,p3);
      break;
      case 3:
       A = omega_var_hyp(p2,p3,a);
       omegaL= (1.-omega0)*exp(A);
       fric = (1.+(p2+p3/2.*(1.+tanh(1./a-3.5))))*omegaL*p1/HAdehyp(a,omega0,A)*0.0974655;
       hade1 = HA2dehyp(a,omega0,A);
       hade2 = HA2de2hyp(a,omega0,A);
      break;
    }

  	F[0] = -G[1]/a;
	  F[1] =1./a*(-(2.+fric-hade1)*G[1]-hade2*G[0]);

	return GSL_SUCCESS;
}


void IOW::initnorm_ide(double A, double omega0, double par1, double par2, double par3, int par4)
{
			//  Solution for wCDM linear growth @ a=1 for our normalisation of the linear power spectrum

					struct param_type2 params_my = {omega0,par1,par2,par3,par4};

					double a1 = 0.0001;
			  	double G1[2] = {a1,-a1}; // initial conditions

			  // Solutions of evolution factors @ a=A
			  	gsl_odeiv2_system sys3 = {funcn3, jac, 2, &params_my};
			  	gsl_odeiv2_driver * d1 =
			  	gsl_odeiv2_driver_alloc_y_new (&sys3, gsl_odeiv2_step_rk8pd,
			  								   1e-3, 1e-3, 0.);

			  	int status2 = gsl_odeiv2_driver_apply (d1, &a1, 1. , G1);

			  	dnorm_spt = G1[0];

			  	gsl_odeiv2_driver_free(d1);


					double a = 0.0001;

					// Non-Eds ICs
					double G[2] = { a,-a};

					gsl_odeiv2_system sys = {funcn4_ide, jac, 2, &params_my};

				//  this gives the system, the method, the initial interval step, the absolute error and the relative error.
					gsl_odeiv2_driver * d =
					gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
												 1e-3, 1e-3, 0.);

					int status1 = gsl_odeiv2_driver_apply (d, &a, A, G);

					D_spt = G[0];
          Dl_spt = G[0];

					gsl_odeiv2_driver_free(d);

}




/// GAUGE INVARIANT EXPRESSIONS OF Anselmi +  PIETRONI //////

// set of equations for CDE bispectrum using Pietroni Gauge-inv expressions (3.5 and 3.6 of 1106.0834v3)

int funcn4_ginvbs(double a, const double G[], double F[], void *params)
{
  param_type3 p = *(param_type3 *)(params);
	real k1 = p.kk1;
	real x1= p.xx1;
	real k2 = p.kk2;
	real x2= p.xx2;
	real k3 = p.kk3;
	real x3= p.xx3;
	real omega0 = p.omega00;
  // p1 = cs2, p2 = w
	real p1 = p.par1;
	real p2 = p.par2;
	real p3 = p.par3;
  real karg1 = p.arg1;


// hubble functions and omega_m
  double hade1,hade2,hade3,hade4,A,omegaf,omegaL,h0;
  h0 = 1./2997.9;
  hade1 = HA2de(a,omega0,p2);
  hade2 = HA2de2(a,omega0,p2);
  A = -3.*(1.+p2);
  omegaf = pow(a,A);
  omegaL= (1.-omega0)*omegaf;
  hade3 = p3*pow3(a)*omegaL*hade2/omega0;
  hade4 = pow2(HAde(a,omega0,p2)*h0*a);

  double a1,a2,a3,a4,a5,a6;
  double b1,b2,b3;
  double c1;
  c1 = p1/(hade4*(1.+1e-8+p2));

// different k arguments

  a1 = alpha(k2,karg1,(-k1*x2-k2)/karg1);
  a2 = alpha(karg1,k2,(-k1*x2-k2)/karg1);
  b1 = beta1(k2,karg1,(-k1*x2-k2)/karg1);

  a3 = alpha(k1,karg1,(-k2*x2-k1)/karg1);
  a4 = alpha(karg1,k1,(-k2*x2-k1)/karg1);
  b2 = beta1(k1,karg1,(-k2*x2-k1)/karg1);

  a5 = alpha(k2,k1,x2);
  a6 = alpha(k1,k2,x2);
  b3 = beta1(k1,k2,x2);

  //  a3 = alpha(k2,k2,x4);
  //  b4 = beta1(k2,k2,x4);

	/* 1st order */
	//1. F1/G1(k)
	F[0] = -G[1]/a;
	F[1] =1./a*(-(2.-hade1)*G[1]-hade2*G[0]-hade3*G[2]);
  F[2] =1./a*(-(1+p2)*G[3]+3*(p2-p1)*G[2] + 9.*(1+p2)*(p2-p1)*hade4/pow2(k1)*G[3]);
  F[3] =1./a*(-(2.-3*p1-hade1)*G[3]-hade2*G[0]-(hade3 - c1*pow2(k1))*G[2]);

  //2. F1/G1(-k-p)
  F[4] = -G[5]/a;
	F[5] =1./a*(-(2.-hade1)*G[5]-hade2*G[4]-hade3*G[6]);
  F[6] =1./a*(-(1+p2)*G[7]+3*(p2-p1)*G[6] + 9.*(1+p2)*(p2-p1)*hade4/pow2(karg1)*G[7]  );
  F[7] =1./a*(-(2.-3*p1-hade1)*G[7]-hade2*G[4]-(hade3 - c1*pow2(karg1))*G[6]);

	//3. F1/G1(p)
  F[8] = -G[9]/a;
	F[9] =1./a*(-(2.-hade1)*G[9]-hade2*G[8]-hade3*G[10]);
  F[10] =1./a*(-(1+p2)*G[11]+3*(p2-p1)*G[10] + 9.*(1+p2)*(p2-p1)*hade4/pow2(k2)*G[11] );
  F[11] =1./a*(-(2.-3*p1-hade1)*G[11]-hade2*G[8]-(hade3 - c1*pow2(k2))*G[10]);

/* 2nd order */
  //4. F2/G2(p,-k-p) (P22)
  F[12] =1./a*(-(a1*G[9]*G[4]+a2*G[5]*G[8])/2.-G[13]);
	F[13] =1./a*(-(2.-hade1)*G[13]-hade2*G[12]-hade3*G[14]-b1*G[9]*G[5]);
  F[14] =1./a*(-(1+p1)*(a1*G[11]*G[6]+a2*G[7]*G[10])/2.-(1+p2)*G[15]+3*(p2-p1)*G[14] + 9.*(1+p2)*(p2-p1)*hade4/pow2(k1)*G[15] );
	F[15] =1./a*(-(2.-3*p1-hade1)*G[15]-hade2*G[12]-(hade3 - c1*pow2(k1))*G[14]- (1.-p1)*b1*G[11]*G[7]);
            //      +3.*(p1-p2)/(1.+1e-8+p2)*(a1*G[11]*G[6]+a2*G[7]*G[10])/2. - p1*a/(1.+1e-8+p2)*(a1*G[11]*F[6]+a2*G[7]*F[10])/2.
          //        -4.*(p1+1)*c1/(1.+1e-8+p2)*(-4.*pow2(k2)-4.*k2*k1*x2-pow2(k1))*G[6]*G[10]);

  //5. F2/G2[k,-k-p]

  F[16] =1./a*(-(a3*G[1]*G[4]+a4*G[5]*G[0])/2.-G[17]);
  F[17] =1./a*(-(2.-hade1)*G[17]-hade2*G[16]-hade3*G[18]-b2*G[1]*G[5]);

  F[18] =1./a*(-(1+p1)*(a3*G[3]*G[6]+a4*G[7]*G[2])/2.-(1+p2)*G[19]+3*(p2-p1)*G[18] + 9.*(1+p2)*(p2-p1)*hade4/pow2(k2)*G[19] );
  F[19] =1./a*(-(2.-3*p1-hade1)*G[19]-hade2*G[16]-(hade3 - c1*pow2(k2))*G[18]- (1.-p1)*b2*G[3]*G[7]);
      //          +3.*(p1-p2)/(1.+1e-8+p2)*(a3*G[3]*G[6]+a4*G[7]*G[2])/2. - p1*a/(1.+1e-8+p2)*(a3*G[3]*F[6]+a4*G[7]*F[2])/2.
      //          -4.*(p1+1)*c1/(1.+1e-8+p2)*(-4.*pow2(k1)-4.*k2*k1*x2-pow2(k2))*G[6]*G[2]);


  //6. F2/G2(p,k)
	F[20] =1./a*(-(a5*G[9]*G[0]+a6*G[1]*G[8])/2.-G[21]) ;
	F[21] =1./a*(-(2.-hade1)*G[21]-hade2*G[20]-hade3*G[22]-b3*G[9]*G[1]);
  F[22] =1./a*(-(1+p1)*(a5*G[11]*G[2]+a6*G[3]*G[10])/2.-(1+p2)*G[23]+3*(p2-p1)*G[22] + 9.*(1+p2)*(p2-p1)*hade4/pow2(karg1)*G[23]) ;
	F[23] =1./a*(-(2.-3.*p1-hade1)*G[23]-hade2*G[20]-(hade3 - c1*pow2(karg1))*G[22]-(1.-p1)*b3*G[11]*G[3]);
        //            +3.*(p1-p2)/(1.+1e-8+p2)*(a5*G[11]*G[2]+a6*G[3]*G[10])/2. - p1*a/(1.+1e-8+p2)*(a5*G[11]*F[2]+a6*G[3]*F[10])/2.
        //            -4.*(p1+1)*c1/(1.+1e-8+p2)*(k2*k2+k1*k1-2*k2*k1*x2)/2.*G[2]*G[10]);


                    //2. F1/G1(k-p)
                    // F[4] = -G[5]/a;
                  	// F[5] =1./a*(-(2.-hade1)*G[5]-hade2*G[4]-hade3*G[6]);
                    // F[6] =1./a*(-(1+p2)*G[7]+3*(p2-p1)*G[6]);
                    // F[7] =1./a*(-(2.-3*p2-hade1)*G[7]-hade2*G[4]-(hade3*(1.+3.*p1) - c1*(k2*k2+k1*k1-2*k2*k1*x2))*G[6]);

                      // //5. F2/G2(p,p)
                    	// F[16] =1./a*(-a3*G[8]*G[9]-G[17]) ;
                    	// F[17] =1./a*(-(2.-hade1)*G[17]-hade2*G[16]-hade3*G[18]-b4*G[9]*G[9]);
                      // F[18] =1./a*(-(1+p1)*a3*G[10]*G[11]-(1+p2)*G[19]+3.*(p2-p1)*G[18]);
                    	// F[19] =1./a*(-(2.-3.*p2-hade1)*G[19]-hade2*G[16]-(hade3*(1+3.*p1) - c1*(2.*pow2(k2)*(1.+x4)))*G[18] - b4*G[11]*G[11]
                      //               +3.*(p1-p2)/(1.+1e-8+p2)*a3*G[11]*G[10] - p1*a/(1.+1e-8+p2)*a3*G[11]*F[10]
                      //               -4.*(p1+1)*c1/(1.+1e-8+p2)*(pow2(k2)*(1.-x4))*G[6]*G[10]);

	return GSL_SUCCESS;
}


//Bispectrum kernel initialisation using Pietroni expressions (3.5 and 3.6 of 1106.0834v3)
// kernel initialization for tree level bispectrum (isosceles and equilateral only)
// angle is angle between k and p vectors (this determines angle between p vectors for isosceles config)

void IOW::initn_cde_btree(double A, double kmag, double pmag, double x, double omega0, double par1, double par2, double par3)
{

  for(int i = 0; i<1; i++){
				real x1,x2,x3,x4,k1,k2,k3,karg1,karg2,karg3,karg4;

				/* k - magnitude sampling */
			//	double R = (YMIN) * exp(j*log(YMAX/YMIN)/(n2*1.-1.)); //exponential sampling

				/*Parameter selection for kernels */

				k1 = kmag;
		    k2 = pmag;//k*R;
				k3 = k2;

				x1 = XMIN;
        x2 = x; //-kmag/(2.*pmag);//x128[i]; // x128 goes from - 0.99995 to 0.99995
				x3 = -x2;

        karg1 = sqrt(k2*k2+k1*k1+2*k2*k1*x2);

				double a = 0.0001;

// Non-Eds ICs
			  double G[24] = { a,-a,a,-a,a,-a,a,-a,a,-a,a,-a, 0.,0.,0.,0., 0.,0.,0.,0., 0.,0.,0.,0.};

			/*Parameters passed to system of equations */
// EDIT : If more than one gravity parameter is used, add them after p1
				struct param_type3 my_params1 = {k1,x1,k2,x2,k3,x3,omega0, par1, par2, par3,karg1,1.,1.,1.};
				gsl_odeiv2_system sys = {funcn4_ginvbs, jac, 24, &my_params1};

			//  this gives the system, the method, the initial interval step, the absolute error and the relative error.
				gsl_odeiv2_driver * d =
				gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
										   1e-3, 1e-3, 1e-3);

				int status1 = gsl_odeiv2_driver_apply (d, &a, A, G);


				/*Allocation of array values */

			//F1(k;a), G1(k;a)
			F1_nk = G[0] ;
			G1_nk = G[1] ;

      //F1q(k;a), G1q(k;a)
      F1q_nk = G[2];
      G1q_nk = G[3];

			// F1(-k-p;a), G1(-k-p;a)
			F1kmp_nk[0] = G[4];
			G1kmp_nk[0] = G[5];

      F1kmpq_nk[0] = G[6];
      G1kmpq_nk[0] = G[7];

			//F1(p;a), G1(p;a)
			F1p_nk[0] = G[8];
			G1p_nk[0] = G[9];

      //F1q(p;a), G1q(p;a)
      F1pq_nk[0] = G[10];
      G1pq_nk[0] = G[11];

			/*2nd order*/

			//F2/G2(p,-k-p) (P22)
			F2_nk[0] =  G[12];
			G2_nk[0] =  G[13];

      F2q_nk[0] =  G[14];
			G2q_nk[0] =  G[15];

      //F2/G2(k,-k-p)
      F2D_nk[0] =  G[16];
      G2D_nk[0] =  G[17];

      F2Dq_nk[0] =  G[18];
      G2Dq_nk[0] =  G[19];

      //F2/G2(p,k)
      F2A_nk[0] =  G[20];
      G2A_nk[0] =  G[21];

      F2Aq_nk[0]= G[22];
      G2Aq_nk[0]= G[23];


			gsl_odeiv2_driver_free(d);
}

double a = 0.0001;
  //  Solution for wCDM linear growth @ a=1 for our normalisation of the linear power spectrum
    	double G[2] = {a,-a}; // initial conditions
    	struct param_type2 params_my = {omega0,par1,par2,0.,1};

    // Solutions of evolution factors @ a=A
    	gsl_odeiv2_system sys3 = {funcn3, jac, 2, &params_my};
    	gsl_odeiv2_driver * d1 =
    	gsl_odeiv2_driver_alloc_y_new (&sys3, gsl_odeiv2_step_rk8pd,
    								   1e-3, 1e-3, 0.0);

    	int status2 = gsl_odeiv2_driver_apply (d1, &a, 1. , G);

    	dnorm_spt = G[0];
      gnorm_spt = G[1];


    	gsl_odeiv2_driver_free(d1);

}

// set of equations for CDE 1-LOOP PS using Pietroni Gauge-inv expressions (3.5 and 3.6 of 1106.0834v3)

int funcn4_ginvps(double a, const double G[], double F[], void *params)
{
  param_type3 p = *(param_type3 *)(params);
	real k1 = p.kk1;
	real x1= p.xx1;
	real k2 = p.kk2;
	real x2= p.xx2;
	real k3 = p.kk3;
	real x3= p.xx3;
	real omega0 = p.omega00;
  // p1 = cs2, p2 = w , p3 was being used to switch off problematic terms so set to 1. normally
	real p1 = p.par1;
	real p2 = p.par2;
	real p3 = p.par3;
  real karg1 = p.arg1;
  real karg2 = p.arg2;

// hubble functions and omega_m
  double hade1,hade2,hade3,hade4,A,omegaf,omegaL,h0;
  h0 = 1./2997.9;
  hade1 = HA2de(a,omega0,p2);
  hade2 = HA2de2(a,omega0,p2);
  A = -3.*(1.+p2);
  omegaf = pow(a,A);
  omegaL= (1.-omega0)*omegaf;
  hade3 = p3*pow3(a)*omegaL*hade2/omega0;
  hade4 = pow2(HAde(a,omega0,p2)*h0*a);

  double a1,a2,a3,a4,a5,a6,a7,a8,a9,a10;
  double b1,b2,b3,b4,b5;
  double c1;
  c1 = p1/(hade4*(1.+1e-8+p2));

// different k arguments

  a1 = alpha(k2,karg1,(k1*x2-k2)/karg1);
  a2 = alpha(karg1,k2,(k1*x2-k2)/karg1);
  b1 = beta1(k2,karg1,(k1*x2-k2)/karg1);

  a3 = alpha(k1,k2,x2);
  a4 = alpha(k2,k1,x2);
  b2 = beta1(k1,k2,x2);

  a5 = alpha(k2,k1,-x2);
  a6 = alpha(k1,k2,-x2);
  b3 = beta1(k1,k2,-x2);

  a7 = alpha(k2,karg2,(k1*x2+k2)/karg2);
  a8 = alpha(k2,karg2,-(k1*x2+k2)/karg2);
  a9 = alpha(karg2,k2,(k1*x2+k2)/karg2);
  a10 = alpha(karg2,k2,-(k1*x2+k2)/karg2);

  b4 = beta1(k2,karg2,(k1*x2+k2)/karg2);
  b5 = beta1(karg2,k2,-(k1*x2+k2)/karg2);


	/* 1st order */
	//1. F1/G1(k)
	F[0] = -G[1]/a;
	F[1] =1./a*(-(2.-hade1)*G[1]-hade2*G[0]-hade3*G[2]);
  F[2] =1./a*(-(1+p2)*G[3]+3*(p2-p1)*G[2] + 9.*(1+p2)*(p2-p1)*hade4/pow2(k1)*G[3]);
  F[3] =1./a*(-(2.-3*p1-hade1)*G[3]-hade2*G[0]-(hade3 - c1*pow2(k1))*G[2]);

  //2. F1/G1(k-p)
  F[4] = -G[5]/a;
	F[5] =1./a*(-(2.-hade1)*G[5]-hade2*G[4]-hade3*G[6]);
  F[6] =1./a*(-(1+p2)*G[7]+3*(p2-p1)*G[6] + 9.*(1+p2)*(p2-p1)*hade4/pow2(karg1)*G[7]  );
  F[7] =1./a*(-(2.-3*p1-hade1)*G[7]-hade2*G[4]-(hade3 - c1*pow2(karg1))*G[6]);

	//3. F1/G1(p)
  F[8] = -G[9]/a;
	F[9] =1./a*(-(2.-hade1)*G[9]-hade2*G[8]-hade3*G[10]);
  F[10] =1./a*(-(1+p2)*G[11]+3*(p2-p1)*G[10] + 9.*(1+p2)*(p2-p1)*hade4/pow2(k2)*G[11] );
  F[11] =1./a*(-(2.-3*p1-hade1)*G[11]-hade2*G[8]-(hade3 - c1*pow2(k2))*G[10]);

/* 2nd order */
  //4. F2/G2(p,k-p) (P22)
  F[12] =1./a*(-(a1*G[9]*G[4]+a2*G[5]*G[8])/2.-G[13]);
	F[13] =1./a*(-(2.-hade1)*G[13]-hade2*G[12]-hade3*G[14]-b1*G[9]*G[5]);
  F[14] =1./a*(-(1+p1)*(a1*G[11]*G[6]+a2*G[7]*G[10])/2.-(1+p2)*G[15]+3*(p2-p1)*G[14] + 9.*(1+p2)*(p2-p1)*hade4/pow2(k1)*G[15] );
	F[15] =1./a*(-(2.-3*p1-hade1)*G[15]-hade2*G[12]-(hade3 - c1*pow2(k1))*G[14]- (1.-p1)*b1*G[11]*G[7]);
            //      +3.*(p1-p2)/(1.+1e-8+p2)*(a1*G[11]*G[6]+a2*G[7]*G[10])/2. - p1*a/(1.+1e-8+p2)*(a1*G[11]*F[6]+a2*G[7]*F[10])/2.
          //        -4.*(p1+1)*c1/(1.+1e-8+p2)*(-4.*pow2(k2)-4.*k2*k1*x2-pow2(k1))*G[6]*G[10]);


// 2nd order kernels used in 3rd order equations :

  //5. F2/G2[p,k]

  F[16] =1./a*(-(a3*G[1]*G[4]+a4*G[5]*G[0])/2.-G[17]);
  F[17] =1./a*(-(2.-hade1)*G[17]-hade2*G[16]-hade3*G[18]-b2*G[1]*G[5]);

  F[18] =1./a*(-(1+p1)*(a3*G[3]*G[6]+a4*G[7]*G[2])/2.-(1+p2)*G[19]+3*(p2-p1)*G[18] + 9.*(1+p2)*(p2-p1)*hade4/pow2(karg2)*G[19] );
  F[19] =1./a*(-(2.-3*p1-hade1)*G[19]-hade2*G[16]-(hade3 - c1*pow2(karg2))*G[18]- (1.-p1)*b2*G[3]*G[7]);
      //          +3.*(p1-p2)/(1.+1e-8+p2)*(a3*G[3]*G[6]+a4*G[7]*G[2])/2. - p1*a/(1.+1e-8+p2)*(a3*G[3]*F[6]+a4*G[7]*F[2])/2.
      //          -4.*(p1+1)*c1/(1.+1e-8+p2)*(-4.*pow2(k1)-4.*k2*k1*x2-pow2(k2))*G[6]*G[2]);


  //6. F2/G2(p,-k)
  F[20] =1./a*(-(a5*G[9]*G[0]+a6*G[1]*G[8])/2.-G[21]) ;
  F[21] =1./a*(-(2.-hade1)*G[21]-hade2*G[20]-hade3*G[22]-b3*G[9]*G[1]);
  F[22] =1./a*(-(1+p1)*(a5*G[11]*G[2]+a6*G[3]*G[10])/2.-(1+p2)*G[23]+3*(p2-p1)*G[22] + 9.*(1+p2)*(p2-p1)*hade4/pow2(karg1)*G[23]) ;
  F[23] =1./a*(-(2.-3.*p1-hade1)*G[23]-hade2*G[20]-(hade3 - c1*pow2(karg1))*G[22]-(1.-p1)*b3*G[11]*G[3]);
        //            +3.*(p1-p2)/(1.+1e-8+p2)*(a5*G[11]*G[2]+a6*G[3]*G[10])/2. - p1*a/(1.+1e-8+p2)*(a5*G[11]*F[2]+a6*G[3]*F[10])/2.
        //            -4.*(p1+1)*c1/(1.+1e-8+p2)*(k2*k2+k1*k1-2*k2*k1*x2)/2.*G[2]*G[10]);

// HAVE OMITTED F2(P,-P) because of singularity in 1/k^2 term - should not affect final result

//9. 3rd order  ~ F3/G3(k,p,-p)
F[24] = - 1./(3.*a)*(a7*G[20]*G[9] // a(p,k-p) F2(k,-p) G1(p)
					        + a8*G[16]*G[9] // a(-p,k+p) F2(k,p) G1(p)

					        + a9*G[21]*G[8]  // a(k-p,p) G2(k,-p) F1(p)
					        + a10*G[17]*G[8]  // a(k+p,-p) G2(k,p) F1(p)

					+3.*G[25]) ;


F[25] =1./(3.*a)*(-3.*(2.-hade1)*G[25]-3.*hade2*G[24] - 3*hade3*G[26]
				 -2*b4*G[9]*G[21] // b(k-p,p) G2(k,-p) G1(p)
				 -2*b5*G[9]*G[17] // b(k+p,-p) G2(k,p) G1(p)
            );

F[26] =  1./(3.*a)*(-(1+p1)*(a7*G[22]*G[11]
                        +   a8*G[18]*G[11]
                        +   a9*G[23]*G[10]
                        +   a10*G[19]*G[10])

            -3.*(1+p2)*G[27]+3.*(p2-p1)*G[26] + 3.*9.*(1+p2)*(p2-p1)*hade4/pow2(k1)*G[27]) ;


F[27] =1./(3.*a)*(-3.*(2.-3*p1-hade1)*G[27]-3.*hade2*G[24]-(hade3 - c1*pow2(k1))*G[26]
            +(1.-p1)*(-2*b4*G[11]*G[23]
                      -2*b5*G[11]*G[19])
          );


	return GSL_SUCCESS;
}


void IOW::initn_cde_ps(double A, double YMIN, double YMAX, double k, double omega0, double par1, double par2, double par3)
{
//#pragma omp parallel for schedule(dynamic)
for (int j = 0; j< n2; ++j)
{
  /* k - magnitude sampling */
  double R = (YMIN) * exp(j*log(YMAX/YMIN)/(n2*1.-1.)); //exponential sampling

	for ( int i = 0; i < n1; ++i)
	{
				real x1,x2,x3,k1,k2,k3,karg1,karg2;

				/*Parameter selection for kernels */

				k1 = k;
		    k2 = k*R;
				k3 = k2;

				x1 = XMIN;
        x2 = x128[i]; // x128 goes from - 0.99995 to 0.99995
				x3 = -x2;

        karg1 = sqrt(k2*k2+k1*k1-2.*k2*k1*x2);
        karg2 = sqrt(k2*k2+k1*k1+2.*k2*k1*x2);

				double a = 0.0001;


// Non-Eds ICs
			  double G[28] = { a,-a,a,-a,a,-a,a,-a,a,-a,a,-a, 0.,0.,0.,0., 0.,0.,0.,0., 0.,0.,0.,0.,0.,0.,0.,0.};

			/*Parameters passed to system of equations */
// EDIT : If more than one gravity parameter is used, add them after p1
				struct param_type3 my_params1 = {k1,x1,k2,x2,k3,x3,omega0, par1, par2, par3,karg1,karg2,1.,1.};

				gsl_odeiv2_system sys = {funcn4_ginvps, jac, 28, &my_params1};

			//  this gives the system, the method, the initial interval step, the absolute error and the relative error.
				gsl_odeiv2_driver * d =
				gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
										   1e-3, 1e-3, 1e-2);

				int status1 = gsl_odeiv2_driver_apply (d, &a, A, G);


				/*Allocation of array values */

			//F1(k;a), G1(k;a)
			F1_nk = G[0] ;
			G1_nk = G[1] ;

      //F1q(k;a), G1q(k;a)
      F1q_nk = G[2];
      G1q_nk = G[3];

			// F1(k-p;a), G1(k-p;a)
			F1kmp_nk[i*n2 + j] = G[4];
			G1kmp_nk[i*n2 + j] = G[5];

			//F1(p;a), G1(p;a)
			F1p_nk[i*n2 + j] = G[8];
			G1p_nk[i*n2 + j] = G[9];

      //F1q(p;a), G1q(p;a)
      F1pq_nk[i*n2 + j] = G[10];
      G1pq_nk[i*n2 + j] = G[11];

			/*2nd order*/

			//F2/G2(p,k-p) (P22)
			F2_nk[i*n2 + j] =  G[12];
			G2_nk[i*n2 + j] =  G[13];

      F2q_nk[i*n2 + j] = G[14];
      G2q_nk[i*n2 + j] = G[15];

      //3rd order used in P13 : F3/G3[k,p,-p]
      F3_nk[i*n2 + j] = G[24];
      G3_nk[i*n2 + j] = G[25];

      F3q_nk[i*n2 + j] = G[26];
      G3q_nk[i*n2 + j] = G[27];

			gsl_odeiv2_driver_free(d);
		}
	}

  //  Solution for wCDM linear growth @ a=1 for our normalisation of the linear power spectrum
    	double a = 0.0001;
    	double G[2] = {a,-a}; // initial conditions
    	struct param_type2 params_my = {omega0,par1,par2,0.,1};

    // Solutions of evolution factors @ a=A
    	gsl_odeiv2_system sys3 = {funcn3, jac, 2, &params_my};
    	gsl_odeiv2_driver * d1 =
    	gsl_odeiv2_driver_alloc_y_new (&sys3, gsl_odeiv2_step_rk8pd,
    								   1e-6, 1e-6, 0.0);

    	int status2 = gsl_odeiv2_driver_apply (d1, &a, 1. , G);

    	dnorm_spt = G[0];
      gnorm_spt = G[1];


    	gsl_odeiv2_driver_free(d1);

}



int funcn4_ginvlin(double a, const double G[], double F[], void *params)
{
  param_type3 p = *(param_type3 *)(params);
	real k1 = p.kk1;
	real x1= p.xx1;
	real k2 = p.kk2;
	real x2= p.xx2;
	real k3 = p.kk3;
	real x3= p.xx3;
	real omega0 = p.omega00;
  // p1 = cs2, p2 = w or mg params
  // p3 chooses cde(1) or mg (0)
	real p1 = p.par1;
	real p2 = p.par2;
  real p3 = p.par3;
	int choice = p.arg1;

if (choice==0) {
  double hub = HA2(a,omega0);

	/* 1st order */
	//1. F1/G1(k)
	F[0] = -G[1]/a;
	F[1] =1./a*(-(2.-hub)*G[1]-hub*G[0]*mu(a,k1,omega0,p1,p2,p3));
  F[2] = 0.;
  F[3] = 0.;
}
else{
// hubble functions and omega_m
  double hade1,hade2,hade3,hade4,A,omegaf,omegaL,h0,c1;
  h0 = 1./2997.9;
  hade1 = HA2de(a,omega0,p2);
  hade2 = HA2de2(a,omega0,p2);
  A = -3.*(1.+p2);
  omegaf = pow(a,A);
  omegaL= (1.-omega0)*omegaf;
  hade3 = pow3(a)*omegaL*hade2/omega0;
  hade4 = pow2(HAde(a,omega0,p2)*h0*a);
  c1 = p1/(hade4*(1.+1e-8+p2));

	//1. F1/G1(k)
	F[0] = -G[1]/a;
	F[1] =1./a*(-(2.-hade1)*G[1]-hade2*G[0]-hade3*G[2]);
  F[2] =1./a*(-(1+p2)*G[3]+3*(p2-p1)*G[2] + 9.*(1+p2)*(p2-p1)*hade4/pow2(k1)*G[3]);
  F[3] =1./a*(-(2.-3*p1-hade1)*G[3]-hade2*G[0]-(hade3 - c1*pow2(k1))*G[2]);

}
	return GSL_SUCCESS;
}


// Include MG extension
void IOW::initn_lin(int choice, double A, double k, double omega0, double par1, double par2, double par3)
{
				double a = 0.0001;

// Non-Eds ICs
			  double G[4] = { a,-a,a,-a};
        double arg1 = choice;

			/*Parameters passed to system of equations */
// EDIT : If more than one gravity parameter is used, add them after p1
				struct param_type3 my_params1 = {k,1.,1.,1.,1.,1.,omega0, par1, par2, par3,arg1,1.,1.,1.};

				gsl_odeiv2_system sys = {funcn4_ginvlin, jac, 4, &my_params1};

			//  this gives the system, the method, the initial interval step, the absolute error and the relative error.
				gsl_odeiv2_driver * d =
				gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
										   1e-4, 1e-3, 1e-3);

				int status1 = gsl_odeiv2_driver_apply (d, &a, A, G);

				/*Allocation of array values */

			//F1(k;a), G1(k;a)
			F1_nk = G[0] ;
			G1_nk = G[1] ;
      F1q_nk = G[2] ;
      G1q_nk = G[3] ;

			gsl_odeiv2_driver_free(d);

}


// Bilinear search in x128. Used in :
// A13 MC integral
// AB term for AT limit location
int searchnearest(const double anArray[], double key){
	    int low = 128;
	    int high = 255;
			double keym = fabs(key);
	    while (low < high) {
	        int mid = (low + high) / 2;
	        double d1 = fabs(anArray[mid] - keym);
	        double d2 = fabs(anArray[mid+1] - keym);
	        if (d2 <= d1)
	        {
	            low = mid+1;
	        }
	        else
	        {
	            high = mid;
	        }
	    }
			if(key<0.) {
				return 255-high;
			}
			if (keym<1e-7) {
				return 127;
			}
			else{
				return high;
			}
	}



  /*Selection functions for numerical functions */

int num_selec_maga(double ymin, double ymax, double y){
  	int mag_int;
  //	 mag_int = (n2*1.-1)*(y-ymin)/ymax + 1./n2; // linear
  //	 mag_int = sqrt((y-ymin)/ymax)*(n2-1); //quadratic
  		mag_int = (int)round((n2-1)*log(y/ymin)/log(ymax/ymin)); //exponential
  		return mag_int;
  }

/* Function that gives the standalone numerical kernels */
// a*b = 1,2,3,4,5,6 : F1,G1,F2,G2,F3,G3 respectively.
real IOW::Kernel(double r, double x, double rmin, double rmax,  int a) const{
	int y = num_selec_maga(rmin, rmax, r);
	int z = searchnearest(x128,x);
	switch (a) {
		case 1:
			return	F1_nk/dnorm_spt;// F1;
		case 2:
			return   dnorm_spt;
		//F2/G2(p,k-p) (P22)
		case 3:
			return	 fdgp_spt;// F2[z*n2+y];
		case 4:
			return	 G2_nk[z*n2 + y];
		case 5:
			return	 F3_nk[z*n2 + y];
		case 6:
			return	 G3_nk[z*n2 + y];
			//F2/G2(p,-k)
		case 7:
			return   F2D_nk[z*n2 + y];//F2B_nk[z*n2 + y];
		case 8:
			return	 G2B_nk[z*n2 + y];
	}
}



/*Kernels for multipole factors*/

/* multipole factors : Example for TNS (case 6) (2l+1)/2 * integ dmu  e^(-(k*f*sigma_v*mu)^2) * mu^(2n) * P_l(mu)*/

static double kernel1(real k, real sigma_v, real F0, real anw, int n, int a, real u1) {
switch(a) {
  case 1:
	return pow(u1,2*n) * (1.-sigma_v*pow2(u1*F0));//RESUM
  case 2:
  return pow(u1,2*n) * (1.-sigma_v*pow2(u1*F0)) * exp(-0.5*pow2(k)*anw*(1+F0*(2+F0)*pow2(u1)));//RESUM
  case 3:
  return pow(u1,2*n) * 0.5 * pow2(k)*anw*(1+F0*(2+F0)*pow2(u1)) * exp(-0.5*pow2(k)*anw*(1+F0*(2+F0)*pow2(u1))); //RESUM
  case 4:
  return pow(u1,2*n) * pow2(1+pow2(u1)*F0); // MATSUBARA
  case 5:
  return pow(u1,2*n) * (1+pow2(u1)*F0); // MATSUBARA
  case 6:
  return pow(u1,2*n) *exp(-k*k*u1*u1*sigma_v*sigma_v) ; // TNS Gaussian
  case 7:
  return pow(u1,2*n) * 1./(1.+(k*k*u1*u1*sigma_v*sigma_v)/2.)  ; // TNS Lorentzian
}
}

static double kernel2(real k, real sigma_v, real F0, real anw, int n, int a, real u1) {
   switch(a) {
     case 1:
     return 3.*pow(u1,2*n+2)*(1.-sigma_v*pow2(u1*F0)) - pow(u1,2*n)*(1.-sigma_v*pow2(u1*F0));
     case 2:
     return (3.*pow(u1,2*n+2)*(1.-sigma_v*pow2(u1*F0)) - pow(u1,2*n)*(1.-sigma_v*pow2(u1*F0)))* exp(-0.5*pow2(k)*anw*(1.+F0*(2+F0)*pow2(u1)));
     case 3:
     return (3.*pow(u1,2*n+2) - pow(u1,2*n))* 0.5 * pow2(k)*anw*(1.+F0*(2+F0)*pow2(u1)) * exp(-0.5*pow2(k)*anw*(1.+F0*(2+F0)*pow2(u1)));
     case 4:
     return 3.*pow(u1,2*n+2)* pow2(1.+pow2(u1)*F0) - pow(u1,2*n)* pow2(1.+pow2(u1)*F0);
     case 5:
     return 3.*pow(u1,2*n+2)* (1.+pow2(u1)*F0) - pow(u1,2*n)* (1.+pow2(u1)*F0);
     case 6:
     return (3.*pow(u1,2*n+2)*exp(-k*k*u1*u1*sigma_v*sigma_v) - pow(u1,2*n)*exp(-k*k*u1*u1*sigma_v*sigma_v));
     case 7:
     return (3.*pow(u1,2*n+2)* 1./(1.+(k*k*u1*u1*sigma_v*sigma_v)/2.)  - pow(u1,2*n)* 1./(1.+(k*k*u1*u1*sigma_v*sigma_v)/2.));
        }
}

static double kernel3(real k, real sigma_v, real F0, real anw, int n, int a, real u1) {
  switch(a){
    case 1:
    return  (35*pow(u1,4)-30*u1*u1+3)*pow(u1,2*n)*(1-sigma_v*pow2(u1*F0));
    case 2:
    return  ((35*pow(u1,4)-30*u1*u1+3)*pow(u1,2*n)*(1-sigma_v*pow2(u1*F0)))* exp(-0.5*pow2(k)*anw*(1+F0*(2+F0)*pow2(u1)));
    case 3:
    return  ((35*pow(u1,4)-30*u1*u1+3)*pow(u1,2*n))*0.5*pow2(k)*anw*(1+F0*(2+F0)*pow2(u1))*exp(-0.5*pow2(k)*anw*(1+F0*(2+F0)*pow2(u1)));
    case 4:
    return  ((35*pow(u1,4)-30*u1*u1+3)*pow(u1,2*n)) * pow2(1+pow2(u1)*F0);
    case 5:
    return  ((35*pow(u1,4)-30*u1*u1+3)*pow(u1,2*n)) * (1+pow2(u1)*F0);
    case 6:
    return  (35*pow(u1,4)-30*u1*u1+3)*pow(u1,2*n)*exp(-k*k*u1*u1*sigma_v*sigma_v);
    case 7:
    return  (35*pow(u1,4)-30*u1*u1+3)*pow(u1,2*n)* 1./(1.+(k*k*u1*u1*sigma_v*sigma_v)/2.);
}
}

// a =1 : monopole
// a =2 : quadrupole
// a =3 : hexdecapole
double factL(real k, real sigma_v, real F0, real anw, int n, int e, int a ){
  switch(e) {
        case 1:
			       return 0.5*Integrate(bind(kernel1,k, sigma_v, F0, anw, n, a, _1), -1. , 1., 1e-3);
        case 2:
            return 5./4.*Integrate(bind(kernel2, k, sigma_v, F0, anw, n, a, _1), -1. , 1., 1e-3);
        case 3:
            return 9./16.*Integrate(bind(kernel3, k, sigma_v, F0, anw, n, a, _1), -1. , 1., 1e-3);
	default:
            warning("FACTL: invalid indices, a = %d \n", a);
            return 0;
    }
}


/* Fingers of god term in exponential form */

// a*b =1 : LCDM
// a*b =2 : nDGP
// a*b =3 : numerical (arbitrary model)

double DFOG(real k, real u,  double sigma_v, int a) {
	double F0;
	switch (a) {
		case 1:
			F0= fl_spt;
			break;
		case 2:
			F0= fdgp_spt;
			break;
		case 3:
			F0= 1;//-G1_nk/F1_nk;
			break;
	}
	return exp(-k*k*u*u*F0*F0*sigma_v*sigma_v);
}






/* GSM EQs.A7 and Squares of A8 appearing in A13 */
// splitting of Eq.7 - used in precomputing components of sig_12^2
double mybessel1(double kr){
	if (kr<0.3) {
	return (sin(kr)/kr-3.*(1./3.-pow2(kr)/30.+pow4(kr)/840.)) ;
	}
	return SphericalBesselJ0(kr)-3.*SphericalBesselJ1(kr)/kr;
}

double mybessel2(double kr){
	if (kr<0.3) {
	return (1./3.-pow2(kr)/30.+pow4(kr)/840.);
	}
	return SphericalBesselJ1(kr)/kr;
}


double J1(double u, double kr){
	if (kr<0.3) {
	return u*(sin(kr)/kr-2.*(1./3.-pow2(kr)/30.+pow4(kr)/840.))+(1-u)*(1./3.-pow2(kr)/30.+pow4(kr)/840.);
	}
	else
	return u*(SphericalBesselJ0(kr)-2.*SphericalBesselJ1(kr)/kr)+(1-u)*SphericalBesselJ1(kr)/kr;
}

double I11(double y1){
	return  2.*SphericalBesselJ1(y1);
}

double I12(double y1){
	if (y1<0.3) {
		return 2.*(sin(y1)/y1 - 2.*(1./3.-pow2(y1)/30.+pow4(y1)/840.));
	}
	return  2.*(2.*y1*cos(y1)+(-2.+y1*y1)*sin(y1))/pow3(y1);
}

double I30(double y1){
	if (y1<0.3) {
		return 4*(1./3.-pow2(y1)/30.+pow4(y1)/840.);
	}
	  return 4.*(sin(y1) - y1*cos(y1))/pow3(y1); //
}

double I31(double y1){
	if (y1<0.3) {
		return 4*(y1*(1./15.-pow2(y1)/210.+pow4(y1)/7560.));
	}
	return -4.*(3*y1*cos(y1)+(-3.+y1*y1)*sin(y1))/pow4(y1);
}

double I13(double y1){
	if (y1<0.3) {
  return 2*(y1*(1./5.- pow2(y1)/42.+pow4(y1)/1080.));
		}
	return  2*(y1*(6.-y1*y1)*cos(y1)+3.*(-2.+y1*y1)*sin(y1))/pow4(y1);
	}



/* Use standard library implementations */
double BesselJ0(double x) { return j0(x); }
double BesselJ1(double x) { return j1(x); }
double BesselJn(int n, double x) { return jn(n, x); }

double SphericalBesselJ0(double x) {
    if(fabs(x) < 1e-4)
        return 1 - pow2(x)/6. + pow4(x)/120. - pow6(x)/5040.;
    else
        return sin(x)/x;
}

double SphericalBesselJ1(double x) {
    if(fabs(x) < 1e-4)
        return x/3. - pow3(x)/30. + pow5(x)/840. - pow7(x)/45360.;
    else
        return (sin(x) - x*cos(x))/pow2(x);
}

double SphericalBesselJ2(double x) {
    if(fabs(x) < 1e-4)
        return pow2(x)/15. - pow4(x)/210. + pow6(x)/7560.;
    else
        return ((3. - pow2(x))*sin(x) - 3.*x*cos(x))/pow3(x);
}

double SphericalBesselJ3(double x) {
    if(fabs(x) < 1e-4)
        return pow3(x)/105. - pow5(x)/1890. + pow7(x)/83160.;
    else
        return ((15. - 6.*pow2(x))*sin(x) - (15.*x - pow3(x))*cos(x))/pow4(x);
}

double SphericalBesselJ4(double x) {
    if(fabs(x) < 1e-4)
        return pow4(x)/945. - pow6(x)/20790.;
    else
        return ((105. - 45.*pow2(x) + pow4(x))*sin(x) - (105.*x + 10.*pow3(x))*cos(x))/pow5(x);
}

///// Bessel Roots for FT of Power Spectra ///////
int static besselsetup = 0;

int fillbesselroots() {
  if(besselsetup == 1) {
    return 0;
    }
  //get spherical bessel function roots.
  double xval;
  double delta = 0.001;
  double fval;
  double fvallast;
  int ello;
  int whichroot;
  double currroot;

  for(ello=0;ello<=2;ello+=2) {
    whichroot=0;
    fvallast=gsl_sf_bessel_jl(ello,0.0); // Bessel Function of ello^th order and argument 0.
    xval = delta;
    currroot = 0.0;
    while(currroot < 2500.0) { // fill up to 1500 roots up to arg=2500
      assert(whichroot < BMAX);
      fval=gsl_sf_bessel_jl(ello,xval); // evaluate bessel at 0.001 to start
      if((fval > 0.0 && fvallast < 0.0) || (fval < 0.0 && fvallast > 0.0)) {
        besselroots[whichroot][ello] = xval - 0.5*delta;
        currroot = xval - 0.5*delta;
        whichroot += 1;
        if(whichroot >= BMAX) {
          printf("ahh bessel error! %e\n",currroot);
          exit(1);
          }
        }
      xval += delta;
      fvallast = fval;
      }
    assert(whichroot < BMAX-1);
    }
    besselsetup = 1;
  return 0;
  }


double getbesselroot(int whichroot, int ell) {
  assert(whichroot <= BMAX-1);
  return(besselroots[whichroot][ell]);
  }


  // wigner symbol under sterling approx (see A3 of arxiv 0310125) for CMB bispectrum
 double wigner3j(int l1, int l2, int l3){
    int L1 = l1 + l2 + l3;
    if(L1 % 2 == 0 ) {
    double L = (double)L1/2;
    double larg = 1./((L*1.+1.) * (L*1.-l1*1.+1.) * (L*1.-l2*1.+1.) * (L*1.-l3*1.+1.)) ;
    return pow(-1,L)*sqrt(exp(3.)/2./M_PI) * pow(larg,0.25) * pow((L*1. - l1*1. + 0.5)/(L*1. - l1*1. + 1.), L*1. - l1*1. + 0.25)
                                                             * pow((L*1. - l2*1. + 0.5)/(L*1. - l2*1. + 1.), L*1. - l2*1. + 0.25)
                                                             * pow((L*1. - l3*1. + 0.5)/(L*1. - l3*1. + 1.), L*1. - l3*1. + 0.25);

                    }
    else {
    return 0.;}
  }




////////////////////////////////////////////////////////////
// Static variables

////////////////////////////////////////////////////////////
// Constructor, desctructor and initializer

integral::integral(  )
{
    this->gl_intg_res = 0.;
    clear(  );
}

integral::~integral(  )
{

}

void integral::clear(  )
{
    x_buf.clear(  );
    y_buf.clear(  );
    return;
}

void integral::clear( const unsigned & size )
{
    x_buf.resize( size );
    y_buf.resize( size );
    return;
}

////////////////////////////////////////////////////////////
// Read data

void integral::read( const double & x, const double & y )
{
    x_buf.push_back( x );
    y_buf.push_back( y );
    return;
}

////////////////////////////////////////////////////////////
// Integrate--just in case that the x are not evenly spaced

double integral::result(  )
{
    double res( 0. );
#pragma omp parallel for reduction( + : res )
    for( unsigned i = 1; i < x_buf.size(  ); ++ i )
	res += 0.5 * ( x_buf.at( i ) - x_buf.at( i - 1 ) )
		   * ( y_buf.at( i ) + y_buf.at( i - 1 ) );
    return res;
}

////////////////////////////////////////////////////////////
// Gauss-Legendre

void integral::gl_clear(  )
{
    this->gl_intg_res = 0.;
    return;
}

double integral::gl_xi( const int & i )
{
    return x[ i ];
}

void integral::gl_read( const int & i, const double & kernel )
{
    gl_intg_res += w[ i ] * kernel;
    return;
}

double integral::gl_result(  )
{
    return gl_intg_res;
}



// Legendre polynomials for GSM multipoles
double legendre(int order, double x ){
    switch( order )
    {
    case 0:
	return 1.;
    case 2:
	return 0.5 * ( 3. * pow( x, 2 ) - 1. );
    case 4:
	return ( 35. * pow( x, 4 ) - 30. * pow( x, 2 )
		 + 3. ) / 8.;
     default:
      warning("Legendre: invalid indices");
      return 0;
      }
}

///// initialization for old
//Clustering DE : SEE 1611.00542 for expressions used here
// int funcn4(double a, const double G[], double F[], void *params)

// int IOW::initn_cde(double A, double YMIN, double YMAX, double k, double omega0, double par1, double par2, double par3)
// {
// //#pragma omp parallel for schedule(dynamic)
// 	for ( int i = 0; i < n1; ++i)
// 	{
// 		for (int j = 0; j< n2; ++j)
// 		{
// 				real x1,x2,x3,k1,k2,k3,karg1,karg2,karg3,karg4;
//
// 				/* k - magnitude sampling */
// 				double R = (YMIN) * exp(j*log(YMAX/YMIN)/(n2*1.-1.)); //exponential sampling
//
// 				/*Parameter selection for kernels */
//
// 				k1 = k;
// 		    k2 = k*R;
// 				k3 = k2;
//
// 				x1 = XMIN;
//         x2 = x128[i]; // x128 goes from - 0.99995 to 0.99995
// 				x3 = -x2;
//
//         karg1 = sqrt(k2*k2+k1*k1-2*k2*k1*x2);
//
// 				double a = 0.0001;
//
//
// // Non-Eds ICs
// 			  double G[24] = { a,-a,a,-a,a,-a,a,-a,a,-a,a,-a, 0.,0.,0.,0., 0.,0.,0.,0., 0.,0.,0.,0.};
//
// 			/*Parameters passed to system of equations */
// // EDIT : If more than one gravity parameter is used, add them after p1
// 				struct param_type3 my_params1 = {k1,x1,k2,x2,k3,x3,omega0, par1, par2, par3,karg1,1.,1.,1.};
//
// 				gsl_odeiv2_system sys = {funcn4, jac, 24, &my_params1};
//
// 			//  this gives the system, the method, the initial interval step, the absolute error and the relative error.
// 				gsl_odeiv2_driver * d =
// 				gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
// 										   1e-3, 1e-3, 1e-3);
//
// 				int status1 = gsl_odeiv2_driver_apply (d, &a, A, G);
//
//
// 				/*Allocation of array values */
//
// 			//F1(k;a), G1(k;a)
// 			F1_nk = G[0] ;
// 			G1_nk = G[1] ;
//
//       //F1q(k;a), G1q(k;a)
//       F1q_nk = G[2];
//       G1q_nk = G[3];
//
// 			// F1(k-p;a), G1(k-p;a)
// 			F1kmp_nk[i*n2 + j] = G[4];
// 			G1kmp_nk[i*n2 + j] = G[5];
//
// 			//F1(p;a), G1(p;a)
// 			F1p_nk[i*n2 + j] = G[8];
// 			G1p_nk[i*n2 + j] = G[9];
//
//       //F1q(p;a), G1q(p;a)
//       F1pq_nk[i*n2 + j] = G[10];
//       G1pq_nk[i*n2 + j] = G[11];
//
// 			/*2nd order*/
//
// 			//F2/G2(p,k-p) (P22)
// 			F2_nk[i*n2 + j] =  G[12];
// 			G2_nk[i*n2 + j] =  G[13];
//
//
//       //F2/G2(p,p)
//       F2D_nk[i*n2 + j] =  G[16];
//       G2D_nk[i*n2 + j] =  G[17];
//
//       F2Dq_nk[i*n2 + j] =  G[18];
//       G2Dq_nk[i*n2 + j] =  G[19];
//
//       //F2/G2(p,k)
//       F2A_nk[i*n2 + j] =  G[20];
//       G2A_nk[i*n2 + j] =  G[21];
//
//
// 			gsl_odeiv2_driver_free(d);
// 		}
// 	}
//
//   //  Solution for wCDM linear growth @ a=1 for our normalisation of the linear power spectrum
//     	double a = 0.0001;
//     	double G[2] = {a,-a}; // initial conditions
//     	struct param_type2 params_my = {omega0,par1,par2,0.,1};
//
//     // Solutions of evolution factors @ a=A
//     	gsl_odeiv2_system sys3 = {funcn3, jac, 2, &params_my};
//     	gsl_odeiv2_driver * d1 =
//     	gsl_odeiv2_driver_alloc_y_new (&sys3, gsl_odeiv2_step_rk8pd,
//     								   1e-3, 1e-3, 0.0);
//
//     	int status2 = gsl_odeiv2_driver_apply (d1, &a, 1. , G);
//
//     	dnorm_spt = G[0];
//       gnorm_spt = G[1];
//
//
//     	gsl_odeiv2_driver_free(d1);
//
// 	return 0;
// }






// //  Clustering DE : SEE 1611.00542
// int funcn4(double a, const double G[], double F[], void *params)
// {
//   param_type3 p = *(param_type3 *)(params);
// 	real k1 = p.kk1;
// 	real x1= p.xx1;
// 	real k2 = p.kk2;
// 	real x2= p.xx2;
// 	real k3 = p.kk3;
// 	real x3= p.xx3;
// 	real omega0 = p.omega00;
//   // p1 = cs2, p2 = w
// 	real p1 = p.par1;
// 	real p2 = p.par2;
// 	real p3 = p.par3;
//   real karg1 = p.arg1;
//   real karg2 = p.arg2;
//   real karg3 = p.arg3;
//   real karg4 = p.arg4;
//
//
// // hubble functions and omega_m
//   double hade1,hade2,hade3,hade4,A,omegaf,omegaL,h0;
//   h0 = 1./2997.9;
//   hade1 = HA2de(a,omega0,p2);
//   hade2 = HA2de2(a,omega0,p2);
//   A = -3.*(1.+p2);
//   omegaf = pow(a,A);
//   omegaL= (1.-omega0)*omegaf;
//   hade3 = pow3(a)*omegaL*hade2/omega0;
//   hade4 = HAde(a,omega0,p2)*h0;
//
//   double a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14;
//   double b1,b2,b3,b4,b5,b6,b7;
//
// // different k arguments
//
//   a1 = alpha(k2,karg1,(k1*x2-k2)/karg1);
//   a2 = alpha(karg1,k2,(k1*x2-k2)/karg1);
//   a3 = alpha(k1,karg1,(k2*x2-k1)/karg1);
//   a4 = alpha(karg1,k1,(k2*x2-k1)/karg1);
//   a5 = alpha(k2,k1,x2);
//   a6 = alpha(k1,k2,x2);
//   a7 = alpha(k3,k1,x3);
//   a8 = alpha(k1,k3,x3);
//
//   a9 = alpha(k1,karg3,(k2*x2+k3*x3)/karg3);
//   a10 = alpha(k3,karg4,(k1*x3+k2*x1)/karg4);
//   a11 = alpha(k2,karg2,(k3*x1+k1*x2)/karg2);
//   a12 = alpha(karg3,k1,(k2*x2+k3*x3)/karg3);
//   a13 = alpha(karg4,k3,(k1*x3+k2*x1)/karg4);
//   a14 = alpha(karg2,k2,(k3*x1+k1*x2)/karg2);
//
//   b1 = beta1(k2,karg1,(k1*x2-k2)/karg1);
//   b2 = beta1(k1,karg1,(k2*x2-k1)/karg1);
//   b3 = beta1(k1,k2,x2);
//   b4 = beta1(k1,k3,x3);
//   b5 = beta1(k1,karg3,(k2*x2+k3*x3)/karg3);
//   b6 = beta1(k3,karg4,(k1*x3+k2*x1)/karg4);
//   b7 = beta1(k2,karg2,(k3*x1+k1*x2)/karg2);
//
//
// 	/* 1st order */
// 	//1. F1/G1(k)
// 	F[0] = -G[1]/a;
// 	F[1] =1./a*(-(2.-hade1)*G[1]-hade2*G[0]-hade3*G[2]);
//   F[2] =1./a*(-(1+p2)*G[3]+3*(p2-p1)*G[2]);
//   F[3] =1./a*(-(2.-3*p2-hade1)*G[3]-hade2*G[0]-(hade3*(1.+3.*p1) + p1*pow2(k1)/(pow2(hade4*a)*(1.000001+p2)))*G[2]);
//
//
// 	//2. F1/G1(k-p)
//   F[4] = -G[5]/a;
// 	F[5] =1./a*(-(2.-hade1)*G[5]-hade2*G[4]-hade3*G[6]);
//   F[6] =1./a*(-(1+p2)*G[7]+3*(p2-p1)*G[6]);
//   F[7] =1./a*(-(2.-3*p2-hade1)*G[7]-hade2*G[4]-(hade3*(1.+3.*p1) + 0*p1*(k2*k2+k1*k1-2*k2*k1*x2)/(pow2(hade4*a)*(1.000001+p2)))*G[6]);
//
// 	//3. F1/G1(p)
//   F[8] = -G[9]/a;
// 	F[9] =1./a*(-(2.-hade1)*G[9]-hade2*G[8]-hade3*G[10]);
//   F[10] =1./a*(-(1+p2)*G[11]+3*(p2-p1)*G[10]);
//   F[11] =1./a*(-(2.-3*p2-hade1)*G[11]-hade2*G[8]-(hade3*(1.+3.*p1) + 0*p1*pow2(k2)/(pow2(hade4*a)*(1.000001+p2)))*G[10]);
//
// 	/* 2nd order */
//   //4. F2/G2(p,k-p) (P22)
//   F[12] =1/a*(-(a1*G[9]*G[4]+a2*G[5]*G[8])/2.-G[13]);
// 	F[13] =1/a*(-(2.-hade1)*G[13]-hade2*G[12]-hade3*G[14]-b1*G[9]*G[5]);
//   F[14] =1/a*(-(1+p1)*(a1*G[11]*G[6]+a2*G[7]*G[10])/2.-(1+p2)*G[15]+3*(p2-p1)*G[14]);
// 	F[15] =1/a*(-(2.-3*p2-hade1)*G[15]-hade2*G[12]-(hade3*(1.+3.*p1) + 0*p1*pow2(k1)/(pow2(hade4*a)*(1.000001+p2)))*G[14]- b1*G[11]*G[7]);
//
//
//   //5. F2/G2(-k,k-p)
//   F[16] =1/a*(-(a3*G[1]*G[4]+a4*G[5]*G[0])/2.-G[17]);
// 	F[17] =1/a*(-(2.-hade1)*G[17]-hade2*G[16]-hade3*G[18] - b2*G[5]*G[1]);
//   F[18] =1/a*(-(1+p1)*(a3*G[3]*G[6]+a4*G[7]*G[2])/2.-(1+p2)*G[19]+3*(p2-p1)*G[18]);
//   F[19] =1/a*(-(2.- 3*p2 -hade1)*G[19]-hade2*G[16]-(hade3*(1.+3.*p1) + 0*p1*pow2(k2)/(pow2(hade4*a)*(1.000001+p2)))*G[18] - b2*G[7]*G[3]);
//
//   //6. F2/G2(p,-p)
// 	F[20] =1./a*(-G[21]) ; // set 1st term to 0 because x1 is not identically -1 (needed to keep evolution at 0 for GR)
// 	F[21] =1./a*(-(2.-hade1)*G[21]-hade2*G[20]-hade3*G[22]);
//   F[22] =1/a*(-(1+p2)*G[23]+3*(p2-p1)*G[22]);
// 	F[23] =1/a*(-(2.-3*p2-hade1)*G[23]-hade2*G[20]-hade3*G[22]*(1.+3.*p1));
//
//
//   //7. F2/G2(p,k)
// 	F[24] =1/a*(-(a5*G[9]*G[0]+a6*G[1]*G[8])/2.-G[25]) ;
// 	F[25] =1/a*(-(2.-hade1)*G[25]-hade2*G[24]-hade3*G[26]-b3*G[9]*G[1]);
//   F[26] =1/a*(-(1+p1)*(a5*G[11]*G[2]+a6*G[3]*G[10])/2.-(1+p2)*G[27]+3*(p2-p1)*G[26]) ;
// 	F[27] =1/a*(-(2.-3*p2-hade1)*G[27]-hade2*G[24]-(hade3*(1.+3.*p1) + 0*p1*(k2*k2+k1*k1+2*k2*k1*x2)/(pow2(hade4*a)*(1.000001+p2)))*G[26]-b3*G[11]*G[3]);
//
//
//   //8. F2/G2(-p,k)=F2/G2(p,-k)
//   F[28] =1/a*(-(a7*G[9]*G[0]+a8*G[1]*G[8])/2.-G[29]) ;
//   F[29] =1/a*(-(2.-hade1)*G[29]-hade2*G[28]-hade3*G[30]- b4*G[9]*G[1]);
//   F[30] =1/a*(-(1+p1)*(a7*G[11]*G[2]+a8*G[3]*G[10])/2.-(1+p2)*G[31]+3*(p2-p1)*G[30]) ;
//   F[31] =1/a*(-(2.-3*p2-hade1)*G[31]-hade2*G[28]-(hade3*(1.+3.*p1) + 0*p1*(k2*k2+k1*k1-2*k2*k1*x2)/(pow2(hade4*a)*(1.000001+p2)))*G[30]- b4*G[11]*G[3]);
//
//
//   //9. 3rd order  ~ F3/G3(k,p,-p)
// 	F[32] = - 1/(3.*a)*(a9*G[20]*G[1]
// 	          + a10*G[24]*G[9]
// 						+ a11*G[28]*G[9]
// 						+ a13*G[25]*G[8]
// 						+ a14*G[29]*G[8]
// 						+ a12*G[21]*G[0]
// 						+3.*G[33]) ;
//
//
// 	F[33] =1/(3.*a)*(-3.*(2.-hade1)*G[33]-3.*hade2*G[32] - 3*hade3*G[34]-
// 					 -2*b5*G[1]*G[21]
// 					 -2*b6*G[9]*G[25]
// 					 -2*b7*G[9]*G[29]
//               );
//
//   F[34] = - 1/(3.*a)*((1+p1)*(a10*G[22]*G[3]
//               +  a11*G[26]*G[11]
//               +  a12*G[30]*G[11]
//               + a13*G[27]*G[10]
//               + a14*G[31]*G[10]
//               + a12*G[23]*G[2])
//               -(1+p2)*G[35]+3*(p2-p1)*G[34]) ;
//
//
//   F[35] =1/(3.*a)*(-3.*(2.-3*p2-hade1)*G[35]-3.*hade2*G[32] - 3*(hade3*(1.+3.*p1) + p1*pow2(k1)/(pow2(hade4*a)*(1.000001+p2)))*G[34]-
//               -2*b5*G[3]*G[23]
//               -2*b6*G[11]*G[27]
//               -2*b7*G[11]*G[31]
//             );
//
// 	return GSL_SUCCESS;
// }



//  Clustering DE : SEE 1611.00542 for expressions used here
// int funcn4(double a, const double G[], double F[], void *params)
// {
//   param_type3 p = *(param_type3 *)(params);
// 	real k1 = p.kk1;
// 	real x1= p.xx1;
// 	real k2 = p.kk2;
// 	real x2= p.xx2;
// 	real k3 = p.kk3;
// 	real x3= p.xx3;
// 	real omega0 = p.omega00;
//   // p1 = cs2, p2 = w
// 	real p1 = p.par1;
// 	real p2 = p.par2;
// 	real p3 = p.par3;
//   real karg1 = p.arg1;
//
//
// // hubble functions and omega_m
//   double hade1,hade2,hade3,hade4,A,omegaf,omegaL,h0;
//   h0 = 1./2997.9;
//   hade1 = HA2de(a,omega0,p2);
//   hade2 = HA2de2(a,omega0,p2);
//   A = -3.*(1.+p2);
//   omegaf = pow(a,A);
//   omegaL= (1.-omega0)*omegaf;
//   hade3 = p3*pow3(a)*omegaL*hade2/omega0;
//   hade4 = pow2(HAde(a,omega0,p2)*h0*a);
//
//   double a1,a2,a3,a4,a5,a6;
//   double b1,b2,b3;
//   double c1;
//   c1 = p1/(hade4*(1.+1e-8+p2));
//
// // different k arguments
//
//   a1 = alpha(k2,karg1,(-k1*x2-k2)/karg1);
//   a2 = alpha(karg1,k2,(-k1*x2-k2)/karg1);
//   b1 = beta1(k2,karg1,(-k1*x2-k2)/karg1);
//
//   a3 = alpha(k1,karg1,(-k2*x2-k1)/karg1);
//   a4 = alpha(karg1,k1,(-k2*x2-k1)/karg1);
//   b2 = beta1(k1,karg1,(-k2*x2-k1)/karg1);
//
//   a5 = alpha(k2,k1,x2);
//   a6 = alpha(k1,k2,x2);
//   b3 = beta1(k1,k2,x2);
//
//   //  a3 = alpha(k2,k2,x4);
//   //  b4 = beta1(k2,k2,x4);
//
// 	/* 1st order */
// 	//1. F1/G1(k)
// 	F[0] = -G[1]/a;
// 	F[1] =1./a*(-(2.-hade1)*G[1]-hade2*G[0]-hade3*(1.+3.*p1)*G[2]);
//   F[2] =1./a*(-(1+p2)*G[3]+3*(p2-p1)*G[2]);
//   F[3] =1./a*(-(2.-3*p2-hade1)*G[3]-hade2*G[0]-(hade3*(1.+3.*p1) - c1*pow2(k1))*G[2]);
//
//   //2. F1/G1(-k-p)
//   F[4] = -G[5]/a;
// 	F[5] =1./a*(-(2.-hade1)*G[5]-hade2*G[4]-hade3*(1.+3.*p1)*G[6]);
//   F[6] =1./a*(-(1+p2)*G[7]+3*(p2-p1)*G[6]);
//   F[7] =1./a*(-(2.-3*p2-hade1)*G[7]-hade2*G[4]-(hade3*(1.+3.*p1) - c1*pow2(karg1))*G[6]);
//
// 	//3. F1/G1(p)
//   F[8] = -G[9]/a;
// 	F[9] =1./a*(-(2.-hade1)*G[9]-hade2*G[8]-hade3*(1.+3.*p1)*G[10]);
//   F[10] =1./a*(-(1+p2)*G[11]+3*(p2-p1)*G[10]);
//   F[11] =1./a*(-(2.-3*p2-hade1)*G[11]-hade2*G[8]-(hade3*(1.+3.*p1) - c1*pow2(k2))*G[10]);
//
// /* 2nd order */
//   //4. F2/G2(p,-k-p) (P22)
//   F[12] =1./a*(-(a1*G[9]*G[4]+a2*G[5]*G[8])/2.-G[13]);
// 	F[13] =1./a*(-(2.-hade1)*G[13]-hade2*G[12]-hade3*(1.+3.*p1)*G[14]-b1*G[9]*G[5]);
//
//
//   F[14] =1./a*(-(1+p1)*(a1*G[11]*G[6]+a2*G[7]*G[10])/2.-(1+p2)*G[15]+3*(p2-p1)*G[14]);
// 	F[15] =1./a*(-(2.-3*p2-hade1)*G[15]-hade2*G[12]-(hade3*(1.+3.*p1) - c1*pow2(k1))*G[14]- b1*G[11]*G[7]
//                   +3.*(p1-p2)/(1.+1e-8+p2)*(a1*G[11]*G[6]+a2*G[7]*G[10])/2. - p1*a/(1.+1e-8+p2)*(a1*G[11]*F[6]+a2*G[7]*F[10])/2.
//                   -4.*(p1+1)*c1/(1.+1e-8+p2)*(-4.*pow2(k2)-4.*k2*k1*x2-pow2(k1))*G[6]*G[10]);
//
//   //5. F2/G2[k,-k-p]
//
//   F[16] =1./a*(-(a3*G[1]*G[4]+a4*G[5]*G[0])/2.-G[17]);
//   F[17] =1./a*(-(2.-hade1)*G[17]-hade2*G[16]-hade3*(1.+3.*p1)*G[18]-b2*G[1]*G[5]);
//   F[18] =1./a*(-(1+p1)*(a3*G[3]*G[6]+a4*G[7]*G[2])/2.-(1+p2)*G[19]+3*(p2-p1)*G[18]);
//   F[19] =1./a*(-(2.-3*p2-hade1)*G[19]-hade2*G[16]-(hade3*(1.+3.*p1) - c1*pow2(k2))*G[18]- b2*G[3]*G[7]
//                 +3.*(p1-p2)/(1.+1e-8+p2)*(a3*G[3]*G[6]+a4*G[7]*G[2])/2. - p1*a/(1.+1e-8+p2)*(a3*G[3]*F[6]+a4*G[7]*F[2])/2.
//                 -4.*(p1+1)*c1/(1.+1e-8+p2)*(-4.*pow2(k1)-4.*k2*k1*x2-pow2(k2))*G[6]*G[2]);
//
//   //6. F2/G2(p,k)
// 	F[20] =1./a*(-(a5*G[9]*G[0]+a6*G[1]*G[8])/2.-G[21]) ;
// 	F[21] =1./a*(-(2.-hade1)*G[21]-hade2*G[20]-hade3*(1.+3.*p1)*G[22]-b3*G[9]*G[1]);
//
//   F[22] =1./a*(-(1+p1)*(a5*G[11]*G[2]+a6*G[3]*G[10])/2.-(1+p2)*G[23]+3*(p2-p1)*G[22]) ;
// 	F[23] =1./a*(-(2.-3.*p2-hade1)*G[23]-hade2*G[20]-(hade3*(1.+3.*p1) - c1*pow2(karg1))*G[22]-b3*G[11]*G[3]
//                     +3.*(p1-p2)/(1.+1e-8+p2)*(a5*G[11]*G[2]+a6*G[3]*G[10])/2. - p1*a/(1.+1e-8+p2)*(a5*G[11]*F[2]+a6*G[3]*F[10])/2.
//                     -4.*(p1+1)*c1/(1.+1e-8+p2)*(k2*k2+k1*k1-2*k2*k1*x2)/2.*G[2]*G[10]);
//
//
//                     //2. F1/G1(k-p)
//                     // F[4] = -G[5]/a;
//                   	// F[5] =1./a*(-(2.-hade1)*G[5]-hade2*G[4]-hade3*G[6]);
//                     // F[6] =1./a*(-(1+p2)*G[7]+3*(p2-p1)*G[6]);
//                     // F[7] =1./a*(-(2.-3*p2-hade1)*G[7]-hade2*G[4]-(hade3*(1.+3.*p1) - c1*(k2*k2+k1*k1-2*k2*k1*x2))*G[6]);
//
//                       // //5. F2/G2(p,p)
//                     	// F[16] =1./a*(-a3*G[8]*G[9]-G[17]) ;
//                     	// F[17] =1./a*(-(2.-hade1)*G[17]-hade2*G[16]-hade3*G[18]-b4*G[9]*G[9]);
//                       // F[18] =1./a*(-(1+p1)*a3*G[10]*G[11]-(1+p2)*G[19]+3.*(p2-p1)*G[18]);
//                     	// F[19] =1./a*(-(2.-3.*p2-hade1)*G[19]-hade2*G[16]-(hade3*(1+3.*p1) - c1*(2.*pow2(k2)*(1.+x4)))*G[18] - b4*G[11]*G[11]
//                       //               +3.*(p1-p2)/(1.+1e-8+p2)*a3*G[11]*G[10] - p1*a/(1.+1e-8+p2)*a3*G[11]*F[10]
//                       //               -4.*(p1+1)*c1/(1.+1e-8+p2)*(pow2(k2)*(1.-x4))*G[6]*G[10]);
//
// 	return GSL_SUCCESS;
// }


// Loops to construct numerical kernels
// Kmax and Kmin indicates the maximum and minimum values of K
// k loop dictates the angular parameter x =k1.k/(k^2*r)
// j loop dictates the magnitude of k1 (k1=rk)
//
// Default initialization at only scale factor A
// k is your k vector dependence (See mu(a,k))
// YMAX = QMAXp/kmin and YMIN=QMINp/kmax where kmin and kmax are the k-limits and QMINp and QMAXp are k1 limits (integrated over)
// set k = 1 and YMIN = QMINp/kmax, YMAX = QMAXp/kmin for no k -dependence (initialize once)
// set kmin = kmax = k for k dependence (initialize for each k)
//
// EDIT : Add in new gravity parameters to function's input as required (par1 is the nDGP parameter omega_rc as default)
// EDIT : if more than 1 gravity parameter is needed, one must also edit initn function in SPT.h file.


// int IOW::initn_cde(double A, double YMIN, double YMAX, double k, double omega0, double par1, double par2, double par3)
// {
// //#pragma omp parallel for schedule(dynamic)
// 	for ( int i = 0; i < n1; ++i)
// 	{
// 		for (int j = 0; j< n2; ++j)
// 		{
// 				real x1,x2,x3,k1,k2,k3,karg1,karg2,karg3,karg4;
//
// 				/* k - magnitude sampling */
// 				  double R = (YMIN) * exp(j*log(YMAX/YMIN)/(n2*1.-1.)); //exponential sampling
// 				//  	double R = (j*1.)/(n2*1.-1)*(YMAX-YMIN) + YMIN; // linear sampling
// 				//	double R = pow2(j*1./(n2*1.-1))*YMAX+YMIN; //quadratic sampling
//
// 				//k,k1,-k1
//
// 				/*Parameter selection for kernels */
//
// 				k1 = k;
// 		    k2 = k*R;
// 				k3 = k2;
//
// 				x1 = XMIN;
//         x2 = x128[i]; // x128 goes from - 0.99995 to 0.99995
// 				x3 = -x2;
//
//         karg1 = sqrt(k2*k2+k1*k1-2*k2*k1*x2);
//         karg2 = sqrt(k3*k3+2*k3*k1*x3+k1*k1);
//         karg3 = sqrt(k2*k2+2*k2*k3*x1+k3*k3);
//         karg4 = sqrt(k1*k1+2*k1*k2*x2+k2*k2);
//
// 			// Initial scale factor for solving system of equations
// 			//Sets 1-3
// 				double a = 0.0001;
//
// 				// Einstein de Sitter initial condtions for 2nd order kernels
//
// 				//4. F2/G2(p,k-p) (P22)
// 				// double EdsF2=a*a*(10./14.*alphas(k2,karg1,(k1*x2-k2)/karg1)+ 2./7.*beta1(k2,karg1,(k1*x2-k2)/karg1));
// 				// double EdsG2=-a*a*(6./14.*alphas(k2,karg1,(k1*x2-k2)/karg1)+ 4./7.*beta1(k2,karg1,(k1*x2-k2)/karg1));
//         //
// 				// //5. F2/G2(-k,k-p)
// 				// double EdsF2C=a*a*(10./14.*alphas(k1,karg1,(k2*x2-k1)/karg1)+ 2./7.*beta1(k1,karg1,(k2*x2-k1)/karg1));
// 				// double EdsG2C=-a*a*(6./14.*alphas(k1,karg1,(k2*x2-k1)/karg1) + 4./7.*beta1(k1,karg1,(k2*x2-k1)/karg1));
//         //
// 				// //6. F2/G2(p,-p)
// 				// double EdsF2A= 0.;
// 				// double EdsG2A= 0.; //Setting manually to 0 because xmin is not identically -1
//         //
// 				// //7. F2/G2(p,k)
// 				// double EdsF2B= a*a*(10./14.*alphas(k2,k1,x2)+ 2./7.*beta1(k2,k1,x2));
// 				// double EdsG2B= -a*a*(6./14.*alphas(k2,k1,x2) + 4./7.*beta1(k2,k1,x2));
// 				// //8. F2/G2(-p,k)=F2/G2(p,-k)
// 				// double Edsf2d= a*a*(10./14.*alphas(k3,k1,x3)+ 2./7.*beta1(k3,k1,x3));
// 				// double Edsg2d= -a*a*(6./14.*alphas(k3,k1,x3)+ 4./7.*beta1(k3,k1,x3));
//         //
// 				// // Einstein de Sitter initial condtions for 3rd order kernels
// 				// //9. 3rd order  ~ F3/G3(k,p,-p)
// 				// double Edsf3 = pow3(a)/3.*(2./63.*beta1(k1,karg3,(k2*x2+k3*x3)/karg3)*(beta1(k2,k3,x1)+6./4.*alphas(k2,k3,x1))
// 				// 							+1./18.*alpha(k1,karg3,(k2*x2+k3*x3)/karg3)*(beta1(k2,k3,x1)+10./2.*alphas(k2,k3,x1))
// 				// 						    +1./9.*alpha(karg3,k1,(k2*x2+k3*x3)/karg3)*(beta1(k2,k3,x1)+6./4.*alphas(k2,k3,x1))
//         //
// 				// 						   +2./63.*beta1(k3,karg4,(k1*x3+k2*x1)/karg4)*(beta1(k1,k2,x2)+6./4.*alphas(k1,k2,x2))
// 				// 						   +1./18.*alpha(k3,karg4,(k1*x3+k2*x1)/karg4)*(beta1(k1,k2,x2)+10./2.*alphas(k1,k2,x2))
// 				// 						   +1./9.*alpha(karg4,k3,(k1*x3+k2*x1)/karg4)*(beta1(k1,k2,x2)+6./4.*alphas(k1,k2,x2))
//         //
// 				// 						   +2./63.*beta1(k2,karg2,(k1*x2+k3*x1)/karg2)*(beta1(k1,k3,x3)+6./4.*alphas(k1,k3,x3))
// 				// 						   +1./18.*alpha(k2,karg2,(k1*x2+k3*x1)/karg2)*(beta1(k1,k3,x3)+10./2.*alphas(k1,k3,x3))
// 				// 						   +1./9.*alpha(karg2,k2,(k1*x2+k3*x1)/karg2)*(beta1(k1,k3,x3)+6./4.*alphas(k1,k3,x3)));
//         //
// 				// double Edsg3 = -pow3(a)/3.*(2./21.*beta1(k1,karg3,(k2*x2+k3*x3)/karg3)*(beta1(k2,k3,x1)+6./4.*alphas(k2,k3,x1))
// 				// 					   +1./42.*alpha(k1,karg3,(k2*x2+k3*x3)/karg3)*(beta1(k2,k3,x1)+10./2.*alphas(k2,k3,x1))
// 				// 					   +1./21.*alpha(karg3,k1,(k2*x2+k3*x3)/karg3)*(beta1(k2,k3,x1)+6./4.*alphas(k2,k3,x1))
//         //
// 				// 					   +2./21.*beta1(k3,karg4,(k1*x3+k2*x1)/karg4)*(beta1(k1,k2,x2)+6./4.*alphas(k1,k2,x2))
// 				// 					   +1./42.*alpha(k3,karg4,(k1*x3+k2*x1)/karg4)*(beta1(k1,k2,x2)+10./2.*alphas(k1,k2,x2))
// 				// 					   +1./21.*alpha(karg4,k3,(k1*x3+k2*x1)/karg4)*(beta1(k1,k2,x2)+6./4.*alphas(k1,k2,x2))
//         //
// 				// 					   +2./21.*beta1(k2,karg2,(k1*x2+k3*x1)/karg2)*(beta1(k1,k3,x3)+6./4.*alphas(k1,k3,x3))
// 				// 					   +1./42.*alpha(k2,karg2,(k1*x2+k3*x1)/karg2)*(beta1(k1,k3,x3)+10./2.*alphas(k1,k3,x3))
// 				// 					   +1./21.*alpha(karg2,k2,(k1*x2+k3*x1)/karg2)*(beta1(k1,k3,x3)+6./4.*alphas(k1,k3,x3)));
//         //
//
// //			  double G[18] = { a,-a,a,-a,a,-a,  EdsF2 ,EdsG2,  EdsF2C, EdsG2C, EdsF2A, EdsG2A , EdsF2B,  EdsG2B, Edsf2d, Edsg2d, Edsf3 ,Edsg3}; // initial conditions
// // Non-Eds ICs
// 			  double G[36] = { a,-a,a,-a,a,-a,a,-a,a,-a,a,-a, 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
//
// 			/*Parameters passed to system of equations */
// // EDIT : If more than one gravity parameter is used, add them after p1
// 				struct param_type3 my_params1 = {k1,x1,k2,x2,k3,x3,omega0, par1, par2, par3,karg1,karg2,karg3,karg4};
//
// 				gsl_odeiv2_system sys = {funcn4, jac, 36, &my_params1};
//
// 			//  this gives the system, the method, the initial interval step, the absolute error and the relative error.
// 				gsl_odeiv2_driver * d =
// 				gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
// 										   1e-6, 1e-6, 0.0);
//
// 				int status1 = gsl_odeiv2_driver_apply (d, &a, A, G);
//
//
// 				/*Allocation of array values */
//
// 			//F1(k;a), G1(k;a)
// 			F1_nk = G[0] ;
// 			G1_nk = G[1] ;
//
//       //F1q(k;a), G1q(k;a)
//       F1q_nk = G[2];
//       G1q_nk = G[3];
//
// 			// F1(k-p;a), G1(k-p;a)
// 			F1kmp_nk[i*n2 + j] = G[4];
// 			G1kmp_nk[i*n2 + j] = G[5];
//
// 			//F1(p;a), G1(p;a)
// 			F1p_nk[i*n2 + j] = G[8];
// 			G1p_nk[i*n2 + j] = G[9];
//
// 			/*2nd order*/
//
// 			//F2/G2(p,k-p) (P22)
// 			F2_nk[i*n2 + j] =  G[12];
// 			G2_nk[i*n2 + j] =  G[13];
//
// 			//F2/G2(-k,k-p)
// 			F2C_nk[i*n2 + j] =  G[16];
// 			G2C_nk[i*n2 + j] =  G[17];
//
//       //F2/G2(p,-p)
//       F2D_nk[i*n2 + j] =  G[20];
//       G2D_nk[i*n2 + j] =  G[21];
//
// 			//F2/G2(p,k)
// 			F2A_nk[i*n2 + j] =  G[24];
// 			G2A_nk[i*n2 + j] =  G[25];
//
// 			//F2/G2(p,-k)
// 			F2B_nk[i*n2 + j] =  G[28];
// 			G2B_nk[i*n2 + j] =  G[29];
//
// 			/* 3rd order */
// 			//F3/G3[k,p,-p]
// 			F3_nk[i*n2 + j]= G[32];
// 			G3_nk[i*n2 + j]= G[33];
//
// 			gsl_odeiv2_driver_free(d);
// 		}
// 	}
//
//   //  Solution for wCDM linear growth @ a=1 for our normalisation of the linear power spectrum
//     	double a = 0.0001;
//     	double G[2] = {a,-a}; // initial conditions
//     	struct param_type2 params_my = {omega0,par1,par2,0.,1};
//
//     // Solutions of evolution factors @ a=A
//     	gsl_odeiv2_system sys3 = {funcn3, jac, 2, &params_my};
//     	gsl_odeiv2_driver * d1 =
//     	gsl_odeiv2_driver_alloc_y_new (&sys3, gsl_odeiv2_step_rk8pd,
//     								   1e-6, 1e-6, 0.0);
//
//     	int status2 = gsl_odeiv2_driver_apply (d1, &a, 1. , G);
//
//     	dnorm_spt = G[0];
//       gnorm_spt = G[1];
//
//
//     	gsl_odeiv2_driver_free(d1);
//
// 	return 0;
// }











#if 0
/* Based on the Netlib routine (D)GAMMA by W. J. Cody and L. Stoltz.  Original
 * source and documentation available at
 *   http://netlib.org/specfun/gamma */
double Gamma(double x) {
    /* Constants */
    const double sqrtpi = 0.9189385332046727417803297;
    const double pi = 3.1415926535897932384626434;

    /* Numerator and denominator coefficients for rational minimax approximation over (1,2). */
    const double P[] = { -1.71618513886549492533811e0,
        2.47656508055759199108314e1, -3.79804256470945635097577e2,
        6.29331155312818442661052e2, 8.66966202790413211295064e2,
        -3.14512729688483675254357e4, -3.61444134186911729807069e4,
        6.64561438202405440627855e4 };
    const double Q[] = { -3.08402300119738975254353e1,
        3.15350626979604161529144e2, -1.01515636749021914166146e3,
        -3.10777167157231109440444e3, 2.25381184209801510330112e4,
        4.75584627752788110767815e3, -1.34659959864969306392456e5,
        -1.15132259675553483497211e5 };

    /* Coefficients for minimax approximation over (12,infty). */
    const double C[] = { -1.910444077728e-3,8.4171387781295e-4,
        -5.952379913043012e-4, 7.93650793500350248e-4,
        -2.777777777777681622553e-3, 8.333333333333333331554247e-2,
        5.7083835261e-3 };

    /* Machine dependent parameters (reasonable values here) */
    const double xbig = 171.624;
    const double xminin = 2.23e-308;
    const double eps = 2.22e-16;
    const double xinf = 1.79e308;

    double fact, res, y;
    int i;

    fact = 1;
    y = x;

    if(x <= 0) {
        y = -x;
        if(y - (int)y == 0) {
            /* x is a negative integer */
            return xinf;
        }
        else {
            /* Use the reflection formula Gamma(1-x) Gamma(x) = pi/sin(pi x) */
            fact = pi/sin(pi*x);
            y = 1 - x;
        }
    }

    /* y is now positive, and we seek Gamma(y) */

    if(y < eps) {
        /* 0 < y < eps: use limiting formula Gamma(y) -> 1/y */
        if(y >= xminin)
            res = 1/y;
        else
            res = xinf;
    }
    else if(y < 12) {
        /* eps < y < 12: use rational function approximation and recursion formula */
        int n = ((int)y) - 1;
        double z = y - (n + 1);
        double xnum, xden;

        /* Evaluate minimax approximation for Gamma(1+z) with 0 < z < 1 */
        xnum = 0;
        xden = 1;
        for(i = 0; i < 8; i++) {
            xnum = (xnum + P[i])*z;
            xden = xden*z + Q[i];
        }
        res = 1 + xnum/xden;

        /* Adjust result for y < 1 or y > 2 */
        if(n == -1)
            res /= y;
        else
            for(i = 0; i < n; i++)
                res *= (z + 1 + i);
    }
    else {
        /* y >= 12: use asymptotic expansion */
        if(y > xbig)
            res = xinf;
        else {
            double ysq = y*y;
            double sum = C[6];
            for(i = 0; i < 6; i++)
                sum = sum/ysq + C[i];
            sum = sum/y - y + sqrtpi + (y-0.5)*log(y);
            res = exp(sum);
        }
    }

    if(fact != 1)
        res = fact/res;
    return res;
}
#endif

#if 0

/******************************************************************************
 * Implementation of complete and incomplete gamma functions, following
 *   N. M. Temme, "A set of algorithms for the incomplete gamma function", 1994.
 * Adapted from the Pascal source code presented in that paper.  Results are
 * accurate to 9 significant digits.
 ******************************************************************************/

struct MachineConstants {
    double machtol, dwarf, giant;
    double sqrtgiant, sqrtdwarf, lndwarf, lnmachtol, sqrtminlnmachtol, oneoversqrt2mt,
           explow, sqrtminexplow, exphigh, sqrttwopi, lnsqrttwopi, sqrtpi, oneoversqrtpi;

    MachineConstants() {
        machtol = DBL_EPSILON;
        dwarf = DBL_MIN;
        giant = DBL_MAX;
        sqrtgiant = sqrt(giant);
        sqrtdwarf = sqrt(dwarf);
        lndwarf = log(dwarf);
        lnmachtol = log(machtol);
        sqrtminlnmachtol = sqrt(-lnmachtol);
        oneoversqrt2mt = 1/sqrt(2*machtol);
        explow = lndwarf;
        sqrtminexplow = sqrt(-explow);
        exphigh = log(giant);
        sqrttwopi = sqrt(2*M_PI);
        lnsqrttwopi = log(sqrttwopi);
        sqrtpi = sqrt(M_PI);
        oneoversqrtpi = 1/sqrtpi;
    }
};

static MachineConstants constants;

/* Compute the rational function
 *   a_m x^m + ... + a_1 x + a_0
 *   ---------------------------
 *   b_n x^n + ... + b_1 x + b_0 */
static double ratfun(double x, int m, double a[], int n, double b[]) {
    int k;
    double num = a[m], den = b[n];
    for(k = m-1; k >= 0; k--)
        num = num*x + a[k];
    for(k = n-1; k >= 0; k--)
        den = den*x + b[k];
    return num/den;
}

/* Compute e^x - 1 */
static double exmin1(double x) {
    static double ak[4] = { 9.999999998390e-1, 6.652950247674e-2, 2.331217139081e-2, 1.107965764952e-3 };
    static double bk[4] = { 1.000000000000e+0,-4.334704979491e-1, 7.338073943202e-2,-5.003986850699e-3 };

    if(x < constants.lnmachtol)
        return -1.;
    else if(x > constants.exphigh)
        return constants.giant;
    else if(x < -0.69 || x > 0.41)
        return exp(x) - 1.;
    else if(fabs(x) < constants.machtol)
        return x;
    else
        return x * ratfun(x, 3, ak, 3, bk);
}

/* Compute ln(1+x) - x */
static double auxln(double x) {
    static double ak[5] = {-4.999999994526e-1,-5.717084236157e-1,-1.423751838241e-1,-8.310525299547e-4, 3.899341537646e-5 };
    static double bk[4] = { 1.000000000000e+0, 1.810083408290e+0, 9.914744762863e-1, 1.575899184525e-1 };

    if(x <= -1.)
        return -constants.giant;
    else if(x < -0.70 || x > 1.36)
        return log(1+x) - x;
    else if(fabs(x) < constants.machtol)
        return -0.5*x*x;
    else {
        if(x > 0)
            return x*x*ratfun(x, 4, ak, 3, bk);
        else {
            double z = -x/(1+x);
            if(z > 1.36)
                return -(log(1+z) - z) + x*z;
            else
                return -z*z*ratfun(z, 4, ak, 3, bk) + x*z;
        }
    }
}

/* Compute Gamma^*(x) */
static double gammastar(double x) {
    static double ak12[3] = { 1.000000000949e+0, 9.781658613041e-1, 7.806359425652e-2 };
    static double bk12[2] = { 1.000000000000e+0, 8.948328926305e-1 };
    static double ak[4] = { 5.115471897484e-2, 4.990196893575e-1, 9.404953102900e-1, 9.999999625957e-1 };
    static double bk[4] = { 1.544892866413e-2, 4.241288251916e-1, 8.571609363101e-1, 1.000000000000e+0 };

    if(x > 1e10) {
        if(x > 1./(12.*constants.machtol))
            return 1.;
        else
            return 1. + 1./(12.*x);
    }
    else if(x >= 12.)
        return ratfun(1/x, 2, ak12, 1, bk12);
    else if(x >= 1.)
        return ratfun(x, 3, ak, 3, bk);
    else if(x > constants.dwarf)
        return gammastar(x+1) * sqrt(1+1/x) * exp(-1 + x*log(1+1/x));
    else
        return 1./(constants.sqrttwopi*constants.sqrtdwarf);
}

double Gamma(double x) {
    static double ak[5] = { 1.000000000000e+0,-3.965937302325e-1, 2.546766167439e-1,-4.880928874015e-2, 9.308302710346e-3 };
    static double bk[5] = {-1.345271397926e-1, 1.510518912977e+0,-6.508685450017e-1, 9.766752854610e-2,-5.024949667262e-3 };
    double a, g, s, dw;
    int j, k, m;

    if(x <= constants.dwarf)
        return 1./constants.dwarf;
    else {
        k = (int) round(x);
        m = (int) trunc(x);
        dw = (k == 0) ? constants.dwarf : (1+x)*constants.machtol;
        if(fabs(k - x) < dw && x <= 15.) {
            /* x = k  ==>  Gamma(x) = (k-1)! */
            g = 1.;
            for(j = 1; j <= k-1; j++)
                g *= j;
            return g;
        }
        else if(fabs((x-m)-0.5) < (1+x)*constants.machtol && x <= 15.) {
            /* x = m + 1/2  ==> use recursion and Gamma(0.5) = sqrt(pi) */
            g = constants.sqrtpi;
            for(j = 1; j <= m; j++)
                g *= (j-0.5);
            return g;
        }
        else if(x < 1.)
            return ratfun(x+2, 4, ak, 4, bk) / (x*(x+1));
        else if(x < 2.)
            return ratfun(x+1, 4, ak, 4, bk) / x;
        else if(x < 3.)
            return ratfun(x, 4, ak, 4, bk);
        else if(x < 10.) {
            g = 1.;
            while(x >= 3.) {
                x -= 1;
                g *= x;
            }
            return g * ratfun(x, 4, ak, 4, bk);
        }
        else if(x < constants.exphigh) {
            a = 1/(x*x);
            g = (1. + a*(-3.33333333333e-2 + a*9.52380952381e-3)) / (12.*x);
            a = -x + (x-0.5)*log(x) + g + constants.lnsqrttwopi;
            return (a < constants.exphigh) ?  exp(a) : constants.giant;
        }
        else {
        }
    }
}

/* Compute the function g(x) in the representation
 *   1/Gamma(1+x) = 1 + x*(x-1)*g(x) */
static double auxgam(double x) {
    static double ak[4] = {-5.772156647338e-1,-1.087824060619e-1, 4.369287357367e-2,-6.127046810372e-3 };
    static double bk[5] = { 1.000000000000e+0, 3.247396119172e-1, 1.776068284106e-1, 2.322361333467e-2, 8.148654046054e-3 };

    if(x <= -1.)
        return -0.5;
    else if(x < 0.)
        return -(1 + (x+1)*(x+1)*ratfun(x+1, 3, ak, 4, bk)) / (1-x);
    else if(x <= 1.)
        return ratfun(x, 3, ak, 4, bk);
    else if(x <= 2.)
        return ((x-2)*ratfun(x-1, 3, ak, 4, bk) - 1) / (x*x);
    else
        return (1/Gamma(x+1) - 1.) / (x*(x-1));
}

/* Compute ln Gamma(x) */
static double lngamma(double x) {
    static double ak4[5] =  {-2.12159572323e5, 2.30661510616e5, 2.74647644705e4,-4.02621119975e4,-2.29660729780e3 };
    static double bk4[5] =  {-1.16328495004e5,-1.46025937511e5,-2.42357409629e4,-5.70691009324e2, 1.00000000000e0 };
    static double ak15[5] = {-7.83359299449e1,-1.42046296688e2, 1.37519416416e2, 7.86994924154e1, 4.16438922228 };
    static double bk15[5] = { 4.70668766060e1, 3.13399215894e2, 2.63505074721e2, 4.33400022514e1, 1.00000000000 };
    static double ak0[5] =  {-2.66685511495,  -2.44387534237e1,-2.19698958928e1, 1.11667541262e1, 3.13060547623 };
    static double bk0[5] =  { 6.07771387771e-1,1.19400905721e1, 3.14690115749e1, 1.52346874070e1, 1.00000000000 };

    double a, g, y;
    if(x > 12.) {
        g = 1./(12.*x);
        a = -x + (x-0.5)*log(x) + constants.lnsqrttwopi;
        if(a + g == a)
            return a;
        y = 1./(x*x);
        return a + g*(1. + y*(-3.33333333333e-2 + y*9.52380952381e-3));
    }
    else if(x >= 4.)
        return ratfun(x, 4, ak4, 4, bk4);
    else if(x > 1.5)
        return (x-2) * ratfun(x, 4, ak15, 4, bk15);
    else if(x >= 0.5)
        return (x-1) * ratfun(x, 4, ak0, 4, bk0);
    else if(x > constants.machtol)
        return -log(x) + x*ratfun(x+1, 4, ak0, 4, bk0);
    else if(x > constants.dwarf)
        return -log(x);
    else
        return -constants.lndwarf;
}

/* Compute the normalized (upper) incomplete gamma function
 *   Q(a,x) = \Gamma(a,x) / \Gamma(a) */
static double incomgam(double a, double x, double eps = 1e-10) {
    double lnx, mu, auxlnmu, dp, pqasymp;
    if(a == 0. && x == 0.)
        return 0.5;
    else if(x == 0.)
        return 1.;
    else if(a == 0.)
        return 0.;
    else {
        lnx = (x <= constants.dwarf) ? constants.lndwarf : log(x);

        /* function dax */
        mu = (x-a)/a;
        auxlnmu = auxln(mu);
        dp = a*auxlnmu - 0.5*log(2*M_PI*a);
        if(dp < constants.explow)
            dp = 0.;
        else
            dp = exp(dp) / gammastar(a);

        if(dp >= constants.dwarf) {
            if(a > 25. && fabs(mu) < 0.2)
                return pqasymp();
            else if(a > alfa(x))
                return ptaylor();
            else if(x < 1.)
                return qtaylor();
            else
                return qfraction();
        }
        else
            return (a > x) ? 1. : 0.;
    }
}

double Gamma(double x) {
    return tgamma(x);
}

double LowerGamma(double a, double x) {
    if(x < 0.1)
        return pow(x,a)*(1/a - x/(a+1) + 0.5*x*x/(a+2));
    else
        return gsl_sf_gamma(a) - gsl_sf_gamma_inc(a, x);
}

double UpperGamma(double a, double x) {
    return gsl_sf_gamma_inc(a, x);
}
#endif // 0

double Gamma(double x) {
    /* Use POSIX call for simplicity */
    return tgamma(x);
}

double LogGamma(double x) {
    /* Use POSIX call for simplicity */
    return lgamma(x);
}


#if 0
/* Modified from http://people.sc.fsu.edu/~jburkardt/c_src/asa239/asa239.c .
 * Original comments follow:
 *    Purpose:
 *      ALNORM computes the cumulative density of the standard normal distribution.
 *    Licensing:
 *      This code is distributed under the GNU LGPL license.
 *    Modified:
 *      13 November 2010
 *    Author:
 *      Original FORTRAN77 version by David Hill.
 *      C version by John Burkardt.
 *    Reference:
 *      David Hill,
 *      Algorithm AS 66:
 *      The Normal Integral,
 *      Applied Statistics,
 *      Volume 22, Number 3, 1973, pages 424-427.
 *    Parameters:
 *     -Input, double X, is one endpoint of the semi-infinite interval
 *      over which the integration takes place.
 *     -Input, int UPPER, determines whether the upper or lower
 *      interval is to be integrated:
 *      1  => integrate from X to + Infinity;
 *      0 => integrate from - Infinity to X.
 *     -Output, double ALNORM, the integral of the standard normal
 *      distribution over the desired interval.  */
static double asa239_alnorm(double x, int upper) {
  double a1 = 5.75885480458;
  double a2 = 2.62433121679;
  double a3 = 5.92885724438;
  double b1 = -29.8213557807;
  double b2 = 48.6959930692;
  double c1 = -0.000000038052;
  double c2 = 0.000398064794;
  double c3 = -0.151679116635;
  double c4 = 4.8385912808;
  double c5 = 0.742380924027;
  double c6 = 3.99019417011;
  double con = 1.28;
  double d1 = 1.00000615302;
  double d2 = 1.98615381364;
  double d3 = 5.29330324926;
  double d4 = -15.1508972451;
  double d5 = 30.789933034;
  double ltone = 7.0;
  double p = 0.398942280444;
  double q = 0.39990348504;
  double r = 0.398942280385;
  int up;
  double utzero = 18.66;
  double value;
  double y;
  double z;

  up = upper;
  z = x;

  if ( z < 0.0 )
  {
    up = !up;
    z = - z;
  }

  if ( ltone < z && ( ( !up ) || utzero < z ) )
  {
    if ( up )
    {
      value = 0.0;
    }
    else
    {
      value = 1.0;
    }
    return value;
  }

  y = 0.5 * z * z;

  if ( z <= con )
  {
    value = 0.5 - z * ( p - q * y
      / ( y + a1 + b1
      / ( y + a2 + b2
      / ( y + a3 ))));
  }
  else
  {
    value = r * exp ( - y )
      / ( z + c1 + d1
      / ( z + c2 + d2
      / ( z + c3 + d3
      / ( z + c4 + d4
      / ( z + c5 + d5
      / ( z + c6 ))))));
  }

  if ( !up )
  {
    value = 1.0 - value;
  }

  return value;
}
#endif



/* Modified from http://people.sc.fsu.edu/~jburkardt/c_src/asa239/asa239.c .
 * Original comments follow:
 *    Purpose:
 *      GAMMAD computes the Incomplete Gamma Integral
 *    Licensing:
 *      This code is distributed under the GNU LGPL license.
 *    Modified:
 *      13 November 2010
 *    Author:
 *      Original FORTRAN77 version by B Shea.
 *      C version by John Burkardt.
 *    Reference:
 *      B Shea,
 *      Algorithm AS 239:
 *      Chi-squared and Incomplete Gamma Integral,
 *      Applied Statistics,
 *      Volume 37, Number 3, 1988, pages 466-473.
 *    Parameters:
 *     -Input, double X, P, the parameters of the incomplete
 *      gamma ratio.  0 <= X, and 0 < P.
 *     -Output, int IFAULT, error flag.
 *      0, no error.
 *      1, X < 0 or P <= 0.
 *     -Output, double GAMMAD, the value of the incomplete
 *      Gamma integral. */
/* Note: alnorm(x) = 0.5*(1 + erf(x/sqrt(2)) */
static double asa239_gammad(double x, double p, int *ifault) {
    const double elimit = -88.0;
    const double oflo = 1.0E+37;
    const double plimit = 1000.0;
    const double tol = 1.0E-14;
    const double xbig = 1.0E+08;
    double a;
    double an;
    double arg;
    double b;
    double c;
    double pn1, pn2, pn3, pn4, pn5, pn6;
    double rn;
    double value = 0.0;

    value = 0.0;
    *ifault = 0;

    /* Check the input. */
    if(x < 0.0 || p <= 0.0) {
        *ifault = 1;
        return 0;
    }

    if(x == 0.0)
        return 0;

    /* If P is large, use a normal approximation. */
    if(p > plimit) {
        pn1 = 3*sqrt(p) * ( pow(x/p, 1/3.) + 1/(9*p) - 1 );
        value = 0.5 + 0.5*erf(pn1/M_SQRT2);
        return value;
    }

    /* If X is large set value = 1. */
    if(xbig < x)
        return 1;

    if(x <= 1.0 || x < p) {
        /* Use Pearson's series expansion. */

        arg = p*log(x) - x - lgamma(p + 1);
        c = 1.0;
        value = 1.0;
        a = p;

        while(c > tol) {
            a = a + 1;
            c = c * x / a;
            value += c;
        }

        arg += log(value);

        if(arg >= elimit)
            value = exp(arg);
        else
            value = 0.0;
    }
    else {
        /* Use a continued fraction expansion. */

        arg = p*log(x) - x - lgamma(p);
        a = 1 - p;
        b = a + x + 1;
        c = 0.0;
        pn1 = 1.0;
        pn2 = x;
        pn3 = x + 1;
        pn4 = x * b;
        value = pn3 / pn4;

        for ( ; ; ) {
            a += 1;
            b += 2;
            c += 1;
            an = a * c;
            pn5 = b * pn3 - an * pn1;
            pn6 = b * pn4 - an * pn2;

            if(pn6 != 0.0) {
                rn = pn5 / pn6;

                if(fabs(value - rn) <= fmin(tol, tol * rn))
                    break;
                value = rn;
            }

            pn1 = pn3;
            pn2 = pn4;
            pn3 = pn5;
            pn4 = pn6;

            /* Re-scale terms in continued fraction if terms are large. */
            if(fabs(pn5) >= oflo) {
                pn1 = pn1 / oflo;
                pn2 = pn2 / oflo;
                pn3 = pn3 / oflo;
                pn4 = pn4 / oflo;
            }
        }

        arg += log(value);

        if(arg >= elimit)
            value = 1 - exp(arg);
        else
            value = 1;
    }

    return value;
}

double LowerGamma(double a, double x) {
    int ifault;
    return Gamma(a) * asa239_gammad(x, a, &ifault);
}

double UpperGamma(double a, double x) {
    int ifault;
    return Gamma(a) * (1 - asa239_gammad(x, a, &ifault));
}
