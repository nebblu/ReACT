#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <cfloat>
#include <cmath>

#include "Common.h"
#include "Quadrature.h"
#include "SpecialFunctions.h"
#include "BSPTN.h"

//#include<omp.h>

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_bessel.h>
#include <math.h>       /* pow */
#include <iostream>
#include <stdlib.h>
#include <functional>


using std::cref;
using std::bind;

// 1st ORDER
double F1b_k1;
double G1b_k1;
double F1b_k2;
double G1b_k2;
double F1b_k3;
double G1b_k3;

double F1b_1;
double G1b_1;

double F1b_2;
double G1b_2;

double F1b_3;
double G1b_3;


// 2nd order
double F2b_k12;
double G2b_k12;
double F2b_k13;
double G2b_k13;
double F2b_k23;
double G2b_k23;

double F2b_p2mp;
double G2b_p2mp;

double F2b_p3mp;
double G2b_p3mp;

double F2b_12a;
double G2b_12a;

double F2b_13a;
double G2b_13a;

double F2b_23a;
double G2b_23a;

// 3rd order

double F3b_12mp;
double G3b_12mp;

double F3b_13mp;
double G3b_13mp;

double F3b_21mp;
double G3b_21mp;

double F3b_23mp;
double G3b_23mp;

double F3b_31mp;
double G3b_31mp;

double F3b_32mp;
double G3b_32mp;

double F3b_1pp;
double G3b_1pp;

double F3b_2pp;
double G3b_2pp;

double F3b_3pp;
double G3b_3pp;

// 4th order
double F4b_12pp;
double G4b_12pp;
double F4b_13pp;
double G4b_13pp;
double F4b_23pp;
double G4b_23pp;

/* Euler and Continuity equations for tree level Bispectrum numerical kernels */

//alphai(k1,k2)
inline double alphai(double k1, double k2, double u1){
	return 1.+k2*u1/k1;
}

//beta(k1,k2)
inline double betai(double k1, double k2, double u1){
	return u1*(k1*k1+k2*k2+2.*k1*k2*u1)/(2.*k1*k2);
}

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

  real arg1;
  real arg2;
  real arg3;

  real omega00;
	real par1;
	real par2;
	real par3;
};

int jacb (double a, const double G[], double *dfdy, double dfdt[], void *params)
{
	return GSL_SUCCESS;
}


int funcb1dgp(double a, const double G[], double F[], void *params)
{
	param_type3 p = *(param_type3 *)(params);
	real k1 = p.kk1;
	real x1= p.xx1;
	real k2 = p.kk2;
	real x2= p.xx2;
	real k3 = p.kk3;
	real x3= p.xx3;
//  real k12 = p.arg1;
//  real k13 = p.arg2;
//  real k23 = p.arg3;

	real omega0 = p.omega00;
	real p1 = p.par1;
	real p2 = p.par2;
	real p3 = p.par3;

  double a1,a2,a3,a4,a5,a6;
  double b1,b2,b3;

    a1 = alphai(k2,k3,x1);//-(k1*x2+k2)/k3);
    a2 = alphai(k3,k2,x1);//-(k1*x2+k2)/k3);
    b1 = betai(k2,k3,x1);//-(k1*x2+k2)/k3);

    a3 = alphai(k1,k3,x3);//-(k2*x2+k1)/k3);
    a4 = alphai(k3,k1,x3);//-(k2*x2+k1)/k3);
    b2 = betai(k1,k3,x3);//-(k2*x2+k1)/k3);

    a5 = alphai(k2,k1,x2);
    a6 = alphai(k1,k2,x2);
    b3 = betai(k1,k2,x2);


    /// Select gravity ////
      double hub = HA2(a,omega0);
      double hub1= HA(a,omega0);
      double hubsqr = pow2(hub1);
      double acub = pow3(a);
      // DGP

      double betadgp = 1./(3.*(mu(a,1.,omega0,p1,p2,p3,3)-1.));
      double betadgpcub = pow3(betadgp);
      double mua = 1.+1./(3.*betadgp);
      double gam2= -1./(hubsqr*24.*betadgpcub*p1)*pow2(omega0/acub);
      double gamk[3];
      gamk[0] = (1.-pow2(x1));
      gamk[1] = (1.-pow2(x2));
      gamk[2] = (1.-pow2(x3));

	/* 1st order */
	//1. F1/G1(k1)
	F[0] = -G[1]/a;
	F[1] =1./a*(-(2.-hub)*G[1]-hub*G[0]*mua);

	//2. F1/G1(k3)
	F[2] = -G[3]/a;
	F[3] =1./a*(-(2.-hub)*G[3]-hub*G[2]*mua);

	//3. F1/G1(k2)
	F[4] = -G[5]/a;
	F[5] =1./a*(-(2.-hub)*G[5]-hub*G[4]*mua);

	/* 2nd order */

	//4. F2/G2(k2,k3) (P22)
	F[6] =1./a*(-(a1*G[5]*G[2]+a2*G[3]*G[4])/2.-G[7]);
	F[7] =1./a*(-(2.-hub)*G[7]-hub*G[6]*mua - gam2*gamk[0]*G[4]*G[2] - b1*G[5]*G[3]);


	//5. F2/G2(k1,k3)
	F[8] =1./a*(-(a3*G[1]*G[2]+a4*G[3]*G[0])/2.-G[9]);
	F[9] =1./a*(-(2.-hub)*G[9]-hub*G[8]*mua - gam2*gamk[2]*G[2]*G[0] - b2*G[3]*G[1]);

	//7. F2/G2(k1,k2)
	F[10] =1./a*(-(a5*G[5]*G[0]+a6*G[1]*G[4])/2.-G[11]) ;
	F[11] =1./a*(-(2.-hub)*G[11]-hub*G[10]*mua - gam2*gamk[1]*G[4]*G[0]-b3*G[5]*G[1]);


	return GSL_SUCCESS;
}


// equations for fofr
int funcb1fofr(double a, const double G[], double F[], void *params)
{
	param_type3 p = *(param_type3 *)(params);
	real k1 = p.kk1;
	real x1= p.xx1;
	real k2 = p.kk2;
	real x2= p.xx2;
	real k3 = p.kk3;
	real x3= p.xx3;
  real k12 = p.arg1;
  real k13 = p.arg2;
  real k23 = p.arg3;

	real omega0 = p.omega00;
	real p1 = p.par1;
//	real p2 = p.par2; // unused
//	real p3 = p.par3; // unused

  double a1,a2,a3,a4,a5,a6;
  double b1,b2,b3;

    a1 = alphai(k2,k3,x1);//-(k1*x2+k2)/k3);
    a2 = alphai(k3,k2,x1);//-(k1*x2+k2)/k3);
    b1 = betai(k2,k3,x1);//-(k1*x2+k2)/k3);

    a3 = alphai(k1,k3,x3);//-(k2*x2+k1)/k3);
    a4 = alphai(k3,k1,x3);//-(k2*x2+k1)/k3);
    b2 = betai(k1,k3,x3);//-(k2*x2+k1)/k3);

    a5 = alphai(k2,k1,x2);
    a6 = alphai(k1,k2,x2);
    b3 = betai(k1,k2,x2);

      double hub = HA2(a,omega0);
      double hub1= HA(a,omega0);
      double hubsqr = pow2(hub1);
      double acub = pow3(a);
			double ap6 = pow2(acub);
			double ap12 = pow2(ap6);
 			double h0 = myh0sqr;
			double fofrp = p1/h0;
			double fofrp2 = pow2(fofrp);
			double om0a3 = pow2(omega0/acub);
			double term0 = pow2(3.*omega0-4.);
			double term0sqr = pow2(term0);
			double term1a = omega0-4.*acub*(omega0-1.);
			double term1b = pow2(term1a);
			double term1bsqr= pow2(term1b);
			double term1 = term1b*term1a/(2.*ap6*acub*fofrp*term0);
			double koa1 = pow2(k1/a) ;
			double koa2 = pow2(k2/a) ;
			double koa3 = pow2(k3/a) ;
			double koa23 = pow2(k23/a) ;
			double koa13 = pow2(k13/a) ;
			double koa12 = pow2(k12/a) ;

      // f(R)
			// could optimise in terms of computations
			// 1st order
      double muak1 = 1. + koa1/(3.*(koa1 + term1));
			double muak2 = 1. + koa2/(3.*(koa2 + term1));
			double muak3 = 1. + koa3/(3.*(koa3 + term1));
			double muak23 = 1. + koa23/(3.*(koa23 + term1));
			double muak13 = 1. + koa13/(3.*(koa13 + term1));
			double muak12 = 1. + koa12/(3.*(koa12 + term1));

      double gam223 =  -(9.*koa23*om0a3*term1bsqr*term1a)/
						    				(48.*ap12*acub*fofrp2*hubsqr*term0sqr
			 			   					*(koa23+term1)
			 			   					*(koa2+term1)
			 	 		   					*(koa3+term1));


			double gam213 =  -(9.*koa13*om0a3*term1bsqr*term1a)/
						    				(48.*ap12*acub*fofrp2*hubsqr*term0sqr
			 			   					*(koa13+term1)
			 			   					*(koa1+term1)
			 	 		   					*(koa3+term1));


  		double gam212 =  -(9.*koa12*om0a3*term1bsqr*term1a)/
					    			   	(48.*ap12*acub*fofrp2*hubsqr*term0sqr
		 			   				  	*(koa12+term1)
		 			   				  	*(koa2+term1)
		 	 		   				  	*(koa1+term1));

	/* 1st order */
	//1. F1/G1(k1)
	F[0] = -G[1]/a;
	F[1] =1./a*(-(2.-hub)*G[1]-hub*G[0]*muak1);

	//2. F1/G1(k3)
	F[2] = -G[3]/a;
	F[3] =1./a*(-(2.-hub)*G[3]-hub*G[2]*muak3);

	//3. F1/G1(k2)
	F[4] = -G[5]/a;
	F[5] =1./a*(-(2.-hub)*G[5]-hub*G[4]*muak2);

	/* 2nd order */

	//4. F2/G2(k2,k3) (P22)
	F[6] =1./a*(-(a1*G[5]*G[2]+a2*G[3]*G[4])/2.-G[7]);
	F[7] =1./a*(-(2.-hub)*G[7]-hub*G[6]*muak23 - gam223*G[4]*G[2] - b1*G[5]*G[3]);


	//5. F2/G2(k1,k3)
	F[8] =1./a*(-(a3*G[1]*G[2]+a4*G[3]*G[0])/2.-G[9]);
	F[9] =1./a*(-(2.-hub)*G[9]-hub*G[8]*muak13 - gam213*G[2]*G[0] - b2*G[3]*G[1]);

	//7. F2/G2(k1,k2)
	F[10] =1./a*(-(a5*G[5]*G[0]+a6*G[1]*G[4])/2.-G[11]) ;
	F[11] =1./a*(-(2.-hub)*G[11]-hub*G[10]*muak12 - gam212*G[4]*G[0]-b3*G[5]*G[1]);


	return GSL_SUCCESS;
}


void BSPTN::initnb0_dgp(double A, double k[], double x[], double kargs[], double omega0, double par1, double par2, double par3)
{
				double a = 0.0001;

        // Non-Eds ICs
			  double G[12] = { a,-a,a,-a,a,-a, 0.,0.,0.,0.,0.,0.};


			/*Parameters passed to system of equations */
				struct param_type3 my_params1 = {k[0],x[0],k[1],x[1],k[2],x[2], kargs[0], kargs[1], kargs[2], omega0, par1, par2, par3};


				gsl_odeiv2_system sys = {funcb1dgp, jacb, 12, &my_params1};

			//  this gives the system, the method, the initial interval step, the absolute error and the relative error.
				gsl_odeiv2_driver * d =
				gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
										   1e-6, 1e-6, 1e-6);

				int status1 = gsl_odeiv2_driver_apply (d, &a, A, G);


				/*Allocation of array values */

// For tree level initialisation

  /*1st order */

			//F1(k;a), G1(k;a)
			F1b_k1 = G[0] ;
			G1b_k1 = G[1] ;

			// F1(k3;a), G1(k3;a)
			F1b_k3 = G[2];
			G1b_k3 = G[3];

			//F1(k2;a), G1(k2;a)
			F1b_k2 = G[4];
			G1b_k2 = G[5];


			/*2nd order*/

			//F2/G2(k3,k2) for tree,  F2/G2(k1-p,k2+p) for B222
			F2b_k23 =  G[6];
			G2b_k23 =  G[7];

      //F2/G2(k1,k3) for tree, F2/G2(-p,k2+p) for B222
      F2b_k13 =  G[8];
      G2b_k13 =  G[9];

      //F2/G2(k1,k2) for tree, F2/G2(p,k1-p)
      F2b_k12 =  G[10];
      G2b_k12 =  G[11];


			gsl_odeiv2_driver_free(d);
}


void BSPTN::initnb0_fofr(double A, double k[], double x[], double kargs[], double omega0, double par1, double par2, double par3)
{
				double a = 0.0001;

        // Non-Eds ICs
			  double G[12] = { a,-a,a,-a,a,-a, 0.,0.,0.,0.,0.,0.};


			/*Parameters passed to system of equations */
				struct param_type3 my_params1 = {k[0],x[0],k[1],x[1],k[2],x[2], kargs[0], kargs[1], kargs[2], omega0, par1, par2, par3};


				gsl_odeiv2_system sys = {funcb1fofr, jacb, 12, &my_params1};

			//  this gives the system, the method, the initial interval step, the absolute error and the relative error.
				gsl_odeiv2_driver * d =
				gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
										   1e-6, 1e-6, 1e-6);

				int status1 = gsl_odeiv2_driver_apply (d, &a, A, G);


				/*Allocation of array values */

// For tree level initialisation

  /*1st order */

			//F1(k;a), G1(k;a)
			F1b_k1 = G[0] ;
			G1b_k1 = G[1] ;

			// F1(k3;a), G1(k3;a)
			F1b_k3 = G[2];
			G1b_k3 = G[3];

			//F1(k2;a), G1(k2;a)
			F1b_k2 = G[4];
			G1b_k2 = G[5];

			/*2nd order*/

			//F2/G2(k3,k2) for tree
			F2b_k23 =  G[6];
			G2b_k23 =  G[7];

      //F2/G2(k1,k3) for tree
      F2b_k13 =  G[8];
      G2b_k13 =  G[9];

      //F2/G2(k1,k2) for tree, F2/G2(p,k1-p)
      F2b_k12 =  G[10];
      G2b_k12 =  G[11];


			gsl_odeiv2_driver_free(d);
}



/* Euler and Continuity equations for 1-loop Bispectrum numerical kernels */


struct param_type5 {
  real kk1;
	real kk2;
	real kk3;
  real kk4;
  real kk5;
  real kk6;
  real kk7;
  real kk8;

	real xx1;
	real xx2;
	real xx3;
	real xx4;
	real xx5;
	real xx6;
	real xx7;
	real xx8;

  real arg1;
  real arg2;
  real arg3;
  real arg4;
  real arg5;
  real arg6;
  real arg7;
  real arg8;
  real arg9;
  real arg10;
  real arg11;
  real arg12;
  real arg13;
  real arg14;
  real arg15;
  real arg16;
  real arg17;
  real arg18;
  real arg19;
  real arg20;
  real arg21;
  real arg22;
  real arg23;
  real arg24;
  real arg25;

  real omega00;
	real par1;
	real par2;
	real par3;


};



int funcdgp(double a, const double G[], double F[], void *params)
{
	param_type5 p = *(param_type5 *)(params);
  // variables
	real k1 = p.kk1;  // p
	real k2 = p.kk2; //  k1-p
	real k3 = p.kk3; // k2+p
  real k4 = p.kk4; // k2-p
	real k5 = p.kk5; // k3-p
  real k6 = p.kk6; // k1
  real k7 = p.kk7; // k2
  real k8 = p.kk8; // k3 (= -k1-k2)


	real x1= p.xx1;  // k1-p.k2+p
	real x2= p.xx2;  // p.k1-p
	real x3= p.xx3; // -p.k2+p
	real x4= p.xx4; // p.k2-p
	real x5= p.xx5; // p.k3-p
  real x6= p.xx6; // p.k1
  real x7= p.xx7; // p.k2
  real x8= p.xx8; // p.k3

//  real k13 = p.arg1; // k2 + 2p
  real k16 = p.arg2; // k1 + p
  real k46 = p.arg3; // -k3-p
  real x46 = p.arg4; // k2-p . k1
  real x416 = p.arg5; // k2-p . p+k1
  real x67 = p.arg6; // k1.k2 = x
  real x146 = p.arg7; // p. -k3-p
  real x56 = p.arg8; // k3-p. k1
  real x516 = p.arg9; // k3-p. k1+p
  real x68 = p.arg10; // k1.k3
  real x37 = p.arg11; // k2.k1+p
  real x57 = p.arg12; // k3-p.k2
  real x35 = p.arg13; // k2+p.k3-p
  real x78 = p.arg14; // k2.k3
  real x157 = p.arg15; // - p . k1+p
  real x27 = p.arg16; // k1-p.k2
  real x746 = p.arg17; // k2. -k3-p
  real x48 = p.arg18; // k2-p.k3
  real x418 = p.arg19; // k2-p.p+k3
  real x28 = p.arg20; // k1-p . k3
  real x218 = p.arg21; // k1-p.p+k3
  real x63 = p.arg22; // k1.k2+p
  real x816 = p.arg23; // k3.k1+p
  real x38 = p.arg24; // k2+p.k3
  real x646 = p.arg25; // k1.k3+p

// Cosmological and gravity parameters
 	real omega0 = p.omega00;
	real p1 = p.par1;
//	real p2 = p.par2; // unused
//	real p3 = p.par3; // unused

  double alph[76];
  double beta[38];


/// Select gravity ////

  double hub = HA2(a,omega0);
  // GR
  // double mua = 1.;
  // double gam2 = 1.;
  // double gam3 = 1.;
  // double gam4 = 1.;
  // DGP
  double hub1= HA(a,omega0);
  double hubsqr = pow2(hub1);
  double acub = pow3(a);
  double betadgp = 1.+hub1/sqrt(p1)*(1.+HA1(a,omega0)/(3.*hubsqr));
  double betadgpcub = pow3(betadgp);
  double om0a = pow2(omega0/acub);
  double mua = 1.+1./(3.*betadgp);
  double gam2 = -1./(hubsqr*24.*betadgpcub*p1)*om0a;
  double gam3 = 1./(hubsqr*144.*betadgpcub*pow2(betadgp)*pow2(p1))*om0a*(omega0/acub);
  double gam4 = -1./(hubsqr*1728.*pow2(betadgpcub)*betadgp*pow3(p1))*pow2(om0a);

double gamk[30];
   //gamma2, gamma3 and gamma4 scale dependencies  (sorry but need to put them by hand for optimisation)


gamk[0] = 1.-pow2(x1); //k1-p.k2+p
gamk[1] = 1.-pow2(x2); // p.k1-p
gamk[2] = 1.-pow2(x3); // p.k2+p
gamk[3] = 1.-pow2(x4); // p. k2-p
gamk[4] = 1.-pow2(x5); // p.k3-p
gamk[5] = 1.-pow2(x6); // p.k1
gamk[6] = 1.-pow2(x7); // p.k2
gamk[7] = 1.-pow2(x8); // p.k3
gamk[8] = 1.-pow2(x46); // k1.k2-p
gamk[9] = 1.-pow2(x416); // p+k1.k2-p
gamk[10] = 1.-pow2(x67); // k1.k2
gamk[11] = 1.-pow2(x146); // -p.k3+p
gamk[12] = 1.-pow2(x56); // k1.k3-p
gamk[13] = 1.-pow2(x516); // k1+p.k3-p
gamk[14] = 1.-pow2(x68); // k1.k3
gamk[15] = 1.-pow2(x37); // k2.k1+p
gamk[16] = 1.-pow2(x57); // k2.k3-p
gamk[17] = 1.-pow2(x35); // k2+p.k3-p
gamk[18] = 1.-pow2(x78); // k2.k3
gamk[19] = 1.-pow2(x157); // -p.k1+p
gamk[20] = 1.-pow2(x27); // k2.k1-p
gamk[21] = 1.-pow2(x746); // -k2.k3+p
gamk[22] = 1.-pow2(x48); // k3.k2-p
gamk[23] = 1.-pow2(x418); // k3+p.k2-p
gamk[24] = 1.-pow2(x28); // k3.k1-p
gamk[25] = 1.-pow2(x218); // k3+p.k1-p
gamk[26] = 1.-pow2(x63); // k1.k2+p
gamk[27] = 1.-pow2(x816); // k3.k1+p
gamk[28] = 1.-pow2(x38); //k3.k2+p
gamk[29] = 1.-pow2(x646); //k1.k3+p


	/* 1st order */
	//1. F1/G1(k1) or F1/G1(p)
	F[0] = -G[1]/a;
	F[1] =1./a*(-(2.-hub)*G[1]-hub*G[0]*mua);

	//2. F1/G1(k3) or F1/G1(k2+p)
	F[2] = -G[3]/a;
	F[3] =1./a*(-(2.-hub)*G[3]-hub*G[2]*mua);

	//3. F1/G1(k2) or F1/G1(k1-p)
	F[4] = -G[5]/a;
	F[5] =1./a*(-(2.-hub)*G[5]-hub*G[4]*mua);

	/* 2nd order for B222 */
  alph[0] = alphai(k2,k3,x1);
  alph[1] = alphai(k3,k2,x1);
  beta[0] = betai(k2,k3,x1);

	//4. F2/G2(k2,k3)  or F2/G2(k1-p,k2+p)
	F[6] =1./a*(-(alph[0]*G[5]*G[2]+alph[1]*G[3]*G[4])/2.-G[7]);
	F[7] =1./a*(-(2.-hub)*G[7]-hub*G[6]*mua - gam2*gamk[0]*G[4]*G[2] - beta[0]*G[5]*G[3]);

  alph[2] = alphai(k1,k3,x3);
  alph[3] = alphai(k3,k1,x3);
  beta[1] = betai(k1,k3,x3);

	//5. F2/G2(-k1,k3)  or F2/G2(-p,k2+p)
	F[8] =1./a*(-(alph[2]*G[1]*G[2]+alph[3]*G[3]*G[0])/2.-G[9]);
	F[9] =1./a*(-(2.-hub)*G[9]-hub*G[8]*mua - gam2*gamk[2]*G[2]*G[0] - beta[1]*G[3]*G[1]);

  alph[4] = alphai(k2,k1,x2);
  alph[5] = alphai(k1,k2,x2);
  beta[2] = betai(k1,k2,x2);

	//7. F2/G2(k2,k1) or F2/G2(k1-p,p)
	F[10] =1./a*(-(alph[4]*G[5]*G[0]+alph[5]*G[1]*G[4])/2.-G[11]) ;
	F[11] =1./a*(-(2.-hub)*G[11]-hub*G[10]*mua - gam2*gamk[1]*G[4]*G[0]-beta[2]*G[5]*G[1]);


  /* 1st order for B321 - I */

  // F1/G1(k4) or F1/G1(k2-p)

  F[12] = -G[13]/a;
  F[13] =1./a*(-(2.-hub)*G[13]-hub*G[12]*mua);

  // F1/G1(k5) or F1/G1(k3-p)

  F[14] = -G[15]/a;
  F[15] =1./a*(-(2.-hub)*G[15]-hub*G[14]*mua);

  // F1/G1(k6) or F1/G1(k1)

  F[16] = -G[17]/a;
  F[17] =1./a*(-(2.-hub)*G[17]-hub*G[16]*mua);

  // F1/G1(k7) or F1/G1(k2)

  F[18] = -G[19]/a;
  F[19] =1./a*(-(2.-hub)*G[19]-hub*G[18]*mua);

  // F1/G1(k8) or F1/G1(k3)

  F[20] = -G[21]/a;
  F[21] =1./a*(-(2.-hub)*G[21]-hub*G[20]*mua);

  // 2nd order used in 3rd order equations of B321-I

  alph[6] = alphai(k4,k1,x4);
  alph[7] = alphai(k1,k4,x4);
  beta[3] = betai(k1,k4,x4);

  // F2/G2(k4,k1)  =  (k2-p,p)
  F[22] =1./a*(-(alph[6]*G[13]*G[0]+alph[7]*G[1]*G[12])/2.-G[23]) ;
  F[23] =1./a*(-(2.-hub)*G[23]-hub*G[22]*mua - gam2*gamk[3]*G[12]*G[0]-beta[3]*G[13]*G[1]);

  alph[8] = alphai(k5,k1,x5);
  alph[9] = alphai(k1,k5,x5);
  beta[4] = betai(k1,k5,x5);

  // F2/G2(k5,k1) = (k3-p,p)
  F[24] =1./a*(-(alph[8]*G[15]*G[0]+alph[9]*G[1]*G[14])/2.-G[25]) ;
  F[25] =1./a*(-(2.-hub)*G[25]-hub*G[24]*mua - gam2*gamk[4]*G[14]*G[0]-beta[4]*G[15]*G[1]);

  alph[10] = alphai(k6,k1,x6);
  alph[11] = alphai(k1,k6,x6);
  beta[5] = betai(k1,k6,x6);

  // F2/G2(k6,k1) = (k1,p)
  F[26] =1./a*(-(alph[10]*G[17]*G[0]+alph[11]*G[1]*G[16])/2.-G[27]) ;
  F[27] =1./a*(-(2.-hub)*G[27]-hub*G[26]*mua - gam2*gamk[5]*G[16]*G[0]-beta[5]*G[17]*G[1]);

  alph[12] = alphai(k6,k4,x46);
  alph[13] = alphai(k4,k6,x46);
  beta[6] = betai(k4,k6,x46);

  // F2/G2(k6,k4) = (k1,k2-p)
  F[28] =1./a*(-(alph[12]*G[17]*G[12]+alph[13]*G[13]*G[16])/2.-G[29]) ;
  F[29] =1./a*(-(2.-hub)*G[29]-hub*G[28]*mua - gam2*gamk[8]*G[16]*G[12]-beta[6]*G[17]*G[13]);

  // F3/G3(k6,k1,k4) = (k1,p,k2-p)

 alph[14] = alphai(k4,k16,x416);
 alph[15] = alphai(k1,k46,x146);
 alph[16] = alphai(k6,k7,x67);

 alph[17] = alphai(k16,k4,x416);
 alph[18] = alphai(k46,k1,x146);
 alph[19] = alphai(k7,k6,x67);

 beta[7] = betai(k4,k16,x416);
 beta[8] = betai(k1,k46,x146);
 beta[9] = betai(k7,k6,x67);

 F[30] = - 1./(3.*a)*(alph[14]*G[26]*G[13]
           + alph[15]*G[28]*G[1]
           + alph[16]*G[22]*G[17]

           + alph[17]*G[27]*G[12]
           + alph[18]*G[29]*G[0]
           + alph[19]*G[23]*G[16]

           +3.*G[31]) ;


 F[31] =1./(3.*a)*(-3.*(2.-hub)*G[31]-3.*hub*G[30]*mua

          -2.*beta[9]*G[17]*G[23]

          -2.*beta[8]*G[1]*G[29]

          -2.*beta[7]*G[13]*G[27]

          -2.*gam2*gamk[10]*G[16]*G[22]
          -2.*gam2*gamk[11]*G[0]*G[28]
          -2.*gam2*gamk[9]*G[12]*G[26]

          -gam3*(gamk[10]*gamk[3] + gamk[11]*gamk[8] + gamk[9]*gamk[5])*G[12]*G[16]*G[0]);



  // F2/G2(k5,k6) = (k3-p,k1)
  alph[20] = alphai(k5,k6,x56);
  alph[21] = alphai(k6,k5,x56);
  beta[10] = betai(k5,k6,x56);

  F[32] =1./a*(-(alph[20]*G[15]*G[16]+alph[21]*G[14]*G[17])/2.-G[33]) ;
  F[33] =1./a*(-(2.-hub)*G[33]-hub*G[32]*mua - gam2*gamk[12]*G[16]*G[14]-beta[10]*G[17]*G[15]);


// F3/G3(k6,k1,k5) = (k1,p,k3-p)

 alph[22] = alphai(k5,k16,x516);
 alph[23] = alphai(k6,k8,x68);

 alph[24] = alphai(k16,k5,x516);
 alph[25] = alphai(k8,k6,x68);

 beta[11] = betai(k5,k16,x516);
 beta[12] = betai(k8,k6,x68);


 F[34] = - 1./(3.*a)*(alph[22]*G[26]*G[15]
           + alph[2]*G[32]*G[1]
           + alph[23]*G[24]*G[17]

           + alph[24]*G[27]*G[14]
           + alph[3]*G[33]*G[0]
           + alph[25]*G[25]*G[16]

           +3.*G[35]) ;


 F[35] =1./(3.*a)*(-3.*(2.-hub)*G[35]-3.*hub*G[34]*mua

          -2.*beta[12]*G[17]*G[25]

          -2.*beta[1]*G[1]*G[33]

          -2.*beta[11]*G[15]*G[27]

          -2*gam2*gamk[14]*G[16]*G[24]
          -2*gam2*gamk[2]*G[0]*G[32]
          -2*gam2*gamk[13]*G[14]*G[26]

            -gam3*(gamk[14]*gamk[4] + gamk[2]*gamk[12] + gamk[13]*gamk[5])*G[14]*G[16]*G[0]);


 // F2/G2(k7,k1) = (k2,p)
 alph[26] = alphai(k7,k1,x7);
 alph[27] = alphai(k1,k7,x7);
 beta[13] = betai(k1,k7,x7);

 F[36] =1./a*(-(alph[26]*G[19]*G[0]+alph[27]*G[1]*G[18])/2.-G[37]) ;
 F[37] =1./a*(-(2.-hub)*G[37]-hub*G[36]*mua - gam2*gamk[6]*G[18]*G[0]-beta[13]*G[19]*G[1]);

 // F2/G2(k5,k7) = (k3-p,k2)
 alph[28] = alphai(k5,k7,x57);
 alph[29] = alphai(k7,k5,x57);
 beta[14] = betai(k5,k7,x57);

 F[38] =1./a*(-(alph[28]*G[15]*G[18]+alph[29]*G[19]*G[14])/2.-G[39]) ;
 F[39] =1./a*(-(2.-hub)*G[39]-hub*G[38]*mua - gam2*gamk[16]*G[14]*G[18]-beta[14]*G[15]*G[19]);


 // F3/G3(k7,k1,k5) = (k2,p,k3-p)

 alph[30] = alphai(k5,k3,x35);
 alph[31] = alphai(k1,k16,x157);
 alph[32] = alphai(k7,k8,x78);

 alph[33] = alphai(k3,k5,x35);
 alph[34] = alphai(k16,k1,x157);
 alph[35] = alphai(k8,k7,x78);

 beta[15] = betai(k5,k3,x35);
 beta[16] = betai(k1,k16,x157);
 beta[17] = betai(k8,k7,x78);



  F[40] = - 1./(3.*a)*(alph[30]*G[36]*G[15]
            + alph[31]*G[38]*G[1]
            + alph[32]*G[24]*G[19]

            + alph[33]*G[37]*G[14]
            + alph[34]*G[39]*G[0]
            + alph[35]*G[25]*G[18]

            +3.*G[41]) ;


  F[41] =1./(3.*a)*(-3.*(2.-hub)*G[41]-3.*hub*G[40]*mua

           -2.*beta[17]*G[19]*G[25]

           -2.*beta[16]*G[1]*G[39]

           -2.*beta[15]*G[15]*G[37]

           -2*gam2*gamk[18]*G[18]*G[24]
           -2*gam2*gamk[19]*G[0]*G[38]
           -2*gam2*gamk[17]*G[14]*G[36]

             -gam3*(gamk[18]*gamk[4] + gamk[19]*gamk[16] + gamk[17]*gamk[6])*G[14]*G[18]*G[0]);


 // F2/G2(k2,k7) = (k1-p,k2)
 alph[36] = alphai(k7,k2,x27);
 alph[37] = alphai(k2,k7,x27);
 beta[18] = betai(k2,k7,x27);

 F[42] =1./a*(-(alph[36]*G[19]*G[2]+alph[37]*G[3]*G[18])/2.-G[43]) ;
 F[43] =1./a*(-(2.-hub)*G[43]-hub*G[42]*mua - gam2*gamk[20]*G[18]*G[2]-beta[18]*G[19]*G[3]);


//F3/G3(k1-p,p,k2)

 F[44] = - 1./(3.*a)*(alph[0]*G[36]*G[5]
           + alph[15]*G[42]*G[1]
           + alph[19]*G[10]*G[19]

           + alph[1]*G[37]*G[4]
           + alph[18]*G[43]*G[0]
           + alph[16]*G[11]*G[18]

           +3.*G[45]) ;


 F[45] =1./(3.*a)*(-3.*(2.-hub)*G[45]-3.*hub*G[44]*mua

          -2.*beta[9]*G[19]*G[11]

          -2.*beta[8]*G[1]*G[43]

          -2.*beta[0]*G[5]*G[37]

          -2*gam2*gamk[0]*G[4]*G[36]
          -2*gam2*gamk[11]*G[0]*G[42]
          -2*gam2*gamk[10]*G[18]*G[10]

            -gam3*(gamk[0]*gamk[6] + gamk[11]*gamk[20] + gamk[10]*gamk[1])*G[4]*G[18]*G[0]);

 // F2/G2(k8,k1) = (k3,p)

 alph[38] = alphai(k8,k1,x8);
 alph[39] = alphai(k1,k8,x8);
 beta[19] = betai(k1,k8,x8);

 F[46] =1./a*(-(alph[38]*G[21]*G[0]+alph[39]*G[1]*G[20])/2.-G[47]) ;
 F[47] =1./a*(-(2.-hub)*G[47]-hub*G[46]*mua - gam2*gamk[7]*G[20]*G[0]-beta[19]*G[21]*G[1]);

 // F2/G2(k4,k8) = (k2-p,k3)
 alph[40] = alphai(k8,k4,x48);
 alph[41] = alphai(k4,k8,x48);
 beta[20] = betai(k4,k8,x48);

 F[48] =1./a*(-(alph[40]*G[21]*G[12]+alph[41]*G[13]*G[20])/2.-G[49]) ;
 F[49] =1./a*(-(2.-hub)*G[49]-hub*G[48]*mua - gam2*gamk[22]*G[20]*G[12]-beta[20]*G[21]*G[13]);




 // F3/G3(k8,k1,k4) = (k3,p,k2-p)
 alph[42] = alphai(k4,k46,x418);
 alph[43] = alphai(k46,k4,x418);
 beta[21] = betai(k4,k46,x418);

 F[50] = - 1./(3.*a)*(
             alph[42]*G[46]*G[13]
           + alph[31]*G[48]*G[1]
           + alph[35]*G[22]*G[21]

           + alph[43]*G[47]*G[12]
           + alph[34]*G[49]*G[0]
           + alph[32]*G[23]*G[20]

           +3.*G[51]) ;


 F[51] =1./(3.*a)*(-3.*(2.-hub)*G[51]-3.*hub*G[50]*mua

          -2.*beta[17]*G[21]*G[23]

          -2.*beta[16]*G[1]*G[49]

          -2.*beta[21]*G[13]*G[47]

          -2*gam2*gamk[18]*G[20]*G[22]
          -2*gam2*gamk[19]*G[0]*G[48]
          -2*gam2*gamk[23]*G[12]*G[46]

            -gam3*(gamk[18]*gamk[3] + gamk[19]*gamk[22] + gamk[23]*gamk[7])*G[12]*G[20]*G[0]);

 // F2/G2(k2,k8) = (k1-p,k3)

 alph[44] = alphai(k8,k2,x28);
 alph[45] = alphai(k2,k8,x28);
 beta[22] = betai(k2,k8,x28);


 F[52] =1./a*(-(alph[44]*G[21]*G[2]+alph[45]*G[3]*G[20])/2.-G[53]) ;
 F[53] =1./a*(-(2.-hub)*G[53]-hub*G[52]*mua - gam2*gamk[24]*G[20]*G[2]-beta[22]*G[21]*G[3]);


// F3/G3(k3,p,k1-p)

 alph[46] = alphai(k2,k46,x218);
 alph[47] = alphai(k46,k2,x218);
 beta[23] = betai(k2,k46,x218);

 F[54] = - 1./(3.*a)*(alph[46]*G[46]*G[3]
           + alph[2]*G[52]*G[1]
           + alph[25]*G[10]*G[21]

           + alph[47]*G[47]*G[2]
           + alph[3]*G[53]*G[0]
           + alph[23]*G[11]*G[20]

           +3.*G[55]) ;


 F[55] =1./(3.*a)*(-3.*(2.-hub)*G[55]-3.*hub*G[54]*mua

          -2.*beta[12]*G[21]*G[11]

          -2.*beta[1]*G[1]*G[53]

          -2.*beta[23]*G[3]*G[47]

          -2*gam2*gamk[14]*G[20]*G[10]
          -2*gam2*gamk[2]*G[0]*G[52]
          -2*gam2*gamk[25]*G[2]*G[46]

            -gam3*(gamk[14]*gamk[1] + gamk[2]*gamk[24] + gamk[25]*gamk[7])*G[2]*G[20]*G[0]);


// F2/G2(k6,k7) = F2/G2(k1,k2)

F[56] =1./a*(-(alph[16]*G[17]*G[18]+alph[19]*G[19]*G[16])/2.-G[57]) ;
F[57] =1./a*(-(2.-hub)*G[57]-hub*G[56]*mua - gam2*gamk[10]*G[18]*G[16]-beta[9]*G[17]*G[19]);


// F2/G2(k6,k8) = F2/G2(k1,k3)

F[58] =1./a*(-(alph[23]*G[17]*G[20]+alph[25]*G[21]*G[16])/2.-G[59]) ;
F[59] =1./a*(-(2.-hub)*G[59]-hub*G[58]*mua - gam2*gamk[14]*G[20]*G[16]-beta[12]*G[17]*G[21]);


// F2/G2(k7,k8) = F2/G2(k2,k3)

F[60] =1./a*(-(alph[32]*G[19]*G[20]+alph[35]*G[21]*G[18])/2.-G[61]) ;
F[61] =1./a*(-(2.-hub)*G[61]-hub*G[60]*mua - gam2*gamk[18]*G[20]*G[18]-beta[17]*G[19]*G[21]);


// F2/G2(-k3,p)

alph[48] = alphai(k8,k1,-x8);
alph[49] = alphai(k1,k8,-x8);
beta[25] = betai(k1,k8,-x8);

F[62] =1./a*(-(alph[48]*G[21]*G[0]+alph[49]*G[1]*G[20])/2.-G[63]) ;
F[63] =1./a*(-(2.-hub)*G[63]-hub*G[62]*mua - gam2*gamk[7]*G[20]*G[0]-beta[25]*G[21]*G[1]);


 // F2/G2(-k2,p)
 alph[50] = alphai(k7,k1,-x7);
 alph[51] = alphai(k1,k7,-x7);
 beta[26] = betai(k1,k7,-x7);

 F[64] =1./a*(-(alph[50]*G[19]*G[0]+alph[51]*G[1]*G[18])/2.-G[65]) ;
 F[65] =1./a*(-(2.-hub)*G[65]-hub*G[64]*mua- gam2*gamk[6]*G[18]*G[0]-beta[26]*G[19]*G[1]);


 alph[52] = alphai(k6,k1,-x6);
 alph[53] = alphai(k1,k6,-x6);
 beta[27] = betai(k1,k6,-x6);

 // F2/G2(k6,k1) = (-k1,p)
 F[66] =1./a*(-(alph[52]*G[17]*G[0]+alph[53]*G[1]*G[16])/2.-G[67]) ;
 F[67] =1./a*(-(2.-hub)*G[67]-hub*G[66]*mua - gam2*gamk[5]*G[16]*G[0]-beta[27]*G[17]*G[1]);


 // F3/G3(k6,k7,k1) = (-k1,-k2,p)

 F[68] = - 1./(3.*a)*(alph[39]*G[56]*G[1]
                    + alph[36]*G[66]*G[19]
                    + alph[12]*G[64]*G[17]

           + alph[38]*G[57]*G[0]
           + alph[37]*G[67]*G[18]
           + alph[13]*G[65]*G[16]

           +3.*G[69]) ;

 F[69] =1./(3.*a)*(-3.*(2.-hub)*G[69]-3.*hub*G[68]*mua

          -2.*beta[6]*G[65]*G[17]

          -2.*beta[18]*G[67]*G[19]

          -2.*beta[19]*G[57]*G[1]

          -2.*gam2*gamk[8]*G[16]*G[64]
          -2.*gam2*gamk[20]*G[18]*G[66]
          -2.*gam2*gamk[7]*G[0]*G[56]

            -gam3*(gamk[8]*gamk[6] + gamk[20]*gamk[5] + gamk[7]*gamk[10])*G[0]*G[18]*G[16]);



// F3/G3(k6,k7,k1) = (-k1,-k2,-p)

alph[54] = alphai(k7,k16,x37);
alph[55] = alphai(k6,k3,x63);

alph[56] = alphai(k16,k7,x37);
alph[57] = alphai(k3,k6,x63);

beta[28] = betai(k7,k16,x37);
beta[29] = betai(k6,k3,x63);


F[70] = - 1./(3.*a)*(alph[49]*G[56]*G[1]
                  + alph[54]*G[26]*G[19]
                  + alph[55]*G[36]*G[17]

         + alph[48]*G[57]*G[0]
         + alph[56]*G[27]*G[18]
         + alph[57]*G[37]*G[16]

         +3.*G[71]) ;


F[71] =1./(3.*a)*(-3.*(2.-hub)*G[71]-3.*hub*G[70]*mua

        -2.*beta[29]*G[37]*G[17]

        -2.*beta[28]*G[27]*G[19]

        -2.*beta[25]*G[57]*G[1]

        -2.*gam2*gamk[26]*G[16]*G[36]
        -2.*gam2*gamk[15]*G[18]*G[26]
        -2.*gam2*gamk[7]*G[0]*G[56]

          -gam3*(gamk[26]*gamk[6] + gamk[15]*gamk[5] + gamk[7]*gamk[10])*G[0]*G[18]*G[16]);


F[72] =1./a*(-G[73]) ; // set 1st term to 0 because x1 is not identically -1 (needed to keep evolution at 0 for GR)
F[73] =1./a*(-(2.-hub)*G[73]-hub*G[72]*mua);


F[74] = - 1./(3.*a)*(alph[5]*G[66]*G[1]
                  + alph[31]*G[26]*G[1]
            //      + alph[58]*G[72]*G[17]

         				+ alph[4]*G[67]*G[0]
         				+ alph[34]*G[27]*G[0]
         		//		+ alph[59]*G[73]*G[16]

         				+3.*G[75]) ;


F[75] =1./(3.*a)*(-3.*(2.-hub)*G[75]-3.*hub*G[74]*mua

    //    -2.*beta[30]*G[73]*G[17]

        -2.*beta[16]*G[27]*G[1]

        -2.*beta[2]*G[67]*G[1]

      //  -2*gam2*0.*G[16]*G[72]
        -2*gam2*gamk[19]*G[0]*G[26]
        -2*gam2*gamk[1]*G[0]*G[66]

          -gam3*(gamk[19]*gamk[5] + gamk[1]*gamk[5])*G[0]*G[0]*G[16]);


// F3/G3(k7,k1,k1) = (-k2,p,-p)

alph[60] = alphai(k7,k1*sqrt(2.*(1.+XMIN)),0.);
alph[61] = alphai(k1*sqrt(2.*(1.+XMIN)),k7,0.);
beta[31] = betai(k1*sqrt(2.*(1.+XMIN)),k7,0.);


F[76] = - 1./(3.*a)*(alph[7]*G[64]*G[1]
                   + alph[2]*G[36]*G[1]
        //           + alph[60]*G[72]*G[19]

                  + alph[6]*G[65]*G[0]
                  + alph[3]*G[37]*G[0]
          //        + alph[61]*G[73]*G[18]

                  +3.*G[77]) ;


F[77] =1./(3.*a)*(-3.*(2.-hub)*G[77]-3.*hub*G[76]*mua

        -2.*beta[3]*G[65]*G[1]

        -2.*beta[1]*G[37]*G[1]

    //    -2.*beta[31]*G[73]*G[19]

      //  -2.*gam2*(1.-pow2(XMIN))*G[18]*G[72]
        -2.*gam2*gamk[2]*G[0]*G[36]
        -2.*gam2*gamk[3]*G[0]*G[64]

        -gam3*(gamk[2]*gamk[6] + gamk[3]*gamk[6] + (1.-pow2(XMIN)))*G[0]*G[0]*G[18]);


// F4/G4(-k6,-k7,k1,-k1) = (-k1,-k2,p,-p)

alph[62] = 1.;//alphai(1.-XMAX,k8,0.);
alph[63] = 1.; //alphai(k8,1.-XMAX,0.);
beta[24] = 0.;// betai(1.-XMAX,k8,0.);




F[78] = -1./(24.*a)*(
										6.*alph[16]*G[76]*G[17] // a(-k1,-k2) F3(-k2,p,-p) G1(k1)
                   +6.*alph[19]*G[74]*G[19] // a(-k2,-k1) F3(-k1,p,-p) G1(k2)
                   +6.*alph[9]*G[70]*G[1] // a(p,k3-p) F3(-p,-k1,-k2) G1(p)
                   +6.*alph[15]*G[68]*G[1] // a(-p,k3+p) F3(p,-k1,-k2) G1(p)

                    +6.*alph[19]*G[77]*G[16] // a(-k2,-k1) G3(-k2,p,-p) F1(k1)
                    +6.*alph[16]*G[75]*G[18] // a(-k2,-k1) G3(-k1,p,-p) F1(k2)
                    +6.*alph[8]*G[71]*G[0] // a(p,k3-p) G3(-p,-k1,-k2) F1(p)
                    +6.*alph[18]*G[69]*G[0] // a(p,k3-p) G3(-p,-k1,-k2) F1(p)

                //    +4.*alph[62]*G[56]*G[73] // a(0,k3) F2(k1,k2) G2(p,-p)
                //    +4.*alph[63]*G[72]*G[57] // a(k3,0) F2(p,-p) G2(k1,k2)

                    + 4.*alph[17]*G[64]*G[27] // a(-k1-p,-k2+p,) F2(-k2,p) G2(-p,-k1)
                    + 4.*alph[14]*G[26]*G[65] // a(-k2+p,-k1-p) F2(-k1,-p) G2(-k2,p)

                    + 4.*alph[0]*G[36]*G[67] // a(-k1+p,-k2-p) F2(-k2,-p)G2(-k1,p)
                    + 4.*alph[1]*G[66]*G[37] // a(-k2-p,-k1+p) F2(-k1,p) G2(-k2,-p)

                    +24.*G[79]);

F[79] = 1./a*(-(2.-hub)*G[79]-hub*G[78]*mua

          -1./24.*(
                  12.*beta[9]*G[77]*G[17] // b(-k1,-k2) G3(-k2,p,-p) G1(k1)
                  +12.*beta[9]*G[75]*G[19] // b(-k2,-k1) G3(-k1,p,-p) G1(k2)
                  +12.*beta[4]*G[71]*G[1] // b(p,k3-p) G3(-p,-k1,-k2) G1(p)
                  +12.*beta[8]*G[69]*G[1] // b(-p,k3+p) G3(p,-k1,-k2) G1(p)

              //    + 8.*beta[24]*G[57]*G[73] // b(k3,0) G2(k1,k2) G2(p,-p)
                  + 8.*beta[7]*G[65]*G[27]  // b(k2-p,k1+p) G2(-k2,p) G2(k1,p)
                  + 8.*beta[0]*G[67]*G[37] // b(k2+p,k1-p) G2(k2,p) G2(-k1,p)

//2nd order
            //      + 8.*gam2*(1.-pow2(XMIN))*G[56]*G[72]
                  + 8.*gam2*gamk[9]*G[64]*G[26]
                  + 8.*gam2*gamk[0]*G[66]*G[36]

//3rd order
                  + 12.*gam2*gamk[10]*G[76]*G[16]
                  + 12.*gam2*gamk[10]*G[74]*G[18]
                  + 12.*gam2*gamk[4]*G[70]*G[0]
                  + 12.*gam2*gamk[11]*G[68]*G[0]


									+4.*gam3*((gamk[5]*gamk[9] + gamk[10]*gamk[3] + gamk[8]*gamk[11])*G[16]*G[0]*G[64]
										 			 +(gamk[6]*gamk[9] + gamk[10]*gamk[19] + gamk[15]*gamk[4])*G[18]*G[0]*G[26]
												   +(gamk[0]*gamk[6] + gamk[10]*gamk[1] + gamk[11]*gamk[20])*G[18]*G[0]*G[66]
												 	 +(gamk[0]*gamk[5] + gamk[10]*gamk[2] + gamk[4]*gamk[26])*G[16]*G[0]*G[36]
												 	 +(gamk[10]+ gamk[10] + gamk[10])*G[16]*G[18]*G[72]
												   +(gamk[7]*gamk[4] + gamk[7]*gamk[11])*G[0]*G[0]*G[56])

// 4th order
									 +  gam4*4.*(gamk[5]*gamk[6]*gamk[0]
										 				 + gamk[5]*gamk[6]*gamk[9]

										 				   + gamk[10]*gamk[6]*gamk[2]
															 + gamk[4]*gamk[6]*gamk[26]
															 + gamk[10]*gamk[6]*gamk[3]
															 + gamk[11]*gamk[6]*gamk[8]
															 + gamk[10]*gamk[5]*gamk[19]
															 + gamk[4]*gamk[5]*gamk[15]
															 + gamk[10]*gamk[5]*gamk[1]
															 + gamk[11]*gamk[5]*gamk[20]
															 + gamk[4]*gamk[10]*gamk[7]
															 + gamk[11]*gamk[10]*gamk[7])
															   *G[16]*G[18]*G[0]*G[0]));

// F3/G3(k6,k8,k1) = (-k1,-k3,p)

F[80] = - 1./(3.*a)*(alph[21]*G[62]*G[17] // a(k1,k3-p) F2(-k3,p) G1(k1)
                   + alph[44]*G[66]*G[21] // a(k3,k1-p) F2(-k1,p) G1(k3)
                   + alph[27]*G[58]*G[1] // a(p,k2) F2(k1,k3) G1(p)

          + alph[20]*G[63]*G[16]
          + alph[45]*G[67]*G[20]
          + alph[26]*G[59]*G[0]

          +3.*G[81]) ;

F[81] =1./(3.*a)*(-3.*(2.-hub)*G[81]-3.*hub*G[80]*mua

         -2.*beta[10]*G[63]*G[17]

         -2.*beta[22]*G[67]*G[21]

         -2.*beta[13]*G[59]*G[1]

         -2*gam2*gamk[12]*G[16]*G[62]
         -2*gam2*gamk[24]*G[20]*G[66]
         -2*gam2*gamk[6]*G[0]*G[58]

           -gam3*(gamk[12]*gamk[7] + gamk[24]*gamk[5] + gamk[6]*gamk[14])*G[0]*G[20]*G[16]);


// F3/G3(k8,k1,k1) = (-k3,p,-p)


F[82] = - 1./(3.*a)*(alph[9]*G[62]*G[1] // a(p,k3-p) F2(k3,-p) G1(p)
                 + alph[15]*G[46]*G[1] // a(p,-k3-p) F2(k3,p) G1(p)
          //       + alph[63]*G[72]*G[20] // a(-k3,0.) F2(p,-p) G1(k3)

        + alph[8]*G[63]*G[0]
        + alph[18]*G[47]*G[0]
    //    + alph[62]*G[73]*G[21]

        +3.*G[83]) ;


F[83] =1./(3.*a)*(-3.*(2.-hub)*G[83]-3.*hub*G[82]*mua

       -2.*beta[4]*G[63]*G[1]

       -2.*beta[8]*G[47]*G[1]

       -2.*beta[24]*G[73]*G[21]


  //     -2*gam2*0.*G[20]*G[72]
       -2*gam2*gamk[11]*G[0]*G[46]
       -2*gam2*gamk[4]*G[0]*G[62]

         -gam3*(gamk[11]*gamk[7] + gamk[4]*gamk[7])*G[0]*G[0]*G[20]);


// F3/G3(k6,k8,k1) = (-k1,-k3,-p)

alph[64] = alphai(k6,k46, x646); // a(k1,k3+p)
alph[65] = alphai(k8,k16,x816); // a(k3,k1+p)

alph[66] = alphai(k46,k6,x646); // a(k3+p,k1)
alph[67] = alphai(k16,k8,x816); // a(k1+p,k3)

beta[32] = betai(k46,k6,x646); // a(k3+p,k1)
beta[33] = betai(k8,k16,x816); // a(k3+p,k3)

F[84] = - 1./(3.*a)*(alph[64]*G[46]*G[17] // a(k1,k3+p) F2(k3,p) G1(k1)
                  + alph[65]*G[26]*G[21] // a(k3,k1+p) F2(k1,p) G1(k3)
                  + alph[51]*G[58]*G[1] // a(p,-k2) F2(k1,k3) G1(p)

         + alph[66]*G[47]*G[16]
         + alph[67]*G[27]*G[20]
         + alph[50]*G[59]*G[0]

         +3.*G[85]) ;

F[85] =1./(3.*a)*(-3.*(2.-hub)*G[85]-3.*hub*G[84]*mua

        -2.*beta[32]*G[47]*G[17]

        -2.*beta[33]*G[27]*G[21]

        -2.*beta[26]*G[59]*G[1]

			 -2*gam2*gamk[29]*G[16]*G[46]
       -2*gam2*gamk[27]*G[20]*G[26]
			 -2*gam2*gamk[6]*G[0]*G[58]

         -gam3*(gamk[29]*gamk[7] + gamk[27]*gamk[5]  + gamk[6]*gamk[14])*G[0]*G[20]*G[16]);


// F4/G4(-k6,-k8,k1,-k1) = (-k1,-k3,p,-p)


alph[68] = alphai(k46,k2,x218);
alph[69] = alphai(k2,k46,x218);
beta[35] = betai(k46,k2,x218);


alph[70] = alphai(k1,k3,x3);
alph[71] = alphai(k3,k1,x3);
beta[34] = betai(k3,k1,x3);



F[86] = -1./(24.*a)*(6.*alph[23]*G[82]*G[17] // a(k1,k3) F3(-k3,p,-p) G1(k1)
                    +6.*alph[25]*G[74]*G[21] // a(k3,k1) F3(-k1,p,-p) G1(k3)
                    +6.*alph[7]*G[84]*G[1] // a(p,k2-p) F3(-p,-k1,-k3) G1(p)
                    +6.*alph[70]*G[80]*G[1] // a(-p,k2+p) F3(p,-k1,-k3) G1(p)

                   +6.*alph[25]*G[83]*G[16]
                   +6.*alph[23]*G[75]*G[20]
                   +6.*alph[6]*G[85]*G[0]
                   +6.*alph[71]*G[81]*G[0]



          //         + 4.*alph[60]*G[58]*G[73] // a(0,k3) F2(k1,k3) G2(p,-p)
          //         + 4.*alph[61]*G[72]*G[59] // a(k3,0) F2(p,-p) G2(k1,k3)

                   + 4.*alph[24]*G[62]*G[27] // a(k1+p,k3-p) F2(-k3,p) G2(-p,-k1)
                   + 4.*alph[22]*G[26]*G[63] // a(k3-p,k1+p) F2(-k1,-p) G2(-k3,p)

                   + 4.*alph[69]*G[46]*G[67] // a(k1-p,k3+p) F2(-k3,-p)G2(-k1,p)
                   + 4.*alph[68]*G[66]*G[47] // a(k3+p,k1-p) F2(-k1,p) G2(-k3,-p)

                   +24.*G[87]);

F[87] = 1./a*(-(2.-hub)*G[87]-hub*G[86]*mua

            -1./24.*(
                  12.*beta[12]*G[83]*G[17] // b(k1,k3) G3(-k3,p,-p) G1(k1)
                 +12.*beta[12]*G[75]*G[21] // b(k3,k1) G3(-k1,p,-p) G1(k3)
                 +12.*beta[3]*G[85]*G[1] // b(p,k2-p) G3(-p,-k1,-k3) G1(p)
                 +12.*beta[34]*G[81]*G[1]  // b(-p,k2+p) G3(p,-k1,-k3) G1(p)

                 + 8.*beta[31]*G[59]*G[73] // b(k2,0) G2(k1,k3) G2(p,-p)
                 + 8.*beta[11]*G[63]*G[27]  // b(k3-p,k1+p) G2(-k3,p) G2(k1,p)
                 + 8.*beta[35]*G[67]*G[47] // b(k3+p,k1-p) G2(k3,p) G2(-k1,p)

            //     + 8.*gam2*0.*G[58]*G[72]
                 + 8.*gam2*gamk[13]*G[62]*G[26]
                 + 8.*gam2*gamk[25]*G[66]*G[46]

                 + 12.*gam2*gamk[14]*G[82]*G[16]
                 + 12.*gam2*gamk[14]*G[74]*G[20]
                 + 12.*gam2*gamk[3]*G[84]*G[0]
                 + 12.*gam2*gamk[2]*G[80]*G[0]

								 + 4.*gam3*((gamk[13]*gamk[5] + gamk[14]*gamk[4] + gamk[2]*gamk[12])*G[16]*G[0]*G[62]
													+(gamk[13]*gamk[7] + gamk[14]*gamk[19] + gamk[27]*gamk[3])*G[20]*G[0]*G[26]
													+(gamk[7]*gamk[25] + gamk[14]*gamk[1] + gamk[2]*gamk[24])*G[20]*G[0]*G[66]
													+(gamk[5]*gamk[25] + gamk[14]*gamk[11] + gamk[3]*gamk[29])*G[16]*G[0]*G[46]
													+(gamk[14]+ gamk[14] + gamk[14])*G[16]*G[20]*G[72]
													+(gamk[6]*gamk[3] + gamk[6]*gamk[2])*G[0]*G[0]*G[58])


								 +  gam4*4.*(gamk[5]*gamk[7]*gamk[25]
									 		     + gamk[5]*gamk[7]*gamk[13]

									 				   + gamk[14]*gamk[7]*gamk[11]
														 + gamk[3]*gamk[7]*gamk[29]
														 + gamk[14]*gamk[7]*gamk[4]
														 + gamk[2]*gamk[7]*gamk[12]
														 + gamk[14]*gamk[5]*gamk[19]
														 + gamk[3]*gamk[5]*gamk[27]
														 + gamk[14]*gamk[5]*gamk[1]
														 + gamk[2]*gamk[5]*gamk[24]
														 + gamk[3]*gamk[14]*gamk[6]
														 + gamk[2]*gamk[14]*gamk[6])
														 *G[16]*G[20]*G[0]*G[0]));




// F3/G3(k6,k8,k1) = (-k2,-k3,-p)

alph[72] = alphai(k7,k46,-x746); // a(k2,k3+p)
alph[73] = alphai(k46,k7,-x746);

alph[74] = alphai(k8,k3,x38); // a(k3,k2+p)
alph[75] = alphai(k3,k8,x38);

beta[36] = betai(k7,k46,-x746);
beta[37] = betai(k8,k3,x38);


F[88] = - 1./(3.*a)*(alph[72]*G[46]*G[19] // a(k2,k3+p) F2(k3,p) G1(k2)
                 + alph[74]*G[36]*G[21] // a(k3,k2+p) F2(k2,p) G1(k3)
                 + alph[53]*G[60]*G[1] // a(p,-k1) F2(k2,k3) G1(p)


        + alph[73]*G[47]*G[18]
        + alph[75]*G[37]*G[20]
        + alph[52]*G[61]*G[0]

        +3.*G[89]) ;


F[89] =1./(3.*a)*(-3.*(2.-hub)*G[89]-3.*hub*G[88]*mua

       -2.*beta[36]*G[47]*G[19]

       -2.*beta[37]*G[37]*G[21]

       -2.*beta[27]*G[61]*G[1]


      -2*gam2*gamk[21]*G[18]*G[46]
      -2*gam2*gamk[28]*G[20]*G[36]
      -2*gam2*gamk[5]*G[0]*G[60]

         -gam3*(gamk[21]*gamk[7] + gamk[28]*gamk[6] + gamk[5]*gamk[18])*G[0]*G[20]*G[18]);



// F3/G3(k6,k8,k1) = (-k2,-k3,p)

 F[90] = - 1./(3.*a)*(alph[29]*G[62]*G[19] // a(k2,k3-p) F2(k3,-p) G1(k2)
                  + alph[40]*G[64]*G[21] // a(k3,k2-p) F2(k2,-p) G1(k3)
                  + alph[11]*G[60]*G[1] // a(p,k1) F2(k2,k3) G1(p)

         + alph[28]*G[63]*G[18]
         + alph[41]*G[65]*G[20]
         + alph[10]*G[61]*G[0]

         +3.*G[91]) ;



 F[91] =1./(3.*a)*(-3.*(2.-hub)*G[91]-3.*hub*G[90]*mua

        -2.*beta[14]*G[63]*G[19]

        -2.*beta[20]*G[65]*G[21]

        -2.*beta[5]*G[61]*G[1]

        -2*gam2*gamk[16]*G[18]*G[62]
        -2*gam2*gamk[22]*G[20]*G[64]
        -2*gam2*gamk[5]*G[0]*G[60]

          -gam3*(gamk[16]*gamk[7] + gamk[22]*gamk[6] + gamk[5]*gamk[18])*G[0]*G[20]*G[18]);


// F4/G4(-k6,-k8,k1,-k1) = (-k2,-k3,p,-p)

F[92] = -1./(24.*a)*(6.*alph[32]*G[82]*G[19] // a(k2,k3) F3(-k3,p,-p) G1(k2)
                   +6.*alph[35]*G[76]*G[21] // a(k3,k2) F3(-k2,p,-p) G1(k3)
                   +6.*alph[5]*G[88]*G[1] // a(p,k1-p) F3(-p,-k2,-k3) G1(p)
                   +6.*alph[31]*G[90]*G[1] // a(-p,k1+p) F3(p,-k2,-k3) G1(p)

                  +6.*alph[35]*G[83]*G[18]
                  +6.*alph[32]*G[77]*G[20]
                  +6.*alph[4]*G[89]*G[0]
                  +6.*alph[34]*G[91]*G[0]


            //      + 4.*alph[58]*G[60]*G[73] // a(0,k1) F2(k2,k3) G2(p,-p)
            //      + 4.*alph[59]*G[72]*G[61] // a(k1,0) F2(p,-p) G2(k2,k3)
                  + 4.*alph[33]*G[62]*G[37] // a(k2+p,k3-p) F2(-k3,p) G2(p,k2)
                  + 4.*alph[30]*G[36]*G[63] // a(k3-p,k2+p) F2(k2,p) G2(-k3,p)

                  + 4.*alph[42]*G[46]*G[65] // a(k2-p,k3+p) F2(-k3,-p)G2(-k2,p)
                  + 4.*alph[43]*G[64]*G[47] // a(k3+p,k2-p) F2(-k2,p) G2(-k3,-p)


                  +24.*G[93]);

F[93] = 1./a*(-(2.-hub)*G[93]-hub*G[92]*mua

        -1./24.*(
                 12.*beta[17]*G[83]*G[19] // b(k2,k3) G3(-k3,p,-p) G1(k2)
                +12.*beta[17]*G[77]*G[21] // b(k3,k2) G3(-k2,p,-p) G1(k3)
                +12.*beta[2]*G[89]*G[1]// b(p,k1-p) G3(-p,-k2,-k3) G1(p)
                +12.*beta[16]*G[91]*G[1]  // b(-p,k1+p) G3(p,-k2,-k3) G1(p)

          //      + 8.*beta[30]*G[61]*G[73] // b(k1,0) G2(k2,k3) G2(p,-p)
                + 8.*beta[15]*G[63]*G[37]  // b(k3-p,k2+p) G2(-k3,p) G2(k2,p)
                + 8.*beta[21]*G[65]*G[47] // b(k3+p,k2-p) G2(k3,p) G2(-k2,p)

        //        + 8.*gam2*0.*G[60]*G[72]
                + 8.*gam2*gamk[17]*G[62]*G[36]
                + 8.*gam2*gamk[23]*G[64]*G[46]

                + 12.*gam2*gamk[18]*G[82]*G[18]
                + 12.*gam2*gamk[18]*G[76]*G[20]
                + 12.*gam2*gamk[1]*G[88]*G[0]
                + 12.*gam2*gamk[19]*G[90]*G[0]

								+4.*gam3*((gamk[6]*gamk[17] + gamk[18]*gamk[4] + gamk[16]*gamk[19])*G[18]*G[0]*G[62]
												 +(gamk[17]*gamk[7] + gamk[18]*gamk[2] + gamk[28]*gamk[1])*G[20]*G[0]*G[36]
												 +(gamk[7]*gamk[23] + gamk[18]*gamk[3] + gamk[19]*gamk[22])*G[20]*G[0]*G[64]
												 +(gamk[6]*gamk[23] + gamk[18]*gamk[11] + gamk[1]*gamk[21])*G[18]*G[0]*G[46]
												 +(gamk[18]+ gamk[18] + gamk[18])*G[18]*G[20]*G[72]
												 +(gamk[5]*gamk[1] + gamk[5]*gamk[19])*G[0]*G[0]*G[60])

							  +  gam4*4.*(gamk[6]*gamk[7]*gamk[23] + gamk[6]*gamk[7]*gamk[17]

								 				   + gamk[18]*gamk[7]*gamk[11] + gamk[1]*gamk[7]*gamk[21]
													 + gamk[18]*gamk[7]*gamk[4] + gamk[19]*gamk[7]*gamk[16]
													 + gamk[18]*gamk[6]*gamk[2] + gamk[1]*gamk[6]*gamk[2]
													 + gamk[18]*gamk[6]*gamk[3] + gamk[19]*gamk[6]*gamk[22]
													 + gamk[1]*gamk[18]*gamk[5] + gamk[19]*gamk[18]*gamk[5])
													 *G[18]*G[20]*G[0]*G[0]));



  	return GSL_SUCCESS;
  }


// Bispectrum B411, B321 AND B222 initialisations
// k array holds vector magnitudes
// x array holds angle between vectors
// kargs holds magnitudes of combinations of vectors
void BSPTN::initnb1_dgp(double A, double k[], double x[], double kargs[], double omega0, double par1, double par2, double par3, double epars[])
{
				double a = 0.0001;

       double G[94] = { a,-a,a,-a,a,-a, 0.,0.,0.,0.,0.,0.,a,-a,a,-a,a,-a,a,-a,a,-a,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                            0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

			/*Parameters passed to system of equations */
          struct param_type5 my_params1 = {k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7],x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],
                                           kargs[0],kargs[1],kargs[2],kargs[3],kargs[4],kargs[5],kargs[6],kargs[7],kargs[8],
                                           kargs[9],kargs[10],kargs[11],kargs[12],kargs[13],kargs[14],kargs[15],kargs[16],kargs[17],
                                           kargs[18],kargs[19],kargs[20], kargs[21], kargs[22],kargs[23],kargs[24], omega0, par1, par2, par3};

          gsl_odeiv2_system sys = {funcdgp, jacb, 94, &my_params1};

			//  this gives the system, the method, the initial interval step, the absolute error and the relative error.
				gsl_odeiv2_driver * d =
				gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
													epars[0], epars[1] ,epars[2]);
// smallest possible accuracy and initial step when comparing to analytic result
// must reduce initial step size for earlier times!
// rk8pd is optimal

				int status1 = gsl_odeiv2_driver_apply (d, &a, A, G);


				/*Allocation of array values */

// For tree level initialisation

  /*1st order */

			//F1(k;a), G1(k;a)
			F1b_k1 = G[0] ;
			G1b_k1 = G[1] ;

			// F1(k3;a), G1(k3;a)
			F1b_k3 = G[2];
			G1b_k3 = G[3];

			//F1(k2;a), G1(k2;a)
			F1b_k2 = G[4];
			G1b_k2 = G[5];


// for B321 and B222 initialisaiton

      /* 1st order */

      //F1(k1;a), G1(k1;a)
			F1b_1 = G[16] ;
			G1b_1 = G[17] ;

			// F1(k2;a), G1(k2;a)
			F1b_2 = G[18];
			G1b_2 = G[19];

			//F1(k3;a), G1(k3;a)
			F1b_3 = G[20];
			G1b_3 = G[21];


			/*2nd order*/

			//F2/G2(k3,k2) for tree,  F2/G2(k1-p,k2+p) for B222
			F2b_k23 =  G[6];
			G2b_k23 =  G[7];

      //F2/G2(k1,k3) for tree, F2/G2(-p,k2+p) for B222
      F2b_k13 =  G[8];
      G2b_k13 =  G[9];

      //F2/G2(k1,k2) for tree, F2/G2(p,k1-p)
      F2b_k12 =  G[10];
      G2b_k12 =  G[11];


      // F2/G2(p,k2-p)
      F2b_p2mp = G[22];
      G2b_p2mp = G[23];

      // F2/G2(p,k3-p)
      F2b_p3mp = G[24];
      G2b_p3mp = G[25];

      // F3/G3(k1,p,k2-p)
      F3b_12mp = G[30];
      G3b_12mp = G[31];

      // F3/G3(k1,p,k3-p)
      F3b_13mp = G[34];
      G3b_13mp = G[35];

      // F3/G3(k2,p,k3-p)
      F3b_23mp = G[40];
      G3b_23mp = G[41];

      // F3/G3(k2,p,k1-p)
      F3b_21mp = G[44];
      G3b_21mp = G[45];


      // F3/G3(k3,p,k2-p)
      F3b_32mp = G[50];
      G3b_32mp = G[51];

      // F3/G3(k3,p,k1-p)
      F3b_31mp = G[54];
      G3b_31mp = G[55];

      // F3/G3(k1,k2) for B321-II
      F2b_12a = G[56];
      G2b_12a = G[57];

      // F3/G3(k1,k3) for B321-II
      F2b_13a = G[58];
      G2b_13a = G[59];

      // F3/G3(k2,k3) for B321-II
      F2b_23a = G[60];
      G2b_23a = G[61];

      //F3/G3(k1,p,-p)
      F3b_1pp = G[74];
      G3b_1pp = G[75];

      //F3/G3(k2,p,-p)
      F3b_2pp = G[76];
      G3b_2pp = G[77];

      //F3/G3(k3,p,-p)
      F3b_3pp = G[82];
      G3b_3pp = G[83];

      // F4/G4(-k1,-k2,p,-p)

      F4b_12pp = G[78];
      G4b_12pp = G[79];

      // F4/G4(-k1,-k3,p,-p)

      F4b_13pp = G[86];
      G4b_13pp = G[87];

      // F4/G4(-k2,-k3,p,-p)

      F4b_23pp = G[92];
      G4b_23pp = G[93];

			gsl_odeiv2_driver_free(d);
}




struct param_type6 {
  real kk1;
	real kk2;
	real kk3;
  real kk4;
  real kk5;
  real kk6;
  real kk7;
  real kk8;

	real xx1;
	real xx2;
	real xx3;
	real xx4;
	real xx5;
	real xx6;
	real xx7;
	real xx8;

  real arg1;
  real arg2;
  real arg3;
  real arg4;
  real arg5;
  real arg6;
  real arg7;
  real arg8;
  real arg9;
  real arg10;
  real arg11;
  real arg12;
  real arg13;
  real arg14;
  real arg15;
  real arg16;
  real arg17;
  real arg18;
  real arg19;
  real arg20;
  real arg21;
  real arg22;
  real arg23;
  real arg24;
  real arg25;

		real k1s;
		real k2s;
		real k3s;
		real k4s;
		real k5s;
		real k6s;
		real k7s;
		real k8s;
		real k9s;
		real k10s;
		real k11s;


  real omega00;
	real par1;
	real par2;
	real par3;


};


int funcfr(double a, const double G[], double F[], void *params)
{
	param_type6 p = *(param_type6 *)(params);
  // variables

// magnitudes
	real k1 = p.kk1;  // p
	real k2 = p.kk2; //  k1-p
	real k3 = p.kk3; // k2+p
  real k4 = p.kk4; // k2-p
	real k5 = p.kk5; // k3-p
  real k6 = p.kk6; // k1
  real k7 = p.kk7; // k2
  real k8 = p.kk8; // k3 (= -k1-k2)
//	real k13 = p.arg1; // k2 + 2p
  real k16 = p.arg2; // k1 + p
  real k46 = p.arg3; // -k3-p


// angles
	real x1= p.xx1;  // k1-p.k2+p
	real x2= p.xx2;  // p.k1-p
	real x3= p.xx3; // -p.k2+p
	real x4= p.xx4; // p.k2-p
	real x5= p.xx5; // p.k3-p
  real x6= p.xx6; // p.k1
  real x7= p.xx7; // p.k2
  real x8= p.xx8; // p.k3

  real x46 = p.arg4; // k2-p . k1
  real x416 = p.arg5; // k2-p . p+k1
  real x67 = p.arg6; // k1.k2 = x
  real x146 = p.arg7; // p. -k3-p
  real x56 = p.arg8; // k3-p. k1
  real x516 = p.arg9; // k3-p. k1+p
  real x68 = p.arg10; // k1.k3
  real x37 = p.arg11; // k2.k1+p
  real x57 = p.arg12; // k3-p.k2
  real x35 = p.arg13; // k2+p.k3-p
  real x78 = p.arg14; // k2.k3
  real x157 = p.arg15; // - p . k1+p
  real x27 = p.arg16; // k1-p.k2
  real x746 = p.arg17; // k2. -k3-p
  real x48 = p.arg18; // k2-p.k3
  real x418 = p.arg19; // k2-p.p+k3
  real x28 = p.arg20; // k1-p . k3
  real x218 = p.arg21; // k1-p.p+k3
  real x63 = p.arg22; // k1.k2+p
  real x816 = p.arg23; // k3.k1+p
  real x38 = p.arg24; // k2+p.k3
  real x646 = p.arg25; // k1.k3+p

	double alph[76];
  double beta[38];
	double sqrs[11];
	double sqrsa[11];
	double pifs[11];
	double muak[11];

// Cosmological and gravity parameters
 	real omega0 = p.omega00;
	real p1 = p.par1;
//	real p2 = p.par2; // unused
//	real p3 = p.par3; // unused

	sqrs[0] = p.k1s;
	sqrs[1] = p.k2s;
	sqrs[2] = p.k3s;
	sqrs[3] = p.k4s;
	sqrs[4] = p.k5s;
	sqrs[5] = p.k6s;
	sqrs[6] = p.k7s;
	sqrs[7] = p.k8s;
	sqrs[8] = p.k9s;
	sqrs[9] = p.k10s;
	sqrs[10] = p.k11s;

// time dependent factors
double hub = HA2(a,omega0);
double hub1= HA(a,omega0);
double hubsqr = pow2(hub1);
double asqr = pow2(a);
double acub = asqr*a;

double epsi = (omega0+4.*acub*(1.-omega0))/acub;
double epsi2 = pow2(epsi);
double epsi3 = epsi*epsi2;
double epsi5 = epsi2*epsi3;
double epsi6 = epsi3*epsi3;
double epsi7 = epsi2*epsi5;
double epsi9 = epsi2*epsi7;


double p1p = p1/myh0sqr * pow2(3.*omega0-4.);
double p1p2 = pow2(p1p);
double p1p3 = p1p*p1p2;
double p1p4 = p1p2*p1p2;

double om0a2 = pow2(omega0/acub);
double om0a3 = om0a2*omega0/acub;;
double om0a4 = pow2(om0a2);

double term1 = epsi3/2./p1p;

//  mu and PI
for(int i =0; i<11; i++){
	sqrsa[i] = sqrs[i]/asqr;
	pifs[i] = sqrsa[i] + term1;
	muak[i] = 1.+ sqrsa[i]/ 3. / pifs[i];
	// redefine to reduce computations in gamma functions (sqrsa always divided by pifs)
	sqrsa[i] = sqrs[i]/asqr / pifs[i];
}

double gam2 = -3./16./hubsqr*om0a2*epsi5/p1p2;
double gam3a = -5./32./hubsqr*om0a3*epsi7/p1p3;
double gam3b = -9./10.*gam3a*epsi3/p1p;

double gam4a = -35./256/hubsqr*om0a4*epsi9/p1p4;
double gam4b = gam4a*27./140.*epsi6/p1p2;
double gam4c = gam4a*9./7.*epsi3/p1p;
double gam4d = gam4b*6.;

double gamk[36];
//gamma2, gamma3 and gamma4 scale dependencies  (sorry but need to put them by hand for optimisation)

// arguments commented on right
// optimise computations - lots of repetitions!
gamk[21] = sqrsa[1]/(pifs[6]*pifs[10]); // k2.k3+p
gamk[28] = sqrsa[1]/(pifs[7]*pifs[2]);  //k3.k2+p
gamk[30] = sqrsa[1]/(pifs[0]*pifs[5]); // -p.k1

gamk[6] = sqrsa[2]/(pifs[0]*pifs[6]); // p.k2
gamk[12] = sqrsa[2]/(pifs[5]*pifs[4]); // k1.k3-p
gamk[24] = sqrsa[2]/(pifs[7]*pifs[1]); // k3.k1-p

gamk[27] = sqrsa[3]/(pifs[7]*pifs[9]);  // k3.k1+p
gamk[29] = sqrsa[3]/(pifs[5]*pifs[10]);  //k1.k3+p
gamk[31] = sqrsa[3]/(pifs[0]*pifs[6]); // -p.k2

gamk[15] = sqrsa[4]/(pifs[6]*pifs[9]); // k2.k1+p
gamk[26] = sqrsa[4]/(pifs[5]*pifs[2]); // k1.k2+p
gamk[32] = sqrsa[4]/(pifs[0]*pifs[7]); // -p.k3

gamk[1] = sqrsa[5]/(pifs[0]*pifs[1]); // p.k1-p
gamk[17] = sqrsa[5]/(pifs[2]*pifs[4]);// k2+p.k3-p
gamk[18] = sqrsa[5]/(pifs[6]*pifs[7]);// k2.k3
gamk[19] = sqrsa[5]/(pifs[0]*pifs[9]);// -p.k1+p
gamk[23] = sqrsa[5]/(pifs[10]*pifs[3]);  // k3+p.k2-p
gamk[33] = sqrsa[5]/(pifs[5]*term1); // k1,0

gamk[2] = sqrsa[6]/(pifs[0]*pifs[2]); // -p.k2+p
gamk[3] = sqrsa[6]/(pifs[0]*pifs[3]); // p. k2-p
gamk[13] = sqrsa[6]/(pifs[9]*pifs[4]);// k1+p.k3-p
gamk[14] = sqrsa[6]/(pifs[5]*pifs[7]); // k1.k3
gamk[25] = sqrsa[6]/(pifs[10]*pifs[3]);  // k3+p.k1-p
gamk[34] = sqrsa[6]/(pifs[6]*term1); // k2,0

gamk[0] = sqrsa[7]/(pifs[1]*pifs[2]); //k1-p.k2+p
gamk[4] = sqrsa[7]/(pifs[0]*pifs[4]); // p.k3-p
gamk[9] = sqrsa[7]/(pifs[9]*pifs[3]); // p+k1.k2-p
gamk[10] = sqrsa[7]/(pifs[5]*pifs[6]); // k1.k2
gamk[11] = sqrsa[7]/(pifs[0]*pifs[10]); // -p.k3+p
gamk[35] = sqrsa[7]/(pifs[7]*term1); // k3,0

gamk[5] = sqrsa[9]/(pifs[0]*pifs[5]); // p.k1
gamk[22] = sqrsa[9]/(pifs[7]*pifs[3]); // k3.k2-p
gamk[16] = sqrsa[9]/(pifs[6]*pifs[4]);// k2.k3-p

gamk[7] = sqrsa[10]/(pifs[0]*pifs[7]); // p.k3
gamk[8] = sqrsa[10]/(pifs[5]*pifs[3]); // k1.k2-p
gamk[20] = sqrsa[10]/(pifs[6]*pifs[1]); // k2.k1-p


double gamk3[27];
gamk3[13] = sqrsa[1]/(pifs[6]*pifs[7]*pifs[0]); // -k2,-k3,-p
gamk3[10] = sqrsa[2]/(pifs[5]*pifs[7]*pifs[0]); //-k1,-k3,p
gamk3[12] = sqrsa[3]/(pifs[5]*pifs[7]*pifs[0]); // -k1,-k3,-p
gamk3[7] = sqrsa[4]/(pifs[5]*pifs[6]*pifs[0]); // -k1,-k2,-p


gamk3[2] = sqrsa[5]/(pifs[6]*pifs[0]*pifs[4]); // k2,p,k3-p
gamk3[4] = sqrsa[5]/(pifs[7]*pifs[0]*pifs[3]); // k3,p,k2-p
gamk3[8] = sqrsa[5]/(pifs[5]*pifs[0]*pifs[0]); //-k1,p,-p
gamk3[19] = sqrsa[5]/(pifs[10]*pifs[6]*pifs[0]); // -k3-p,-k2,p
gamk3[20] = sqrsa[5]/(pifs[2]*pifs[7]*pifs[0]); // k2+p,-k3,p
gamk3[25] = sqrsa[5]/(pifs[0]*pifs[0]*pifs[5]); // k1,p,-p
gamk3[26] = sqrsa[5]/(term1*pifs[6]*pifs[7]); // 0.,-k2,-k3

gamk3[1] = sqrsa[6]/(pifs[5]*pifs[0]*pifs[4]); // k1,p,k3-p
gamk3[5] = sqrsa[6]/(pifs[7]*pifs[0]*pifs[1]); // k3,p,k1-p
gamk3[17] = sqrsa[6]/(pifs[9]*pifs[7]*pifs[0]); // -k1-p,-k3,p
gamk3[18] = sqrsa[6]/(pifs[10]*pifs[5]*pifs[0]); // -k3-p,-k1,p
gamk3[23] = sqrsa[6]/(pifs[0]*pifs[0]*pifs[6]); // k2,p,-p
gamk3[24] = sqrsa[6]/(term1*pifs[5]*pifs[7]); // 0.,-k1,-k3
gamk3[9] = sqrsa[6]/(pifs[6]*pifs[0]*pifs[0]); // -k2,p,-p


gamk3[0] = sqrsa[7]/(pifs[5]*pifs[0]*pifs[3]); // k1,p,k2-p
gamk3[3] = sqrsa[7]/(pifs[1]*pifs[0]*pifs[6]); // k1-p,p,k2
gamk3[11] = sqrsa[7]/(pifs[7]*pifs[0]*pifs[0]); // -k3,p,-p
gamk3[15] = sqrsa[7]/(pifs[9]*pifs[6]*pifs[0]); // -k1-p,-k2,p
gamk3[16] = sqrsa[7]/(pifs[2]*pifs[5]*pifs[0]);// -k2-p,-k1,p
gamk3[21] = sqrsa[7]/(pifs[0]*pifs[0]*pifs[7]); // k3,p,-p
gamk3[22] = sqrsa[7]/(term1*pifs[6]*pifs[5]); // 0.,-k1,-k2

gamk3[14] = sqrsa[9]/(pifs[6]*pifs[7]*pifs[0]); // -k2,-k3,p
gamk3[6] = sqrsa[10]/(pifs[5]*pifs[6]*pifs[0]); // -k1,-k2,p

double gamk4[3];
gamk4[0] = sqrsa[7]/(pifs[5]*pifs[6]*pifs[0]*pifs[0]); // -k1,-k2,p,-p
gamk4[1] = sqrsa[6]/(pifs[5]*pifs[7]*pifs[0]*pifs[0]); // -k1,-k3,p,-p
gamk4[2] = sqrsa[5]/(pifs[7]*pifs[6]*pifs[0]*pifs[0]); // -k2,-k3,p,-p


	/* 1st order */
	//1. F1/G1(p)
	F[0] = -G[1]/a;
	F[1] =1./a*(-(2.-hub)*G[1]-hub*G[0]*muak[0]);

	//2.  F1/G1(k2+p)
	F[2] = -G[3]/a;
	F[3] =1./a*(-(2.-hub)*G[3]-hub*G[2]*muak[2]);

	//3. F1/G1(k1-p)
	F[4] = -G[5]/a;
	F[5] =1./a*(-(2.-hub)*G[5]-hub*G[4]*muak[1]);

	/* 2nd order for B222 */
  alph[0] = alphai(k2,k3,x1);
  alph[1] = alphai(k3,k2,x1);
  beta[0] = betai(k2,k3,x1);

	//4. F2/G2(k1-p,k2+p)
	F[6] =1./a*(-(alph[0]*G[5]*G[2]+alph[1]*G[3]*G[4])/2.-G[7]);
	F[7] =1./a*(-(2.-hub)*G[7]-hub*G[6]*muak[7] - gam2*gamk[0]*G[4]*G[2] - beta[0]*G[5]*G[3]);

  alph[2] = alphai(k1,k3,x3);
  alph[3] = alphai(k3,k1,x3);
  beta[1] = betai(k1,k3,x3);

	//5. F2/G2(-p,k2+p)
	F[8] =1./a*(-(alph[2]*G[1]*G[2]+alph[3]*G[3]*G[0])/2.-G[9]);
	F[9] =1./a*(-(2.-hub)*G[9]-hub*G[8]*muak[6] - gam2*gamk[2]*G[2]*G[0] - beta[1]*G[3]*G[1]);

  alph[4] = alphai(k2,k1,x2);
  alph[5] = alphai(k1,k2,x2);
  beta[2] = betai(k1,k2,x2);

	//7. F2/G2(k2,k1) or F2/G2(k1-p,p)
	F[10] =1./a*(-(alph[4]*G[5]*G[0]+alph[5]*G[1]*G[4])/2.-G[11]) ;
	F[11] =1./a*(-(2.-hub)*G[11]-hub*G[10]*muak[5] - gam2*gamk[1]*G[4]*G[0]-beta[2]*G[5]*G[1]);


  /* 1st order for B321 - I */

  //  F1/G1(k2-p)

  F[12] = -G[13]/a;
  F[13] =1./a*(-(2.-hub)*G[13]-hub*G[12]*muak[3]);

  // F1/G1(k3-p)

  F[14] = -G[15]/a;
  F[15] =1./a*(-(2.-hub)*G[15]-hub*G[14]*muak[4]);

  //  F1/G1(k1)

  F[16] = -G[17]/a;
  F[17] =1./a*(-(2.-hub)*G[17]-hub*G[16]*muak[5]);

  //  F1/G1(k2)

  F[18] = -G[19]/a;
  F[19] =1./a*(-(2.-hub)*G[19]-hub*G[18]*muak[6]);

  //  F1/G1(k3)

  F[20] = -G[21]/a;
  F[21] =1./a*(-(2.-hub)*G[21]-hub*G[20]*muak[7]);

  // 2nd order used in 3rd order equations of B321-I

  alph[6] = alphai(k4,k1,x4);
  alph[7] = alphai(k1,k4,x4);
  beta[3] = betai(k1,k4,x4);


  // F2/G2(k2-p,p)
  F[22] =1./a*(-(alph[6]*G[13]*G[0]+alph[7]*G[1]*G[12])/2.-G[23]) ;
  F[23] =1./a*(-(2.-hub)*G[23]-hub*G[22]*muak[6] - gam2*gamk[3]*G[12]*G[0]-beta[3]*G[13]*G[1]);

  alph[8] = alphai(k5,k1,x5);
  alph[9] = alphai(k1,k5,x5);
  beta[4] = betai(k1,k5,x5);

  // F2/G2(k5,k1) = (k3-p,p)
  F[24] =1./a*(-(alph[8]*G[15]*G[0]+alph[9]*G[1]*G[14])/2.-G[25]) ;
  F[25] =1./a*(-(2.-hub)*G[25]-hub*G[24]*muak[7] - gam2*gamk[4]*G[14]*G[0]-beta[4]*G[15]*G[1]);


  alph[10] = alphai(k6,k1,x6);
  alph[11] = alphai(k1,k6,x6);
  beta[5] = betai(k1,k6,x6);

  // F2/G2(k6,k1) = (k1,p)
  F[26] =1./a*(-(alph[10]*G[17]*G[0]+alph[11]*G[1]*G[16])/2.-G[27]) ;
  F[27] =1./a*(-(2.-hub)*G[27]-hub*G[26]*muak[9] - gam2*gamk[5]*G[16]*G[0]-beta[5]*G[17]*G[1]);

  alph[12] = alphai(k6,k4,x46);
  alph[13] = alphai(k4,k6,x46);
  beta[6] = betai(k4,k6,x46);

  // F2/G2(k6,k4) = (k1,k2-p)
  F[28] =1./a*(-(alph[12]*G[17]*G[12]+alph[13]*G[13]*G[16])/2.-G[29]) ;
  F[29] =1./a*(-(2.-hub)*G[29]-hub*G[28]*muak[10] - gam2*gamk[8]*G[16]*G[12]-beta[6]*G[17]*G[13]);




  // F3/G3(k6,k1,k4) = (k1,p,k2-p)

 alph[14] = alphai(k4,k16,x416);
 alph[15] = alphai(k1,k46,x146);
 alph[16] = alphai(k6,k7,x67);

 alph[17] = alphai(k16,k4,x416);
 alph[18] = alphai(k46,k1,x146);
 alph[19] = alphai(k7,k6,x67);

 beta[7] = betai(k4,k16,x416);
 beta[8] = betai(k1,k46,x146);
 beta[9] = betai(k7,k6,x67);

 F[30] = - 1./(3.*a)*(alph[14]*G[26]*G[13]
           + alph[15]*G[28]*G[1]
           + alph[16]*G[22]*G[17]

           + alph[17]*G[27]*G[12]
           + alph[18]*G[29]*G[0]
           + alph[19]*G[23]*G[16]

           +3.*G[31]) ;


 F[31] =1./(3.*a)*(-3.*(2.-hub)*G[31]-3.*hub*G[30]*muak[7]

          -2.*beta[9]*G[17]*G[23]

          -2.*beta[8]*G[1]*G[29]

          -2.*beta[7]*G[13]*G[27]

          -2.*gam2*gamk[10]*G[16]*G[22] // gam2(k1,p + k2- p )
          -2.*gam2*gamk[11]*G[0]*G[28] // gam2(p, k1 + k2-p )
          -2.*gam2*gamk[9]*G[12]*G[26] // gam2(k2-p, k1 + p)

					-gamk3[0]*(3.*gam3a + gam3b*(1./pifs[10] + 1./pifs[9] + 1./pifs[6]))*G[12]*G[16]*G[0]);


  // F2/G2(k3-p,k1)
  alph[20] = alphai(k5,k6,x56);
  alph[21] = alphai(k6,k5,x56);
  beta[10] = betai(k5,k6,x56);

  F[32] =1./a*(-(alph[20]*G[15]*G[16]+alph[21]*G[14]*G[17])/2.-G[33]) ;
  F[33] =1./a*(-(2.-hub)*G[33]-hub*G[32]*muak[2] - gam2*gamk[12]*G[16]*G[14]-beta[10]*G[17]*G[15]);


// F3/G3(k6,k1,k5) = (k1,p,k3-p)

 alph[22] = alphai(k5,k16,x516);
 alph[23] = alphai(k6,k8,x68);

 alph[24] = alphai(k16,k5,x516);
 alph[25] = alphai(k8,k6,x68);

 beta[11] = betai(k5,k16,x516);
 beta[12] = betai(k8,k6,x68);


 F[34] = - 1./(3.*a)*(alph[22]*G[26]*G[15]
           + alph[2]*G[32]*G[1]
           + alph[23]*G[24]*G[17]

           + alph[24]*G[27]*G[14]
           + alph[3]*G[33]*G[0]
           + alph[25]*G[25]*G[16]

           +3.*G[35]) ;


 F[35] =1./(3.*a)*(-3.*(2.-hub)*G[35]-3.*hub*G[34]*muak[6]

          -2.*beta[12]*G[17]*G[25]

          -2.*beta[1]*G[1]*G[33]

          -2.*beta[11]*G[15]*G[27]

          -2.*gam2*gamk[14]*G[16]*G[24] // gam2(k1, p + k3-p)
          -2.*gam2*gamk[2]*G[0]*G[32] // gam2(p, k3-p + k1)
          -2.*gam2*gamk[13]*G[14]*G[26] // gam2( k3-p, p + k1)

				-gamk3[1]*(3.*gam3a + gam3b*(1./pifs[9] + 1./pifs[7] + 1./pifs[2]))*G[14]*G[16]*G[0]);


 // F2/G2(k7,k1) = (k2,p)
 alph[26] = alphai(k7,k1,x7);
 alph[27] = alphai(k1,k7,x7);
 beta[13] = betai(k1,k7,x7);

 F[36] =1./a*(-(alph[26]*G[19]*G[0]+alph[27]*G[1]*G[18])/2.-G[37]) ;
 F[37] =1./a*(-(2.-hub)*G[37]-hub*G[36]*muak[2] - gam2*gamk[6]*G[18]*G[0]-beta[13]*G[19]*G[1]);


 // F2/G2(k5,k7) = (k3-p,k2)
 alph[28] = alphai(k5,k7,x57);
 alph[29] = alphai(k7,k5,x57);
 beta[14] = betai(k5,k7,x57);

 F[38] =1./a*(-(alph[28]*G[15]*G[18]+alph[29]*G[19]*G[14])/2.-G[39]) ;
 F[39] =1./a*(-(2.-hub)*G[39]-hub*G[38]*muak[9] - gam2*gamk[16]*G[14]*G[18]-beta[14]*G[15]*G[19]);


 // F3/G3(k7,k1,k5) = (k2,p,k3-p)

 alph[30] = alphai(k5,k3,x35);
 alph[31] = alphai(k1,k16,x157);
 alph[32] = alphai(k7,k8,x78);

 alph[33] = alphai(k3,k5,x35);
 alph[34] = alphai(k16,k1,x157);
 alph[35] = alphai(k8,k7,x78);

 beta[15] = betai(k5,k3,x35);
 beta[16] = betai(k1,k16,x157);
 beta[17] = betai(k8,k7,x78);



  F[40] = - 1./(3.*a)*(alph[30]*G[36]*G[15]
            + alph[31]*G[38]*G[1]
            + alph[32]*G[24]*G[19]

            + alph[33]*G[37]*G[14]
            + alph[34]*G[39]*G[0]
            + alph[35]*G[25]*G[18]

            +3.*G[41]) ;



  F[41] =1./(3.*a)*(-3.*(2.-hub)*G[41]-3.*hub*G[40]*muak[5]

           -2.*beta[17]*G[19]*G[25]

           -2.*beta[16]*G[1]*G[39]

           -2.*beta[15]*G[15]*G[37]

           -2.*gam2*gamk[18]*G[18]*G[24] // gam2(k2, p + k3-p)
           -2.*gam2*gamk[19]*G[0]*G[38]  // gam2(p, k3-p + k2)
           -2.*gam2*gamk[17]*G[14]*G[36]  // gam2(k3-p, k2 + p)

				-gamk3[2]*(3.*gam3a + gam3b*(1./pifs[2] + 1./pifs[7] + 1./pifs[9]))*G[14]*G[18]*G[0]);


 // F2/G2(k2,k7) = (k1-p,k2)
 alph[36] = alphai(k7,k2,x27);
 alph[37] = alphai(k2,k7,x27);
 beta[18] = betai(k2,k7,x27);

 F[42] =1./a*(-(alph[36]*G[19]*G[2]+alph[37]*G[3]*G[18])/2.-G[43]) ;
 F[43] =1./a*(-(2.-hub)*G[43]-hub*G[42]*muak[10] - gam2*gamk[20]*G[18]*G[2]-beta[18]*G[19]*G[3]);


//F3/G3(k1-p,p,k2)

 F[44] = - 1./(3.*a)*(alph[0]*G[36]*G[5]
           + alph[15]*G[42]*G[1]
           + alph[19]*G[10]*G[19]

           + alph[1]*G[37]*G[4]
           + alph[18]*G[43]*G[0]
           + alph[16]*G[11]*G[18]

           +3.*G[45]) ;


 F[45] =1./(3.*a)*(-3.*(2.-hub)*G[45]-3.*hub*G[44]*muak[7]

          -2.*beta[9]*G[19]*G[11]

          -2.*beta[8]*G[1]*G[43]

          -2.*beta[0]*G[5]*G[37]


          -2*gam2*gamk[0]*G[4]*G[36] // gam2(k1-p, k2, p  )
          -2*gam2*gamk[11]*G[0]*G[42] // gam2(p, k1-p + k2 )
          -2*gam2*gamk[10]*G[18]*G[10] // gam2(k2, k1-p + p  )

				-gamk3[3]*(3.*gam3a + gam3b*(1./pifs[5] + 1./pifs[10] + 1./pifs[2]))*G[4]*G[18]*G[0]);


 // F2/G2(k8,k1) = (k3,p)

 alph[38] = alphai(k8,k1,x8);
 alph[39] = alphai(k1,k8,x8);
 beta[19] = betai(k1,k8,x8);

 F[46] =1./a*(-(alph[38]*G[21]*G[0]+alph[39]*G[1]*G[20])/2.-G[47]) ;
 F[47] =1./a*(-(2.-hub)*G[47]-hub*G[46]*muak[10] - gam2*gamk[7]*G[20]*G[0]-beta[19]*G[21]*G[1]);

 // F2/G2(k4,k8) = (k2-p,k3)
 alph[40] = alphai(k8,k4,x48);
 alph[41] = alphai(k4,k8,x48);
 beta[20] = betai(k4,k8,x48);

 F[48] =1./a*(-(alph[40]*G[21]*G[12]+alph[41]*G[13]*G[20])/2.-G[49]) ;
 F[49] =1./a*(-(2.-hub)*G[49]-hub*G[48]*muak[9] - gam2*gamk[22]*G[20]*G[12]-beta[20]*G[21]*G[13]);



 // F3/G3(k8,k1,k4) = (k3,p,k2-p)
 alph[42] = alphai(k4,k46,x418);
 alph[43] = alphai(k46,k4,x418);
 beta[21] = betai(k4,k46,x418);

 F[50] = - 1./(3.*a)*(
             alph[42]*G[46]*G[13]
           + alph[31]*G[48]*G[1]
           + alph[35]*G[22]*G[21]

           + alph[43]*G[47]*G[12]
           + alph[34]*G[49]*G[0]
           + alph[32]*G[23]*G[20]

           +3.*G[51]) ;


 F[51] =1./(3.*a)*(-3.*(2.-hub)*G[51]-3.*hub*G[50]*muak[5]

          -2.*beta[17]*G[21]*G[23]

          -2.*beta[16]*G[1]*G[49]

          -2.*beta[21]*G[13]*G[47]

          -2.*gam2*gamk[18]*G[20]*G[22] // gam2(k3, p + k2-p)
          -2.*gam2*gamk[19]*G[0]*G[48] // gam2(p, k3 + k2-p)
          -2.*gam2*gamk[23]*G[12]*G[46] // gam2(k2-p, k3 + p)

					-gamk3[4]*(3.*gam3a + gam3b*(1./pifs[10] + 1./pifs[6] + 1./pifs[9]))*G[12]*G[20]*G[0]);

 // F2/G2(k2,k8) = (k1-p,k3)


 alph[44] = alphai(k8,k2,x28);
 alph[45] = alphai(k2,k8,x28);
 beta[22] = betai(k2,k8,x28);


 F[52] =1./a*(-(alph[44]*G[21]*G[2]+alph[45]*G[3]*G[20])/2.-G[53]) ;
 F[53] =1./a*(-(2.-hub)*G[53]-hub*G[52]*muak[2] - gam2*gamk[24]*G[20]*G[2]-beta[22]*G[21]*G[3]);


// F3/G3(k3,p,k1-p)

 alph[46] = alphai(k2,k46,x218);
 alph[47] = alphai(k46,k2,x218);
 beta[23] = betai(k2,k46,x218);

 F[54] = - 1./(3.*a)*(alph[46]*G[46]*G[3]
           + alph[2]*G[52]*G[1]
           + alph[25]*G[10]*G[21]

           + alph[47]*G[47]*G[2]
           + alph[3]*G[53]*G[0]
           + alph[23]*G[11]*G[20]

           +3.*G[55]) ;


 F[55] =1./(3.*a)*(-3.*(2.-hub)*G[55]-3.*hub*G[54]*muak[6]

          -2.*beta[12]*G[21]*G[11]

          -2.*beta[1]*G[1]*G[53]

          -2.*beta[23]*G[3]*G[47]

          -2*gam2*gamk[14]*G[20]*G[10] // gam2(k3, p + k1-p)
          -2*gam2*gamk[2]*G[0]*G[52] // gam2(p, k3 + k1-p)
          -2*gam2*gamk[25]*G[2]*G[46] // gam2(k1-p, k3 + p)


				-gamk3[5]*(3.*gam3a + gam3b*(1./pifs[10] + 1./pifs[2] + 1./pifs[5]))*G[2]*G[20]*G[0]);



// F2/G2(k6,k7) = F2/G2(k1,k2)

F[56] =1./a*(-(alph[16]*G[17]*G[18]+alph[19]*G[19]*G[16])/2.-G[57]) ;
F[57] =1./a*(-(2.-hub)*G[57]-hub*G[56]*muak[7] - gam2*gamk[10]*G[18]*G[16]-beta[9]*G[17]*G[19]);


// F2/G2(k6,k8) = F2/G2(k1,k3)

F[58] =1./a*(-(alph[23]*G[17]*G[20]+alph[25]*G[21]*G[16])/2.-G[59]) ;
F[59] =1./a*(-(2.-hub)*G[59]-hub*G[58]*muak[6] - gam2*gamk[14]*G[20]*G[16]-beta[12]*G[17]*G[21]);


// F2/G2(k7,k8) = F2/G2(k2,k3)

F[60] =1./a*(-(alph[32]*G[19]*G[20]+alph[35]*G[21]*G[18])/2.-G[61]) ;
F[61] =1./a*(-(2.-hub)*G[61]-hub*G[60]*muak[5] - gam2*gamk[18]*G[20]*G[18]-beta[17]*G[19]*G[21]);


// F2/G2(-k3,p)

alph[48] = alphai(k8,k1,-x8);
alph[49] = alphai(k1,k8,-x8);
beta[25] = betai(k1,k8,-x8);

F[62] =1./a*(-(alph[48]*G[21]*G[0]+alph[49]*G[1]*G[20])/2.-G[63]) ;
F[63] =1./a*(-(2.-hub)*G[63]-hub*G[62]*muak[4] - gam2*gamk[32]*G[20]*G[0]-beta[25]*G[21]*G[1]);


 // F2/G2(-k2,p)
 alph[50] = alphai(k7,k1,-x7);
 alph[51] = alphai(k1,k7,-x7);
 beta[26] = betai(k1,k7,-x7);

 F[64] =1./a*(-(alph[50]*G[19]*G[0]+alph[51]*G[1]*G[18])/2.-G[65]) ;
 F[65] =1./a*(-(2.-hub)*G[65]-hub*G[64]*muak[3]- gam2*gamk[31]*G[18]*G[0]-beta[26]*G[19]*G[1]);


 alph[52] = alphai(k6,k1,-x6);
 alph[53] = alphai(k1,k6,-x6);
 beta[27] = betai(k1,k6,-x6);

 // F2/G2(k6,k1) = (-k1,p)
 F[66] =1./a*(-(alph[52]*G[17]*G[0]+alph[53]*G[1]*G[16])/2.-G[67]) ;
 F[67] =1./a*(-(2.-hub)*G[67]-hub*G[66]*muak[1] - gam2*gamk[30]*G[16]*G[0]-beta[27]*G[17]*G[1]);


 // F3/G3(k6,k7,k1) = (-k1,-k2,p)

 F[68] = - 1./(3.*a)*(alph[39]*G[56]*G[1]
                    + alph[36]*G[66]*G[19]
                    + alph[12]*G[64]*G[17]

           + alph[38]*G[57]*G[0]
           + alph[37]*G[67]*G[18]
           + alph[13]*G[65]*G[16]

           +3.*G[69]) ;



 F[69] =1./(3.*a)*(-3.*(2.-hub)*G[69]-3.*hub*G[68]*muak[10]

          -2.*beta[6]*G[65]*G[17]

          -2.*beta[18]*G[67]*G[19]

          -2.*beta[19]*G[57]*G[1]

          -2.*gam2*gamk[8]*G[16]*G[64]
          -2.*gam2*gamk[20]*G[18]*G[66]
          -2.*gam2*gamk[7]*G[0]*G[56]

					-gamk3[6]*(3.*gam3a + gam3b*(1./pifs[7] + 1./pifs[1] + 1./pifs[3]))*G[0]*G[18]*G[16]);


// F3/G3(k6,k7,k1) = (-k1,-k2,-p)

alph[54] = alphai(k7,k16,x37);
alph[55] = alphai(k6,k3,x63);

alph[56] = alphai(k16,k7,x37);
alph[57] = alphai(k3,k6,x63);

beta[28] = betai(k7,k16,x37);
beta[29] = betai(k6,k3,x63);


F[70] = - 1./(3.*a)*(alph[49]*G[56]*G[1]
                  + alph[54]*G[26]*G[19]
                  + alph[55]*G[36]*G[17]

         + alph[48]*G[57]*G[0]
         + alph[56]*G[27]*G[18]
         + alph[57]*G[37]*G[16]

         +3.*G[71]) ;


F[71] =1./(3.*a)*(-3.*(2.-hub)*G[71]-3.*hub*G[70]*muak[4]

        -2.*beta[29]*G[37]*G[17]

        -2.*beta[28]*G[27]*G[19]

        -2.*beta[25]*G[57]*G[1]

        -2.*gam2*gamk[26]*G[16]*G[36]
        -2.*gam2*gamk[15]*G[18]*G[26]
        -2.*gam2*gamk[7]*G[0]*G[56]

				-gamk3[7]*(3.*gam3a + gam3b*(1./pifs[7] + 1./pifs[9] + 1./pifs[2]))*G[0]*G[18]*G[16]);


F[72] =1./a*(-G[73]) ; // set 1st term to 0 because x1 is not identically -1 (needed to keep evolution at 0 for GR)
F[73] =1./a*(-(2.-hub)*G[73]-hub*G[72] - betai(k1,k1,-0.99999999)*G[1]*G[1]);


// F3/G3(k6,k1,k1) = (-k1,p,-p)

alph[58] = alphai(k6,1.-0.99999999,0.);
alph[59] = alphai(1.-0.99999999,k6,0.);
beta[30] = betai(1.-0.99999999,k6,0.);

F[74] = - 1./(3.*a)*(alph[5]*G[66]*G[1]
                  + alph[31]*G[26]*G[1]
                  + alph[58]*G[72]*G[17]

         				+ alph[4]*G[67]*G[0]
         				+ alph[34]*G[27]*G[0]
         				+ alph[59]*G[73]*G[16]

         				+3.*G[75]) ;


F[75] =1./(3.*a)*(-3.*(2.-hub)*G[75]-3.*hub*G[74]*muak[5]

        -2.*beta[30]*G[73]*G[17]

        -2.*beta[16]*G[27]*G[1]

        -2.*beta[2]*G[67]*G[1]

        -2*gam2*gamk[33]*G[16]*G[72]  //gam2 (k1,0)
        -2*gam2*gamk[19]*G[0]*G[26]
        -2*gam2*gamk[1]*G[0]*G[66]

				-gamk3[8]*(3.*gam3a + gam3b*(1./pifs[1] + 1./pifs[9] + 1./term1))*G[0]*G[0]*G[16]);



// F3/G3(k7,k1,k1) = (-k2,p,-p)

alph[60] = alphai(k7,1.-0.99999999,0.);
alph[61] = alphai(1.-0.99999999,k7,0.);
beta[31] = betai(1.-0.99999999,k7,0.);


F[76] = - 1./(3.*a)*(alph[7]*G[64]*G[1]
                   + alph[2]*G[36]*G[1]
                   + alph[60]*G[72]*G[19]

                  + alph[6]*G[65]*G[0]
                  + alph[3]*G[37]*G[0]
                  + alph[61]*G[73]*G[18]

                  +3.*G[77]) ;


F[77] =1./(3.*a)*(-3.*(2.-hub)*G[77]-3.*hub*G[76]*muak[6]

        -2.*beta[3]*G[65]*G[1]

        -2.*beta[1]*G[37]*G[1]

        -2.*beta[31]*G[73]*G[19]

        -2.*gam2*gamk[34]*G[18]*G[72]
        -2.*gam2*gamk[2]*G[0]*G[36]
        -2.*gam2*gamk[3]*G[0]*G[64]

				-gamk3[9]*(3.*gam3a + gam3b*(1./pifs[3] + 1./pifs[2] + 1./term1))*G[0]*G[0]*G[18]);



// F4/G4(-k6,-k7,k1,-k1) = (-k1,-k2,p,-p)

alph[62] =  alphai(1.-0.99999999,k8,0.);
alph[63] =  alphai(k8,1.-0.99999999,0.);
beta[24] =  betai(1.-0.99999999,k8,0.);



F[78] = -1./(24.*a)*(
										6.*alph[16]*G[76]*G[17] // a(-k1,-k2) F3(-k2,p,-p) G1(k1)
                   +6.*alph[19]*G[74]*G[19] // a(-k2,-k1) F3(-k1,p,-p) G1(k2)
                   +6.*alph[9]*G[70]*G[1] // a(p,k3-p) F3(-p,-k1,-k2) G1(p)
                   +6.*alph[15]*G[68]*G[1] // a(-p,k3+p) F3(p,-k1,-k2) G1(p)

                    +6.*alph[19]*G[77]*G[16] // a(-k2,-k1) G3(-k2,p,-p) F1(k1)
                    +6.*alph[16]*G[75]*G[18] // a(-k2,-k1) G3(-k1,p,-p) F1(k2)
                    +6.*alph[8]*G[71]*G[0] // a(p,k3-p) G3(-p,-k1,-k2) F1(p)
                    +6.*alph[18]*G[69]*G[0] // a(p,k3-p) G3(-p,-k1,-k2) F1(p)

                    +4.*alph[62]*G[56]*G[73] // a(0,k3) F2(k1,k2) G2(p,-p)
                    +4.*alph[63]*G[72]*G[57] // a(k3,0) F2(p,-p) G2(k1,k2)

                    + 4.*alph[17]*G[64]*G[27] // a(-k1-p,-k2+p,) F2(-k2,p) G2(-p,-k1)
                    + 4.*alph[14]*G[26]*G[65] // a(-k2+p,-k1-p) F2(-k1,-p) G2(-k2,p)

                    + 4.*alph[0]*G[36]*G[67] // a(-k1+p,-k2-p) F2(-k2,-p)G2(-k1,p)
                    + 4.*alph[1]*G[66]*G[37] // a(-k2-p,-k1+p) F2(-k1,p) G2(-k2,-p)

                    +24.*G[79]);

F[79] = 1./a*(-(2.-hub)*G[79]-hub*G[78]*muak[7]

          -1./24.*(
                  12.*beta[9]*G[77]*G[17] // b(-k1,-k2) G3(-k2,p,-p) G1(k1)
                  +12.*beta[9]*G[75]*G[19] // b(-k2,-k1) G3(-k1,p,-p) G1(k2)
                  +12.*beta[4]*G[71]*G[1] // b(p,k3-p) G3(-p,-k1,-k2) G1(p)
                  +12.*beta[8]*G[69]*G[1] // b(-p,k3+p) G3(p,-k1,-k2) G1(p)

                  + 8.*beta[24]*G[57]*G[73] // b(k3,0) G2(k1,k2) G2(p,-p)
                  + 8.*beta[7]*G[65]*G[27]  // b(k2-p,k1+p) G2(-k2,p) G2(k1,p)
                  + 8.*beta[0]*G[67]*G[37] // b(k2+p,k1-p) G2(k2,p) G2(-k1,p)

//2nd order
                  + 8.*gam2*gamk[35]*G[56]*G[72] // gam2(k3,0)
                  + 8.*gam2*gamk[9]*G[64]*G[26]  // gam2(k2-p,k1+p)
                  + 8.*gam2*gamk[0]*G[66]*G[36]  // gam2(k1-p,k2+p)


                  + 12.*gam2*gamk[10]*G[76]*G[16] // gam2(-k2+p-p,-k1)
                  + 12.*gam2*gamk[10]*G[74]*G[18] // gam2(-k1+p-p,-k2)
                  + 12.*gam2*gamk[4]*G[70]*G[0] // gam2(k3-p,p)
                  + 12.*gam2*gamk[11]*G[68]*G[0] // gam2(k3+p,-p)

//3rd order

							+4.*(gamk3[0]*(3.*gam3a + gam3b*(1./pifs[9] + 1./pifs[6] + 1./pifs[10]))*G[16]*G[0]*G[64]
									   + gamk3[3]*(3.*gam3a + gam3b*(1./pifs[2] + 1./pifs[5] + 1./pifs[10]))*G[18]*G[0]*G[26]
										 + gamk3[15]*(3.*gam3a + gam3b*(1./pifs[3] + 1./pifs[5] + 1./pifs[4]))*G[18]*G[0]*G[66]
										 + gamk3[16]*(3.*gam3a + gam3b*(1./pifs[1] + 1./pifs[6] + 1./pifs[4]))*G[16]*G[0]*G[36]
									   + gamk3[21]*(3.*gam3a + gam3b*(1./term1 + 1./pifs[10] + 1./pifs[4]))*G[0]*G[0]*G[56]
								     + gamk3[22]*(3.*gam3a + gam3b*(1./pifs[5] + 1./pifs[6] + 1./pifs[7]))*G[16]*G[18]*G[72])


//4th order
							+ gamk4[0]*(24.*gam4a + 8.*gam4b*(1./(pifs[7]*term1) + 1./(pifs[1]*pifs[2]) + 1./(pifs[9]*pifs[3]))
																		+ 4.*gam4c*(1./pifs[7] + 1./term1 + 1./pifs[1] + 1./pifs[2] + 1./pifs[9] + 1./pifs[3])
																		+ 4./3.*gam4d*(1./(pifs[4]*pifs[7]) + 1./(pifs[4]*pifs[2]) + 1./(pifs[4]*pifs[9])
																			 			   + 1./(pifs[10]*pifs[7]) + 1./(pifs[10]*pifs[1]) + 1./(pifs[10]*pifs[3])
																							 + 1./(pifs[5]*term1) + 1./(pifs[5]*pifs[1]) +	1./(pifs[5]*pifs[9])
																							 + 1./(pifs[6]*term1) + 1./(pifs[6]*pifs[2]) + 1./(pifs[6]*pifs[3]))
																		+ 12./3.*gam4c*(1./pifs[10] + 1./pifs[4] + 1./pifs[5]  + 1./pifs[6]))*G[16]*G[18]*G[0]*G[0]));


// F3/G3(k6,k8,k1) = (-k1,-k3,p)

F[80] = - 1./(3.*a)*(alph[21]*G[62]*G[17] // a(k1,k3-p) F2(-k3,p) G1(k1)
                   + alph[44]*G[66]*G[21] // a(k3,k1-p) F2(-k1,p) G1(k3)
                   + alph[27]*G[58]*G[1] // a(p,k2) F2(k1,k3) G1(p)

          + alph[20]*G[63]*G[16]
          + alph[45]*G[67]*G[20]
          + alph[26]*G[59]*G[0]

          +3.*G[81]) ;


F[81] =1./(3.*a)*(-3.*(2.-hub)*G[81]-3.*hub*G[80]*muak[2]

         -2.*beta[10]*G[63]*G[17]

         -2.*beta[22]*G[67]*G[21]

         -2.*beta[13]*G[59]*G[1]

         -2.*gam2*gamk[12]*G[16]*G[62]
         -2.*gam2*gamk[24]*G[20]*G[66]
         -2.*gam2*gamk[6]*G[0]*G[58]

				 -gamk3[10]*(3.*gam3a + gam3b*(1./pifs[6] + 1./pifs[2] + 1./pifs[4]))*G[0]*G[20]*G[16]);



// F3/G3(k8,k1,k1) = (-k3,p,-p)


F[82] = - 1./(3.*a)*(alph[9]*G[62]*G[1] // a(p,k3-p) F2(k3,-p) G1(p)
                 + alph[15]*G[46]*G[1] // a(p,-k3-p) F2(k3,p) G1(p)
                 + alph[63]*G[72]*G[20] // a(-k3,0.) F2(p,-p) G1(k3)

        + alph[8]*G[63]*G[0]
        + alph[18]*G[47]*G[0]
        + alph[62]*G[73]*G[21]

        +3.*G[83]) ;


F[83] =1./(3.*a)*(-3.*(2.-hub)*G[83]-3.*hub*G[82]*muak[7]

       -2.*beta[4]*G[63]*G[1]

       -2.*beta[8]*G[47]*G[1]

       -2.*beta[24]*G[73]*G[21]


       -2*gam2*gamk[35]*G[20]*G[72]
       -2.*gam2*gamk[11]*G[0]*G[46]
       -2.*gam2*gamk[4]*G[0]*G[62]

			 -gamk3[11]*(3.*gam3a + gam3b*(1./pifs[4] + 1./pifs[10] + 1./term1))*G[0]*G[20]*G[0]);


// F3/G3(k6,k8,k1) = (-k1,-k3,-p)

alph[64] = alphai(k6,k46, x646); // a(k1,k3+p)
alph[65] = alphai(k8,k16,x816); // a(k3,k1+p)

alph[66] = alphai(k46,k6,x646); // a(k3+p,k1)
alph[67] = alphai(k16,k8,x816); // a(k1+p,k3)

beta[32] = betai(k46,k6,x646); // a(k3+p,k1)
beta[33] = betai(k8,k16,x816); // a(k3+p,k3)

F[84] = - 1./(3.*a)*(alph[64]*G[46]*G[17] // a(k1,k3+p) F2(k3,p) G1(k1)
                  + alph[65]*G[26]*G[21] // a(k3,k1+p) F2(k1,p) G1(k3)
                  + alph[51]*G[58]*G[1] // a(p,-k2) F2(k1,k3) G1(p)

         + alph[66]*G[47]*G[16]
         + alph[67]*G[27]*G[20]
         + alph[50]*G[59]*G[0]

         +3.*G[85]) ;

F[85] =1./(3.*a)*(-3.*(2.-hub)*G[85]-3.*hub*G[84]*muak[3]

        -2.*beta[32]*G[47]*G[17]

        -2.*beta[33]*G[27]*G[21]

        -2.*beta[26]*G[59]*G[1]

			 -2*gam2*gamk[29]*G[16]*G[46]
       -2*gam2*gamk[27]*G[20]*G[26]
			 -2*gam2*gamk[6]*G[0]*G[58]

			 -gamk3[12]*(3.*gam3a + gam3b*(1./pifs[6] + 1./pifs[9] + 1./pifs[10]))*G[0]*G[20]*G[16]);


// F4/G4(-k6,-k8,k1,-k1) = (-k1,-k3,p,-p)


alph[68] = alphai(k46,k2,x218);
alph[69] = alphai(k2,k46,x218);
beta[35] = betai(k46,k2,x218);


alph[70] = alphai(k1,k3,x3);
alph[71] = alphai(k3,k1,x3);
beta[34] = betai(k3,k1,x3);


F[86] = -1./(24.*a)*(6.*alph[23]*G[82]*G[17] // a(k1,k3) F3(-k3,p,-p) G1(k1)
                    +6.*alph[25]*G[74]*G[21] // a(k3,k1) F3(-k1,p,-p) G1(k3)
                    +6.*alph[7]*G[84]*G[1] // a(p,k2-p) F3(-p,-k1,-k3) G1(p)
                    +6.*alph[70]*G[80]*G[1] // a(-p,k2+p) F3(p,-k1,-k3) G1(p)

                   +6.*alph[25]*G[83]*G[16]
                   +6.*alph[23]*G[75]*G[20]
                   +6.*alph[6]*G[85]*G[0]
                   +6.*alph[71]*G[81]*G[0]



                   + 4.*alph[60]*G[58]*G[73] // a(0,k3) F2(k1,k3) G2(p,-p)
                   + 4.*alph[61]*G[72]*G[59] // a(k3,0) F2(p,-p) G2(k1,k3)

                   + 4.*alph[24]*G[62]*G[27] // a(k1+p,k3-p) F2(-k3,p) G2(-p,-k1)
                   + 4.*alph[22]*G[26]*G[63] // a(k3-p,k1+p) F2(-k1,-p) G2(-k3,p)

                   + 4.*alph[69]*G[46]*G[67] // a(k1-p,k3+p) F2(-k3,-p)G2(-k1,p)
                   + 4.*alph[68]*G[66]*G[47] // a(k3+p,k1-p) F2(-k1,p) G2(-k3,-p)

                   +24.*G[87]);


F[87] = 1./a*(-(2.-hub)*G[87]-hub*G[86]*muak[6]

            -1./24.*(
                  12.*beta[12]*G[83]*G[17] // b(k1,k3) G3(-k3,p,-p) G1(k1)
                 +12.*beta[12]*G[75]*G[21] // b(k3,k1) G3(-k1,p,-p) G1(k3)
                 +12.*beta[3]*G[85]*G[1] // b(p,k2-p) G3(-p,-k1,-k3) G1(p)
                 +12.*beta[34]*G[81]*G[1]  // b(-p,k2+p) G3(p,-k1,-k3) G1(p)

                 + 8.*beta[31]*G[59]*G[73] // b(k2,0) G2(k1,k3) G2(p,-p)
                 + 8.*beta[11]*G[63]*G[27]  // b(k3-p,k1+p) G2(-k3,p) G2(k1,p)
                 + 8.*beta[35]*G[67]*G[47] // b(k3+p,k1-p) G2(k3,p) G2(-k1,p)

                 + 8.*gam2*gamk[34]*G[58]*G[72]
                 + 8.*gam2*gamk[13]*G[62]*G[26]
                 + 8.*gam2*gamk[25]*G[66]*G[46]

                 + 12.*gam2*gamk[14]*G[82]*G[16]
                 + 12.*gam2*gamk[14]*G[74]*G[20]
                 + 12.*gam2*gamk[3]*G[84]*G[0]
                 + 12.*gam2*gamk[2]*G[80]*G[0]


	 							+4.*(gamk3[1]*(3.*gam3a + gam3b*(1./pifs[9] + 1./pifs[7] + 1./pifs[2]))*G[16]*G[0]*G[62]
	 									   + gamk3[5]*(3.*gam3a + gam3b*(1./pifs[10] + 1./pifs[5] + 1./pifs[2]))*G[20]*G[0]*G[26]
	 										 + gamk3[17]*(3.*gam3a + gam3b*(1./pifs[4] + 1./pifs[5] + 1./pifs[3]))*G[20]*G[0]*G[66]
	 										 + gamk3[18]*(3.*gam3a + gam3b*(1./pifs[1] + 1./pifs[7] + 1./pifs[3]))*G[16]*G[0]*G[46]
											 + gamk3[23]*(3.*gam3a + gam3b*(1./term1 + 1./pifs[2] + 1./pifs[3]))*G[0]*G[0]*G[58]
									     + gamk3[24]*(3.*gam3a + gam3b*(1./pifs[5] + 1./pifs[7] + 1./pifs[6]))*G[16]*G[20]*G[72])


							 + gamk4[1]*(24.*gam4a+ 8.*gam4b*(1./(pifs[6]*term1) + 1./(pifs[1]*pifs[10]) + 1./(pifs[9]*pifs[4]))
 																		 + 4.*gam4c*(1./pifs[6] + 1./term1 + 1./pifs[1] + 1./pifs[10] + 1./pifs[9] + 1./pifs[4])
																		 + 4./3.*gam4d*(1./(pifs[3]*pifs[6]) + 1./(pifs[3]*pifs[10]) + 1./(pifs[3]*pifs[9])
 																			 			   + 1./(pifs[2]*pifs[6]) + 1./(pifs[2]*pifs[1]) + 1./(pifs[2]*pifs[4])
 																							 + 1./(pifs[5]*term1) + 1./(pifs[5]*pifs[1]) +	1./(pifs[5]*pifs[9])
 																							 + 1./(pifs[7]*term1) + 1./(pifs[7]*pifs[10]) + 1./(pifs[7]*pifs[4]))
 																		 + 12./3.*gam4c*(1./pifs[2] + 1./pifs[3] + 1./pifs[5]  + 1./pifs[7]))*G[16]*G[20]*G[0]*G[0]));





// F3/G3(k6,k8,k1) = (-k2,-k3,-p)

alph[72] = alphai(k7,k46,-x746); // a(k2,k3+p)
alph[73] = alphai(k46,k7,-x746);

alph[74] = alphai(k8,k3,x38); // a(k3,k2+p)
alph[75] = alphai(k3,k8,x38);

beta[36] = betai(k7,k46,-x746);
beta[37] = betai(k8,k3,x38);


F[88] = - 1./(3.*a)*(alph[72]*G[46]*G[19] // a(k2,k3+p) F2(k3,p) G1(k2)
                 + alph[74]*G[36]*G[21] // a(k3,k2+p) F2(k2,p) G1(k3)
                 + alph[53]*G[60]*G[1] // a(p,-k1) F2(k2,k3) G1(p)


        + alph[73]*G[47]*G[18]
        + alph[75]*G[37]*G[20]
        + alph[52]*G[61]*G[0]

        +3.*G[89]) ;


F[89] =1./(3.*a)*(-3.*(2.-hub)*G[89]-3.*hub*G[88]*muak[1]

       -2.*beta[36]*G[47]*G[19]

       -2.*beta[37]*G[37]*G[21]

       -2.*beta[27]*G[61]*G[1]


      -2.*gam2*gamk[21]*G[18]*G[46]
      -2.*gam2*gamk[28]*G[20]*G[36]
      -2.*gam2*gamk[5]*G[0]*G[60]

			-gamk3[13]*(3.*gam3a + gam3b*(1./pifs[5] + 1./pifs[2] + 1./pifs[10]))*G[0]*G[20]*G[18]);


// F3/G3(k6,k8,k1) = (-k2,-k3,p)

 F[90] = - 1./(3.*a)*(alph[29]*G[62]*G[19] // a(k2,k3-p) F2(k3,-p) G1(k2)
                  + alph[40]*G[64]*G[21] // a(k3,k2-p) F2(k2,-p) G1(k3)
                  + alph[11]*G[60]*G[1] // a(p,k1) F2(k2,k3) G1(p)

         + alph[28]*G[63]*G[18]
         + alph[41]*G[65]*G[20]
         + alph[10]*G[61]*G[0]

         +3.*G[91]) ;

 F[91] =1./(3.*a)*(-3.*(2.-hub)*G[91]-3.*hub*G[90]*muak[9]

        -2.*beta[14]*G[63]*G[19]

        -2.*beta[20]*G[65]*G[21]

        -2.*beta[5]*G[61]*G[1]

        -2.*gam2*gamk[16]*G[18]*G[62]
        -2.*gam2*gamk[22]*G[20]*G[64]
        -2.*gam2*gamk[5]*G[0]*G[60]

				-gamk3[14]*(3.*gam3a + gam3b*(1./pifs[5] + 1./pifs[3] + 1./pifs[4]))*G[0]*G[20]*G[18]);


// F4/G4(-k6,-k8,k1,-k1) = (-k2,-k3,p,-p)

F[92] = -1./(24.*a)*(6.*alph[32]*G[82]*G[19] // a(k2,k3) F3(-k3,p,-p) G1(k2)
                   +6.*alph[35]*G[76]*G[21] // a(k3,k2) F3(-k2,p,-p) G1(k3)
                   +6.*alph[5]*G[88]*G[1] // a(p,k1-p) F3(-p,-k2,-k3) G1(p)
                   +6.*alph[31]*G[90]*G[1] // a(-p,k1+p) F3(p,-k2,-k3) G1(p)

                  +6.*alph[35]*G[83]*G[18]
                  +6.*alph[32]*G[77]*G[20]
                  +6.*alph[4]*G[89]*G[0]
                  +6.*alph[34]*G[91]*G[0]


                  + 4.*alph[58]*G[60]*G[73] // a(0,k1) F2(k2,k3) G2(p,-p)
                  + 4.*alph[59]*G[72]*G[61] // a(k1,0) F2(p,-p) G2(k2,k3)
                  + 4.*alph[33]*G[62]*G[37] // a(k2+p,k3-p) F2(-k3,p) G2(p,k2)
                  + 4.*alph[30]*G[36]*G[63] // a(k3-p,k2+p) F2(k2,p) G2(-k3,p)

                  + 4.*alph[42]*G[46]*G[65] // a(k2-p,k3+p) F2(-k3,-p)G2(-k2,p)
                  + 4.*alph[43]*G[64]*G[47] // a(k3+p,k2-p) F2(-k2,p) G2(-k3,-p)


                  +24.*G[93]);

F[93] = 1./a*(-(2.-hub)*G[93]-hub*G[92]*muak[5]

        -1./24.*(
                 12.*beta[17]*G[83]*G[19] // b(k2,k3) G3(-k3,p,-p) G1(k2)
                +12.*beta[17]*G[77]*G[21] // b(k3,k2) G3(-k2,p,-p) G1(k3)
                +12.*beta[2]*G[89]*G[1]// b(p,k1-p) G3(-p,-k2,-k3) G1(p)
                +12.*beta[16]*G[91]*G[1]  // b(-p,k1+p) G3(p,-k2,-k3) G1(p)

                + 8.*beta[30]*G[61]*G[73] // b(k1,0) G2(k2,k3) G2(p,-p)
                + 8.*beta[15]*G[63]*G[37]  // b(k3-p,k2+p) G2(-k3,p) G2(k2,p)
                + 8.*beta[21]*G[65]*G[47] // b(k3+p,k2-p) G2(k3,p) G2(-k2,p)

                + 8.*gam2*gamk[33]*G[60]*G[72]
                + 8.*gam2*gamk[17]*G[62]*G[36]
                + 8.*gam2*gamk[23]*G[64]*G[46]

                + 12.*gam2*gamk[18]*G[82]*G[18]
                + 12.*gam2*gamk[18]*G[76]*G[20]
                + 12.*gam2*gamk[1]*G[88]*G[0]
                + 12.*gam2*gamk[19]*G[90]*G[0]


						+4.*(gamk3[4]*(3.*gam3a + gam3b*(1./pifs[10] + 1./pifs[6] + 1./pifs[2]))*G[18]*G[0]*G[62]
							 + gamk3[2]*(3.*gam3a + gam3b*(1./pifs[2] + 1./pifs[7] + 1./pifs[2]))*G[20]*G[0]*G[36]
							 + gamk3[19]*(3.*gam3a + gam3b*(1./pifs[3] + 1./pifs[7] + 1./pifs[1]))*G[20]*G[0]*G[64]
							 + gamk3[20]*(3.*gam3a + gam3b*(1./pifs[4] + 1./pifs[6] + 1./pifs[1]))*G[18]*G[0]*G[46]
							 + gamk3[25]*(3.*gam3a + gam3b*(1./term1 + 1./pifs[9] + 1./pifs[1]))*G[0]*G[0]*G[60]
							 + gamk3[26]*(3.*gam3a + gam3b*(1./pifs[5] + 1./pifs[7] + 1./pifs[6]))*G[18]*G[20]*G[72])


			 + gamk4[2]*(24.*gam4a+ 8.*gam4b*(1./(pifs[5]*term1) + 1./(pifs[3]*pifs[10]) + 1./(pifs[2]*pifs[4]))
														 + 4.*gam4c*(1./pifs[5] + 1./term1 + 1./pifs[3] + 1./pifs[10] + 1./pifs[2] + 1./pifs[4])
														 + 4./3.*gam4d*(1./(pifs[1]*pifs[5]) + 1./(pifs[1]*pifs[10]) + 1./(pifs[1]*pifs[2])
																			 + 1./(pifs[9]*pifs[5]) + 1./(pifs[9]*pifs[3]) + 1./(pifs[9]*pifs[4])
																			 + 1./(pifs[6]*term1) + 1./(pifs[6]*pifs[3]) +	1./(pifs[6]*pifs[2])
																			 + 1./(pifs[7]*term1) + 1./(pifs[7]*pifs[10]) + 1./(pifs[7]*pifs[4]))
														 + 12./3.*gam4c*(1./pifs[9] + 1./pifs[1] + 1./pifs[6]  + 1./pifs[7]))*G[18]*G[20]*G[0]*G[0]));


  	return GSL_SUCCESS;
  }






	void BSPTN::initnb1_fr(double A, double k[], double x[], double kargs[], double omega0, double par1, double par2, double par3, double epars[])
	{
					double a = 0.0001;

	       double G[94] = { a,-a,a,-a,a,-a, 0.,0.,0.,0.,0.,0.,a,-a,a,-a,a,-a,a,-a,a,-a,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	                            0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

				/*Parameters passed to system of equations */
	          struct param_type6 my_params1 = {k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7],x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],
	                                           kargs[0],kargs[1],kargs[2],kargs[3],kargs[4],kargs[5],kargs[6],kargs[7],kargs[8],
	                                           kargs[9],kargs[10],kargs[11],kargs[12],kargs[13],kargs[14],kargs[15],kargs[16],kargs[17],
	                                           kargs[18],kargs[19],kargs[20], kargs[21], kargs[22],kargs[23],kargs[24],
																						 kargs[25],kargs[26],kargs[27],kargs[28],kargs[29],kargs[30],kargs[31],kargs[32],kargs[33],
																						 kargs[34],kargs[35], omega0, par1, par2, par3};

	          gsl_odeiv2_system sys = {funcfr, jacb, 94, &my_params1};

				//  this gives the system, the method, the initial interval step, the absolute error and the relative error.
					gsl_odeiv2_driver * d =
					gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
														epars[0], epars[1] ,epars[2]);
	// smallest possible accuracy and initial step when comparing to analytic result
	// must reduce initial step size for earlier times!
	// rk8pd is optimal

					int status1 = gsl_odeiv2_driver_apply (d, &a, A, G);


					/*Allocation of array values */


	// For tree level initialisation

	  /*1st order */

				//F1(k;a), G1(k;a)
				F1b_k1 = G[0] ;
				G1b_k1 = G[1] ;

				// F1(k3;a), G1(k3;a)
				F1b_k3 = G[2];
				G1b_k3 = G[3];

				//F1(k2;a), G1(k2;a)
				F1b_k2 = G[4];
				G1b_k2 = G[5];


	// for B321 and B222 initialisaiton

	      /* 1st order */

	      //F1(k1;a), G1(k1;a)
				F1b_1 = G[16] ;
				G1b_1 = G[17] ;

				// F1(k2;a), G1(k2;a)
				F1b_2 = G[18];
				G1b_2 = G[19];

				//F1(k3;a), G1(k3;a)
				F1b_3 = G[20];
				G1b_3 = G[21];


				/*2nd order*/

				//F2/G2(k3,k2) for tree,  F2/G2(k1-p,k2+p) for B222
				F2b_k23 =  G[6];
				G2b_k23 =  G[7];

	      //F2/G2(k1,k3) for tree, F2/G2(-p,k2+p) for B222
	      F2b_k13 =  G[8];
	      G2b_k13 =  G[9];

	      //F2/G2(k1,k2) for tree, F2/G2(p,k1-p)
	      F2b_k12 =  G[10];
	      G2b_k12 =  G[11];


	      // F2/G2(p,k2-p)
	      F2b_p2mp = G[22];
	      G2b_p2mp = G[23];

	      // F2/G2(p,k3-p)
	      F2b_p3mp = G[24];
	      G2b_p3mp = G[25];

	      // F3/G3(k1,p,k2-p)
	      F3b_12mp = G[30];
	      G3b_12mp = G[31];

	      // F3/G3(k1,p,k3-p)
	      F3b_13mp = G[34];
	      G3b_13mp = G[35];

	      // F3/G3(k2,p,k3-p)
	      F3b_23mp = G[40];
	      G3b_23mp = G[41];

	      // F3/G3(k2,p,k1-p)
	      F3b_21mp = G[44];
	      G3b_21mp = G[45];


	      // F3/G3(k3,p,k2-p)
	      F3b_32mp = G[50];
	      G3b_32mp = G[51];

	      // F3/G3(k3,p,k1-p)
	      F3b_31mp = G[54];
	      G3b_31mp = G[55];

	      // F3/G3(k1,k2) for B321-II
	      F2b_12a = G[56];
	      G2b_12a = G[57];

	      // F3/G3(k1,k3) for B321-II
	      F2b_13a = G[58];
	      G2b_13a = G[59];

	      // F3/G3(k2,k3) for B321-II
	      F2b_23a = G[60];
	      G2b_23a = G[61];

	      //F3/G3(k1,p,-p)
	      F3b_1pp = G[74];
	      G3b_1pp = G[75];

	      //F3/G3(k2,p,-p)
	      F3b_2pp = G[76];
	      G3b_2pp = G[77];

	      //F3/G3(k3,p,-p)
	      F3b_3pp = G[82];
	      G3b_3pp = G[83];

	      // F4/G4(-k1,-k2,p,-p)

	      F4b_12pp = G[78];
	      G4b_12pp = G[79];

	      // F4/G4(-k1,-k3,p,-p)

	      F4b_13pp = G[86];
	      G4b_13pp = G[87];

	      // F4/G4(-k2,-k3,p,-p)

	      F4b_23pp = G[92];
	      G4b_23pp = G[93];

				gsl_odeiv2_driver_free(d);
	}
