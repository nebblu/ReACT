#ifndef SPT_H
#define SPT_H

#include "Common.h"

/******************************************************************************
 * SPT
 *
 * 1-loop standard perturbation theory.  Indices (a,b) are used to indicate
 * either the density-density, density-velocity, or velocity-velocity spectra,
 *   P_{11} <--> P_{\delta\delta}
 *   P_{12} <--> P_{\delta\theta}
 *   P_{22} <--> P_{\theta\theta}
 * whereas the labels "22" or "13" refer to different terms in the perturbation
 * series.
 ******************************************************************************/

extern double phpars[13];

class SPT {
public:
    SPT(const Cosmology& C, const PowerSpectrum& P_L, real epsrel = 1e-3);


	/* P_{ab}(k) 1-loop: analytical LCDM and nDGP */
	/* The division by Dl is my own normalisation of P_L - it requires z=0 always in spt2.cpp */
  // values of a :
  // 1: P_dd, 2: P_dt , 3: P_tt ; LCDM
	// 4:  P_dd, 5: P_dt , 6: P_tt ; nDGP
    real PLOOP(real k, int a) const;

	/* P_{ab}(k) 1-loop : numerical for arbitrary model of gravity - kernel dependent */
  // values of a :
  // 1: P_dd, 2: P_dt , 3: P_tt
  // ymin = ymax = k for scale dependent models
  // ymin = QMINp/kmax , ymax = QMAXp/kmin for scale independant models
  // where QMINp = 0.0001, QMAXp = 5. and kmin and kmax are your specified output ranges in examples/spt.cpp
	real PLOOPn(double ymin, double ymax, int a,  real k ) const;
  // optimised PLOOPn
  real PLOOPn2(int a, double vars[], double k, double err) const;
  // halo fit coefficient initialisation
  void phinit(double scalef, double omega0) const;
  // halofit model
  double PHALO(double k) const;

	/*Redshift Space PS */
  // bl is galaxy bias
	/* KAISER RSD */
	// a =1 : LCDM
	// a =2 : nDGP
	// a =3 : numerical (arbitrary model)
    real PRSD(real k,real u, real bl, int a) const;


	/* P(k,u) RSD TNS */
  // sigma_v is the damping factor
  // u is the cosine of the angle between k and LOS
	// a = 1 : LCDM
	// a = 2 : nDGP
    real PTNS(real k, real u, real bl, real sigma_v, int a) const;

	/* P(k) RSD TNS: Arbitrary model   */
  // ymin = ymax = k for scale dependent models
  // ymin = QMINp/kmax , ymax = QMAXp/kmin for scale independant models
  // where QMINp = 0.0001, QMAXp = 5. and kmin and kmax are your specified output ranges in examples/spt.cpp
	real PTNSn(double ymin, double ymax, real bl, real sigma_v, real k, real u) const ;

	/*Multipoles*/

	/*  LCDM  Kaiser Multipoles with DFoG term (set sigma_v=0 for linear kaiser multipoles) */

	// a = 1 : Monopole
	// a = 2 : Quadrupole
	// a = 3 : Hexdecapole

	real KASM(real k, real bl, real sigma_v, int a) const ;


	/*  Arbitrary model Kaiser Multipoles with DFoG term (set sigma_v=0 for linear kaiser multipoles) */
	// a = 1 : Monopole
	// a = 2 : Quadrupole
	// a = 3 : Hexdecapole
	real KASMmg(real k, real bl, real sigma_v,   int a) const ;



	/*  LCDM  TNS  Multipoles */
	// a = 1 : Monopole
	// a = 2 : Quadrupole
	// a = 3 : Hexdecapole
  //linear bias
	real PTNSM(real k, real bl, real sigma_v, int a) const;
  // qbias
  real PTNSMq(real k, double barr[], real sigma_v, int a) const;

  real PTNSMl(real k, double barr[], real sigma_v, int a) const;

  real Lag_bias(int a, real k, real bias[]) const;



	/*  nDGP TNS  Multipoles  */
	// a = 1 : Monopole
	// a = 2 : Quadrupole
	// a = 3 : Hexdecapole
  //linear bias
	real PTNSMnDGP(real k, real bl, real sigma_v, int a) const ;
  // qbias
  real PTNSMnDGPq(real k, double barr[], real sigma_v, int a) const;


	/*  Arbitrary model TNS  Multipoles with DFoG term  */
	// a = 1 : Monopole
	// a = 2 : Quadrupole
	// a = 3 : Hexdecapole
  // ymin = ymax = k for scale dependent models
  // ymin = QMINp/kmax , ymax = QMAXp/kmin for scale independant models
  // where QMINp = 0.0001, QMAXp = 5. and kmin and kmax are your specified output ranges in examples/spt.cpp
	real PTNSMmg(real k, double ymin, double ymax, int a, real bl, real sigma_v) const ;

	//LCDM 1-loop terms

    /* $P_{ab}^{(2,2)}(k)$ */
    real P22(real k, int a=1, int b=1) const;

    /* $P_{ab}^{(1,3)}(k)$ */
    real P13(real k, int a=1, int b =1) const;

    /* $P_{\delta\delta}^{(2,2)}(k)$ */
    real P22_dd(real k) const;

    /* $P_{\delta\theta}^{(2,2)}(k)$ */
    real P22_dt(real k) const;

    /* $P_{\theta\theta}^{(2,2)}(k)$ */
	real P22_tt(real k) const;

    /* $P_{\delta\delta}^{(1,3)}(k)$ */
    real P13_dd(real k) const;

    /* $P_{\delta\theta}^{(1,3)}(k)$ */
    real P13_dt(real k) const;

    /* $P_{\theta\theta}^{(1,3)}(k)$ */
	real P13_tt(real k) const ;


	//DGP 1-loop terms

 /* $P_{ab}^{(2,2)}(k)$ */
    real P22D(real k, int a) const;

    /* $P_{ab}^{(1,3)}(k)$ */
    real P13D(real k, int a) const;

    /* $P_{\delta\delta}^{(2,2)}(k)$ */
    real P22D_dd(real k) const;

    /* $P_{\delta\theta}^{(2,2)}(k)$ */
    real P22D_dt(real k) const;

    /* $P_{\theta\theta}^{(2,2)}(k)$ */
     real P22D_tt(real k) const;

    /* $P_{\delta\delta}^{(1,3)}(k)$ */
    real P13D_dd(real k) const;

    /* $P_{\delta\theta}^{(1,3)}(k)$ */
    real P13D_dt(real k) const;

    /* $P_{\theta\theta}^{(1,3)}(k)$ */
     real P13D_tt(real k) const ;

	//Numerical 1-loop terms
	real P22n(double ymin, double ymax, int a, real k) const;

	real P13n(double ymin, double ymax, int a, real k) const;


  real PRESUM(real k) const;



    /* Linear growth fixing */
    void G(real s8) const;

    /* Log growth free param */

    void remp(real f) const;


    void rempn(real f) const;

    void Gn(real dfid) const ;


	/* A+B+C term analytical*/

  // a = 0:  RSD PS
  // a = 1 : Monopole
  // a = 2 : Quadrupole
  // a = 3 : Hexdecapole

  //U = u for 2d PS
  //U = sigma_v for multipoles
	real AB(real k, real bl, real U, int a) const;


    /* A+B+C term numerical */
    // e = 0:  RSD PS
    // e = 1 : Monopole
    // e = 2 : Quadrupole
    // e = 3 : Hexdecapole

    //U = u for 2D PS
    //U = sigma_v for multipoles
	real ABn(double ymin, double ymax, real bl, real k, real U, int e) const;

// used for scale dep interpolation scheme
  void ABs_selec(double myarray[], real k) const ;


private:
    const Cosmology& C;
    const PowerSpectrum& P_L;
    real epsrel;

};

#endif // SPT_H
