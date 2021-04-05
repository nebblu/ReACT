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


/* REAL SPACE */

	/* P_{ab}(k) 1-loop: analytical (EdS) LCDM and nDGP */
	/* The division by Dl is my own normalisation of P_L - it requires z=0 always in spt2.cpp */
  // values of a :
  // 1: P_dd, 2: P_dt , 3: P_tt ; LCDM
	// 4:  P_dd, 5: P_dt , 6: P_tt ; nDGP
  real PLOOP(real k, int a) const;

	/* P_{ab}(k) 1-loop : numerical for arbitrary model of gravity for massless/no neutrinos- see  1606.02520  */
  // values of a :
  // 0: P_dd linear
  // 1: P_dd, 2: P_dt , 3: P_tt for mg @ 1-loop level
  // 4: P_dd pseudo @ 1-loop level - used for halo model reactions (see HALO.cpp)
  // vars: 0:scale factor , 1: omega_0, 2: mg1, 3 : mg2 , 4 : mg3
  real PLOOPn2(int a, double vars[], int model, double k, double err) const;

  /* P_{ab}(k) 1-loop : numerical for arbitrary model of gravity with massive neutrinos - kernel dependent */
  // 0: P_dd for mg linear
  // 1: P_dd for mg @ 1-loop
  // 2: P_dd pseudo @ 1-loop level - used for halo model reactions (see HALO.cpp)
  real PLOOPn2_nu(int a, double vars[], int model, double k, double err) const;

  // initialise p_loop values over redshifts[] at k=k0 for k_star in reaction code (HALO.cpp)
  void ploop_init(double ploopr[], double ploopp[], double redshifts[], int noz, double vars[], int model, double k0);

// Work in progress -- initialisation of 1-loop k0 predictions for multiple redshifts
  //void ploop_init_nu(double ploopr[], double ploopp[], double redshifts[], int noz, double vars[], double k0);

  // halo fit coefficient initialisation
  void phinit(double scalef, double omega0) const;
  // halofit model
  double PHALO(double k) const;

/*Redshift Space PS */

    /* Funtions to vary growth rate freely */

    void remp(real f) const;
    void rempn(real f) const;
    void rempdgp(real f) const;

  // bl is galaxy bias
  // u is the cosine of the angle between k and the line of sight

	/* KAISER RSD */
	// a =1 : LCDM
	// a =2 : nDGP
  real PRSD(real k,real u, real bl, int a) const;


	/* P(k,u) RSD TNS */
  // sigma_v is the damping factor
  // u is the cosine of the angle between k and LOS
	// a = 1 : LCDM
	// a = 2 : nDGP
  real PTNS(real k, real u, real bl, real sigma_v, int a) const;

	/*Multipoles*/

	/*  LCDM  Kaiser Multipoles with DFoG term (set sigma_v=0 for linear kaiser multipoles) */

	// a = 1 : Monopole
	// a = 2 : Quadrupole
	// a = 3 : Hexdecapole
	real KASM(real k, real bl, real sigma_v, int a) const ;

	/*  LCDM  TNS  Multipoles */
	// a = 1 : Monopole
	// a = 2 : Quadrupole
	// a = 3 : Hexdecapole
  //linear bias
	real PTNSM(real k, real bl, real sigma_v, int a) const;

  // qbias
  real PTNSMq(real k, double barr[], real sigma_v, int a) const;

  // lagrangian bias
  real PTNSMl(real k, double barr[], real sigma_v, int a) const;

  	/*  nDGP TNS  Multipoles  */
  	// a = 1 : Monopole
  	// a = 2 : Quadrupole
  	// a = 3 : Hexdecapole
    //linear bias
  	real PTNSMnDGP(real k, real bl, real sigma_v, int a) const ;
    // qbias
    real PTNSMnDGPq(real k, double barr[], real sigma_v, int a) const;


  // modified gravity TNS and Kaiser model - numerically calculated
  double PRSD_mg(int a, int b, double bias[], double vars[], int model, double sigma_v, double k, double err) const;

/// lagrangian bias terms
  real Lag_bias(int a, real k, real bias[]) const;


// individual loop terms
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


	/* A+B+C term analytical*/
  // a = 1 : Monopole
  // a = 2 : Quadrupole
  // a = 3 : Hexdecapole
	real AB(real k, real bl, real sigmav, int a) const;
  // for 2d spectra
  real AB_mu(real k, real bl, real u, int a) const;


private:
    const Cosmology& C;
    const PowerSpectrum& P_L;
    real epsrel;

};

#endif // SPT_H
