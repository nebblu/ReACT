#if HAVE_CONFIG_H
# include <config.h>
#endif

#include "PowerSpectrum.h"
#include "Quadrature.h"
#include "SpecialFunctions.h"
#include "array.h"
#include "SPT.h"
#include "RegPT.h"
#include "PMG_Interpolation.h"
#include <gsl/gsl_sf_bessel.h>


#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

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


// Initialization of 17 2dsplines over mg param and k for :
// P_dd, P_dt, P_tt, A11, A12, A22, A21, A33, B12,B13, B14, B22, B23, B24, B33, B34, B44
const gsl_interp2d_type *PDDs = gsl_interp2d_bilinear;
const gsl_interp2d_type *PDTs = gsl_interp2d_bilinear;
const gsl_interp2d_type *PTTs = gsl_interp2d_bilinear;
const gsl_interp2d_type *A11s = gsl_interp2d_bilinear;
const gsl_interp2d_type *A12s = gsl_interp2d_bilinear;
const gsl_interp2d_type *A22s = gsl_interp2d_bilinear;
const gsl_interp2d_type *A21s = gsl_interp2d_bilinear;
const gsl_interp2d_type *A33s = gsl_interp2d_bilinear;
const gsl_interp2d_type *B12s = gsl_interp2d_bilinear;
const gsl_interp2d_type *B13s = gsl_interp2d_bilinear;
const gsl_interp2d_type *B14s = gsl_interp2d_bilinear;
const gsl_interp2d_type *B22s = gsl_interp2d_bilinear;
const gsl_interp2d_type *B23s = gsl_interp2d_bilinear;
const gsl_interp2d_type *B24s = gsl_interp2d_bilinear;
const gsl_interp2d_type *B33s = gsl_interp2d_bilinear;
const gsl_interp2d_type *B34s = gsl_interp2d_bilinear;
const gsl_interp2d_type *B44s = gsl_interp2d_bilinear;


// create new accels
gsl_interp_accel *kacc1 = gsl_interp_accel_alloc();
gsl_interp_accel *kacc2 = gsl_interp_accel_alloc();
gsl_interp_accel *kacc3 = gsl_interp_accel_alloc();
gsl_interp_accel *kacc4 = gsl_interp_accel_alloc();
gsl_interp_accel *kacc5 = gsl_interp_accel_alloc();
gsl_interp_accel *kacc6 = gsl_interp_accel_alloc();
gsl_interp_accel *kacc7 = gsl_interp_accel_alloc();
gsl_interp_accel *kacc8 = gsl_interp_accel_alloc();
gsl_interp_accel *kacc9 = gsl_interp_accel_alloc();
gsl_interp_accel *kacc10 = gsl_interp_accel_alloc();
gsl_interp_accel *kacc11 = gsl_interp_accel_alloc();
gsl_interp_accel *kacc12 = gsl_interp_accel_alloc();
gsl_interp_accel *kacc13 = gsl_interp_accel_alloc();
gsl_interp_accel *kacc14 = gsl_interp_accel_alloc();
gsl_interp_accel *kacc15 = gsl_interp_accel_alloc();
gsl_interp_accel *kacc16 = gsl_interp_accel_alloc();
gsl_interp_accel *kacc17 = gsl_interp_accel_alloc();


gsl_interp_accel *mgacc1 = gsl_interp_accel_alloc();
gsl_interp_accel *mgacc2 = gsl_interp_accel_alloc();
gsl_interp_accel *mgacc3 = gsl_interp_accel_alloc();
gsl_interp_accel *mgacc4 = gsl_interp_accel_alloc();
gsl_interp_accel *mgacc5 = gsl_interp_accel_alloc();
gsl_interp_accel *mgacc6 = gsl_interp_accel_alloc();
gsl_interp_accel *mgacc7 = gsl_interp_accel_alloc();
gsl_interp_accel *mgacc8 = gsl_interp_accel_alloc();
gsl_interp_accel *mgacc9 = gsl_interp_accel_alloc();
gsl_interp_accel *mgacc10 = gsl_interp_accel_alloc();
gsl_interp_accel *mgacc11 = gsl_interp_accel_alloc();
gsl_interp_accel *mgacc12 = gsl_interp_accel_alloc();
gsl_interp_accel *mgacc13 = gsl_interp_accel_alloc();
gsl_interp_accel *mgacc14 = gsl_interp_accel_alloc();
gsl_interp_accel *mgacc15 = gsl_interp_accel_alloc();
gsl_interp_accel *mgacc16 = gsl_interp_accel_alloc();
gsl_interp_accel *mgacc17 = gsl_interp_accel_alloc();




// k values to interpolate between: k in [0.001,0.25]
const int kvint = 50;
double kvals[kvint] = {1.000000e-03 , 1.119277e-03 , 1.252781e-03 , 1.402209e-03 , 1.569460e-03 , 1.756660e-03 , 1.966189e-03 , 2.200710e-03 , 2.463204e-03 , 2.757008e-03 , 3.085855e-03 , 3.453926e-03 , 3.865900e-03 , 4.327013e-03 , 4.843125e-03 , 5.420799e-03 , 6.067375e-03 , 6.791073e-03 , 7.601091e-03 , 8.507726e-03 , 9.522501e-03 , 1.065832e-02 , 1.192961e-02 , 1.335253e-02 , 1.494518e-02 , 1.672780e-02 , 1.872304e-02 , 2.095627e-02 , 2.345586e-02 , 2.625361e-02 , 2.938506e-02 , 3.289002e-02 , 3.681304e-02 , 4.120398e-02 , 4.611867e-02 , 5.161956e-02 , 5.777658e-02 , 6.466799e-02 , 7.238139e-02 , 8.101482e-02 , 9.067802e-02 , 1.014938e-01 , 1.135997e-01 , 1.271495e-01 , 1.423155e-01 , 1.592905e-01 , 1.782902e-01 , 1.995561e-01 , 2.233585e-01 , 2.500000e-01};

// Allocate memory for PT values to be splined accross
  double *pdd_pmgs;
  double *pdt_pmgs;
  double *ptt_pmgs;
  double *a11_pmgs;
  double *a12_pmgs;
  double *a22_pmgs;
  double *a21_pmgs;
  double *a33_pmgs;
  double *b12_pmgs;
  double *b13_pmgs;
  double *b14_pmgs;
  double *b22_pmgs;
  double *b23_pmgs;
  double *b24_pmgs;
  double *b33_pmgs;
  double *b34_pmgs;
  double *b44_pmgs;

  // Point to an allocated interpolation object
    gsl_interp2d *pdd_pmg ;
    gsl_interp2d *pdt_pmg ;
    gsl_interp2d *ptt_pmg ;
    gsl_interp2d *a11_pmg ;
    gsl_interp2d *a12_pmg ;
    gsl_interp2d *a22_pmg ;
    gsl_interp2d *a21_pmg ;
    gsl_interp2d *a33_pmg ;
    gsl_interp2d *b12_pmg ;
    gsl_interp2d *b13_pmg ;
    gsl_interp2d *b14_pmg ;
    gsl_interp2d *b22_pmg ;
    gsl_interp2d *b23_pmg ;
    gsl_interp2d *b24_pmg ;
    gsl_interp2d *b33_pmg ;
    gsl_interp2d *b34_pmg ;
    gsl_interp2d *b44_pmg ;



  // Load Linear PS
  PMGI::PMGI(const Cosmology& C, const PowerSpectrum& P_l, real epsrel_)
  : C(C), P_l(P_l)
  {
      epsrel = epsrel_;
  }


static void ptvalues_init_spt(double scalef, double omega0, double anarray1[][50][17], double anarray2[] , const int N1, const PowerSpectrum& P_l, const Cosmology& C){
  // Load classes for RegPT and SPT and kernel intializer
  SPT spt(C, P_l, 1e-3);
  IOW iow;
  iow.inite(scalef, omega0 , 1e-14, 1., 1.); // Initialize normalization of linear growth
   warning("Number of MG values = %d\n", N1);
  //  #pragma omp parallel for
  for(int j=0; j<N1; j++){
    for(int i = 0; i<kvint; i++){
    // still need to extend mgvals to 2d array
    double k = kvals[i];
    double abarray[14];
    iow.initn(scalef, QMINp/k ,QMAXp/k, k, omega0, anarray2[j] , 1., 1.);
    // Initialize AB terms
    spt.ABs_selec(abarray, k);
    anarray1[j][i][0] = spt.PLOOPn(k,k,1,k); //pdd
    anarray1[j][i][1] = spt.PLOOPn(k,k,2,k); //pdt
    anarray1[j][i][2] = spt.PLOOPn(k,k,3,k); //ptt
    anarray1[j][i][3] = abarray[0]; //a11
    anarray1[j][i][4] = abarray[1]; //a12
    anarray1[j][i][5] = abarray[2]; //a22
    anarray1[j][i][6] = abarray[3]; // a21
    anarray1[j][i][7] = abarray[4];  //a33
    anarray1[j][i][8] = abarray[5];  //b12
    anarray1[j][i][9] = abarray[6];  //b13
    anarray1[j][i][10] = abarray[7]; //b14
    anarray1[j][i][11] = abarray[8]; //b22
    anarray1[j][i][12] = abarray[9]; //b23
    anarray1[j][i][13] = abarray[10]; //b24
    anarray1[j][i][14] = abarray[11]; //b33
    anarray1[j][i][15] = abarray[12]; //b34
    anarray1[j][i][15] = abarray[13]; //b44
  }
  }
  return;
}

static void ptvalues_init_rpt(double scalef, double omega0, double anarray1[][50][17], double anarray2[], const int N1, const PowerSpectrum& P_l,const Cosmology& C){
  // Load classes for RegPT and SPT and kernel intializer
  RegPT rpt(C, P_l, 1e-3);
  IOW iow;
  iow.inite(scalef, omega0 , 1e-14, 1., 1.); // Initialize normalization of linear growth
//  #pragma omp parallel for
  for(int j=0; j<N1; j++){
    iow.inite(scalef, omega0, anarray2[j],  1., 1.);
    rpt.sigmad_init(); // initialize sigma_d damping term
    for(int i = 0; i<kvint; i++){
    double abarray[14];
    // still need to extend mgvals to 2d array
    double k = kvals[i];
    iow.initn(scalef, QMINp/k ,QMAXp/k, k, omega0, anarray2[j] , 1., 1.);
    // Initialize AB terms
    rpt.ABr_selec(abarray, k);
    anarray1[j][i][0] = rpt.PLOOPnr(k,k,1,k); //pdd
    anarray1[j][i][1] = rpt.PLOOPnr(k,k,2,k); //pdt
    anarray1[j][i][2] = rpt.PLOOPnr(k,k,3,k); //ptt
    anarray1[j][i][3] = abarray[0]; //a11
    anarray1[j][i][4] = abarray[1]; //a12
    anarray1[j][i][5] = abarray[2]; //a22
    anarray1[j][i][6] = abarray[3]; // a21
    anarray1[j][i][7] = abarray[4];  //a33
    anarray1[j][i][8] = abarray[5];  //b12
    anarray1[j][i][9] = abarray[6];  //b13
    anarray1[j][i][10] = abarray[7]; //b14
    anarray1[j][i][11] = abarray[8]; //b22
    anarray1[j][i][12] = abarray[9]; //b23
    anarray1[j][i][13] = abarray[10]; //b24
    anarray1[j][i][14] = abarray[11]; //b33
    anarray1[j][i][15] = abarray[12]; //b34
    anarray1[j][i][15] = abarray[13]; //b44
  }
  }
  return;
}

const int mgv = 20;
double mgvals_glob[mgv];

// a selects SPT or RegPT
void PMGI::spline_init(double scalef, double omega0, double mgvals[], const int mgsize, int a)const{

    for(int i = 1; i<mgsize; i++){
      if (mgvals[i]<=mgvals[i-1]) {
        warning("PMG: mgvals need to be monotonically inceasing");
        return;
      }
      }

      if (mgsize>mgv) {
        warning("PMG: Only first 20 parms will be used");
      }

      for(int i = 0; i<mgv; i++){
        if (i<mgsize) {
          mgvals_glob[i] = mgvals[i];
        }
        else{
        mgvals_glob[i] = mgvals_glob[i-1]+1;
      }
      }

  // Allocate memory for PT values to be splined accross
    pdd_pmgs = (double*)malloc(kvint * mgv * sizeof(double));
    pdt_pmgs = (double*)malloc(kvint * mgv * sizeof(double));
    ptt_pmgs = (double*)malloc(kvint * mgv * sizeof(double));
    a11_pmgs = (double*)malloc(kvint * mgv * sizeof(double));
    a12_pmgs = (double*)malloc(kvint * mgv * sizeof(double));
    a22_pmgs = (double*)malloc(kvint * mgv * sizeof(double));
    a21_pmgs = (double*)malloc(kvint * mgv * sizeof(double));
    a33_pmgs = (double*)malloc(kvint * mgv * sizeof(double));
    b12_pmgs = (double*)malloc(kvint * mgv * sizeof(double));
    b13_pmgs = (double*)malloc(kvint * mgv * sizeof(double));
    b14_pmgs = (double*)malloc(kvint * mgv * sizeof(double));
    b22_pmgs = (double*)malloc(kvint * mgv * sizeof(double));
    b23_pmgs = (double*)malloc(kvint * mgv * sizeof(double));
    b24_pmgs = (double*)malloc(kvint * mgv * sizeof(double));
    b33_pmgs = (double*)malloc(kvint * mgv * sizeof(double));
    b34_pmgs = (double*)malloc(kvint * mgv * sizeof(double));
    b44_pmgs = (double*)malloc(kvint * mgv * sizeof(double));


    // Point to an allocated interpolation object
      pdd_pmg = gsl_interp2d_alloc(PDDs, kvint, mgv);
      pdt_pmg = gsl_interp2d_alloc(PDTs, kvint, mgv);
      ptt_pmg = gsl_interp2d_alloc(PTTs, kvint, mgv);
      a11_pmg = gsl_interp2d_alloc(A11s, kvint, mgv);
      a12_pmg = gsl_interp2d_alloc(A12s, kvint, mgv);
      a22_pmg = gsl_interp2d_alloc(A22s, kvint, mgv);
      a21_pmg = gsl_interp2d_alloc(A21s, kvint, mgv);
      a33_pmg = gsl_interp2d_alloc(A33s, kvint, mgv);
      b12_pmg = gsl_interp2d_alloc(B12s, kvint, mgv);
      b13_pmg = gsl_interp2d_alloc(B13s, kvint, mgv);
      b14_pmg = gsl_interp2d_alloc(B14s, kvint, mgv);
      b22_pmg = gsl_interp2d_alloc(B22s, kvint, mgv);
      b23_pmg = gsl_interp2d_alloc(B23s, kvint, mgv);
      b24_pmg = gsl_interp2d_alloc(B24s, kvint, mgv);
      b33_pmg = gsl_interp2d_alloc(B33s, kvint, mgv);
      b34_pmg = gsl_interp2d_alloc(B34s, kvint, mgv);
      b44_pmg = gsl_interp2d_alloc(B44s, kvint, mgv);


  double my_array[mgv][kvint][17]; // array to hold PT values

// Initialize PT values for the k grid and mg grid
  switch (a) {
    case 1:
    ptvalues_init_spt(scalef,omega0,my_array,mgvals,mgsize,cref(P_l),cref(C));
    break;
    case 2:
    ptvalues_init_rpt(scalef,omega0,my_array,mgvals,mgsize,cref(P_l),cref(C));
    break;
    default:
    warning("PMG: invalid indices, a = %d\n", a);
    return;
  }


  /* set z grid values */

  for(int i = 0; i<kvint ; i ++){
    for(int j = 0; j<mgv; j++){

      gsl_interp2d_set(pdd_pmg, pdd_pmgs, i, j, my_array[j][i][0]);
      gsl_interp2d_set(pdt_pmg, pdt_pmgs, i, j, my_array[j][i][1]);
      gsl_interp2d_set(ptt_pmg, ptt_pmgs, i, j, my_array[j][i][2]);
      gsl_interp2d_set(a11_pmg, a11_pmgs, i, j, my_array[j][i][3]);
      gsl_interp2d_set(a12_pmg, a12_pmgs, i, j, my_array[j][i][4]);
      gsl_interp2d_set(a22_pmg, a22_pmgs, i, j, my_array[j][i][5]);
      gsl_interp2d_set(a21_pmg, a21_pmgs, i, j, my_array[j][i][6]);
      gsl_interp2d_set(a33_pmg, a33_pmgs, i, j, my_array[j][i][7]);
      gsl_interp2d_set(b12_pmg, b12_pmgs, i, j, my_array[j][i][8]);
      gsl_interp2d_set(b13_pmg, b13_pmgs, i, j, my_array[j][i][9]);
      gsl_interp2d_set(b14_pmg, b14_pmgs, i, j, my_array[j][i][10]);
      gsl_interp2d_set(b22_pmg, b22_pmgs, i, j, my_array[j][i][11]);
      gsl_interp2d_set(b23_pmg, b23_pmgs, i, j, my_array[j][i][12]);
      gsl_interp2d_set(b24_pmg, b24_pmgs, i, j, my_array[j][i][13]);
      gsl_interp2d_set(b33_pmg, b33_pmgs, i, j, my_array[j][i][14]);
      gsl_interp2d_set(b34_pmg, b34_pmgs, i, j, my_array[j][i][15]);
      gsl_interp2d_set(b44_pmg, b44_pmgs, i, j, my_array[j][i][16]);
    }
  }


  /* initialize interpolation */
  gsl_interp2d_init(pdd_pmg, kvals, mgvals_glob, pdd_pmgs, kvint, mgv);
  gsl_interp2d_init(pdt_pmg, kvals, mgvals_glob, pdt_pmgs, kvint, mgv);
  gsl_interp2d_init(ptt_pmg, kvals, mgvals_glob, ptt_pmgs, kvint, mgv);
  gsl_interp2d_init(a11_pmg, kvals, mgvals_glob, a11_pmgs, kvint, mgv);
  gsl_interp2d_init(a12_pmg, kvals, mgvals_glob, a12_pmgs, kvint, mgv);
  gsl_interp2d_init(a22_pmg, kvals, mgvals_glob, a22_pmgs, kvint, mgv);
  gsl_interp2d_init(a21_pmg, kvals, mgvals_glob, a21_pmgs, kvint, mgv);
  gsl_interp2d_init(a33_pmg, kvals, mgvals_glob, a33_pmgs, kvint, mgv);
  gsl_interp2d_init(b12_pmg, kvals, mgvals_glob, b12_pmgs, kvint, mgv);
  gsl_interp2d_init(b13_pmg, kvals, mgvals_glob, b13_pmgs, kvint, mgv);
  gsl_interp2d_init(b14_pmg, kvals, mgvals_glob, b14_pmgs, kvint, mgv);
  gsl_interp2d_init(b22_pmg, kvals, mgvals_glob, b22_pmgs, kvint, mgv);
  gsl_interp2d_init(b23_pmg, kvals, mgvals_glob, b23_pmgs, kvint, mgv);
  gsl_interp2d_init(b24_pmg, kvals, mgvals_glob, b24_pmgs, kvint, mgv);
  gsl_interp2d_init(b33_pmg, kvals, mgvals_glob, b33_pmgs, kvint, mgv);
  gsl_interp2d_init(b34_pmg, kvals, mgvals_glob, b34_pmgs, kvint, mgv);
  gsl_interp2d_init(b44_pmg, kvals, mgvals_glob, b44_pmgs, kvint, mgv);

}

// TNS Spectra using interpolated functions
double PMGI::PMG_TNS(double k, int a, double bl, double sigma_v, double mg1)const{
  //calculate components using splined objects
  double pdd =gsl_interp2d_eval(pdd_pmg, kvals, mgvals_glob, pdd_pmgs, k,mg1, kacc1, mgacc1);
  double pdt =gsl_interp2d_eval(pdt_pmg, kvals, mgvals_glob, pdt_pmgs, k,mg1, kacc2, mgacc2);
  double ptt =gsl_interp2d_eval(ptt_pmg, kvals, mgvals_glob, ptt_pmgs ,k,mg1, kacc3, mgacc3);
  double a11 =gsl_interp2d_eval(a11_pmg, kvals, mgvals_glob, a11_pmgs ,k,mg1, kacc4, mgacc4);
  double a12 =gsl_interp2d_eval(a12_pmg, kvals, mgvals_glob, a12_pmgs ,k,mg1, kacc5, mgacc5);
  double a22 =gsl_interp2d_eval(a22_pmg, kvals, mgvals_glob, a22_pmgs ,k,mg1, kacc6, mgacc6);
  double a21 =gsl_interp2d_eval(a21_pmg, kvals, mgvals_glob, a21_pmgs ,k,mg1, kacc7, mgacc7);
  double a33 =gsl_interp2d_eval(a33_pmg, kvals, mgvals_glob, a33_pmgs ,k,mg1, kacc8, mgacc8);
  double b12 =gsl_interp2d_eval(b12_pmg, kvals, mgvals_glob, b12_pmgs ,k,mg1, kacc9, mgacc9);
  double b13 =gsl_interp2d_eval(b13_pmg, kvals, mgvals_glob, b13_pmgs ,k,mg1, kacc10, mgacc10);
  double b14 =gsl_interp2d_eval(b14_pmg, kvals, mgvals_glob, b14_pmgs ,k,mg1, kacc11, mgacc11);
  double b22 =gsl_interp2d_eval(b22_pmg, kvals, mgvals_glob, b22_pmgs ,k,mg1, kacc12, mgacc12);
  double b23 =gsl_interp2d_eval(b23_pmg, kvals, mgvals_glob, b23_pmgs ,k,mg1, kacc13, mgacc13);
  double b24 =gsl_interp2d_eval(b24_pmg, kvals, mgvals_glob, b24_pmgs ,k,mg1, kacc14, mgacc14);
  double b33 =gsl_interp2d_eval(b33_pmg, kvals, mgvals_glob, b33_pmgs ,k,mg1, kacc15, mgacc15);
  double b34 =gsl_interp2d_eval(b34_pmg, kvals, mgvals_glob, b34_pmgs ,k,mg1, kacc16, mgacc16);
  double b44 =gsl_interp2d_eval(b44_pmg, kvals, mgvals_glob, b44_pmgs ,k,mg1, kacc17, mgacc17);

  //calculate multipole factors
 double u00 = factL(k,sigma_v,fdgp_spt,1.,0,a,6);
 double u01 = factL(k,sigma_v,fdgp_spt,1.,1,a,6);
 double u02 = factL(k,sigma_v,fdgp_spt,1.,2,a,6);
 double u03 = factL(k,sigma_v,fdgp_spt,1.,3,a,6);
 double u04 = factL(k,sigma_v,fdgp_spt,1.,4,a,6);

 return  pow2(bl)*u00*pdd - 2*u01*bl*pdt + u02*ptt
         + pow2(bl)*(u01*a11 + u01*b12 + u02*b22)
         + bl*(u01*a12 + u02*a22 + u01*b13 + u02*b23 + u03*b33)
         + u02*a21 + u03*a33 + u01*b14 + u02*b24 + u03*b34 + u04*b44;

}

// free objects
void PMGI::PMG_FREE()const{
  gsl_interp2d_free(pdd_pmg);
  gsl_interp2d_free(pdt_pmg);
  gsl_interp2d_free(ptt_pmg);
  gsl_interp2d_free(a11_pmg);
  gsl_interp2d_free(a12_pmg);
  gsl_interp2d_free(a22_pmg);
  gsl_interp2d_free(a21_pmg);
  gsl_interp2d_free(a33_pmg);
  gsl_interp2d_free(b12_pmg);
  gsl_interp2d_free(b13_pmg);
  gsl_interp2d_free(b14_pmg);
  gsl_interp2d_free(b22_pmg);
  gsl_interp2d_free(b23_pmg);
  gsl_interp2d_free(b24_pmg);
  gsl_interp2d_free(b33_pmg);
  gsl_interp2d_free(b34_pmg);
  gsl_interp2d_free(b44_pmg);
  gsl_interp_accel_free(kacc1);
  gsl_interp_accel_free(kacc2);
  gsl_interp_accel_free(kacc3);
  gsl_interp_accel_free(kacc4);
  gsl_interp_accel_free(kacc5);
  gsl_interp_accel_free(kacc6);
  gsl_interp_accel_free(kacc7);
  gsl_interp_accel_free(kacc8);
  gsl_interp_accel_free(kacc9);
  gsl_interp_accel_free(kacc10);
  gsl_interp_accel_free(kacc11);
  gsl_interp_accel_free(kacc12);
  gsl_interp_accel_free(kacc13);
  gsl_interp_accel_free(kacc14);
  gsl_interp_accel_free(kacc15);
  gsl_interp_accel_free(kacc16);
  gsl_interp_accel_free(kacc17);
  gsl_interp_accel_free(mgacc1);
  gsl_interp_accel_free(mgacc2);
  gsl_interp_accel_free(mgacc3);
  gsl_interp_accel_free(mgacc4);
  gsl_interp_accel_free(mgacc5);
  gsl_interp_accel_free(mgacc6);
  gsl_interp_accel_free(mgacc7);
  gsl_interp_accel_free(mgacc8);
  gsl_interp_accel_free(mgacc9);
  gsl_interp_accel_free(mgacc10);
  gsl_interp_accel_free(mgacc11);
  gsl_interp_accel_free(mgacc12);
  gsl_interp_accel_free(mgacc13);
  gsl_interp_accel_free(mgacc14);
  gsl_interp_accel_free(mgacc15);
  gsl_interp_accel_free(mgacc16);
  gsl_interp_accel_free(mgacc17);
  free(pdd_pmgs);
  free(pdt_pmgs);
  free(ptt_pmgs);
  free(a11_pmgs);
  free(a12_pmgs);
  free(a22_pmgs);
  free(a21_pmgs);
  free(a33_pmgs);
  free(b12_pmgs);
  free(b13_pmgs);
  free(b14_pmgs);
  free(b22_pmgs);
  free(b23_pmgs);
  free(b24_pmgs);
  free(b33_pmgs);
  free(b34_pmgs);
  free(b44_pmgs);
}
