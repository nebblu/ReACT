#ifndef SCOL_H
#define SCOL_H
#include "Common.h"
// SP_IVP.h
// Spherical Collapse IVP
//#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <cvode/cvode_direct.h>        /* access to CVDls interface            */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>    /* contains the macros ABS, SUNSQR, EXP */

#define Ith(v,i)    NV_Ith_S(v,i-1)         /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) SM_ELEMENT_D(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */

#define NEQ     2                /* number of equations  */
#define Y1      RCONST(0.0)      /* initial y components */
#define RTOL    RCONST(1.0e-5)   /* scalar relative tolerance            */

#define ATOL1   RCONST(1.0e-11)   /* vector absolute tolerance components */
#define ATOL2   RCONST(1.0e-11)

#define T00     RCONST(-10.41431317630211772495840705232694745063781738281250)      /* initial time           */
#define NOUT    1024               /* number of output times */
#define NOUT_DC    10               /* number of output times */

#define HALF    RCONST(0.5)      /* 0.5                    */
#define ONE     RCONST(1.0)      /* 1.0                    */
#define TWO     RCONST(2.0)
#define THREE   RCONST(3.0)       /* 3.0                    */
#define DELTA1  RCONST(3e-5)   /* Delta_i1  = ai */
#define DELTA2  RCONST(0.00012)   /* Delta_i2 = 4 x ai */
#define EPSILON RCONST(1.0e-9)
#define ZERO    RCONST(0.0)
#define NINE    RCONST(9.0)
#define TEN     RCONST(10.0)
#define Gnewton RCONST(4.302e-09)
#define coef    RCONST(8987404.41) // 1/H0^2


// structures to store spherical collapse calculations
typedef struct arrays{
  int count;
  double xx[1002];
  double yy[1002]; } *arrays_T;


typedef struct arrays3D{
  int count;
  double xx[1002];
  double yy[1002];
  double zz[1002]; } *arrays_T3;


typedef struct usdat {
  realtype IC;
  realtype OM;
  realtype OCB;
  realtype Rth;
  realtype T1;
  double par1;
  double par2;
  double par3;
  double maxt;
  int mymodel;
  gsl_spline *spline;
  gsl_interp_accel *acc;
} *UserData;

typedef struct scol_error{
  int errorno = 0;
} scol_error_T;

extern int check_flagscol(void *flagvalue, const char *funcname, int opt); //


class SCOL {
public:

// used functions in example file
   double maxP_zeta(double sig2, double dsig2dR, double OM, double Z);
   double Delta_Lambda(double OM, double Z); //
   // solves for y_enviornment
   int yenv(double OM_REAL,  double OCB_REAL, double XF, double delta_envi, arrays_T xxyy); // gives the environmental dependence of spherical collapse
   // solves for y_halo
   int SphericalCollapse(double *dC, arrays_T3 xxyyzz, UserData data_vec, double TMULT_REAL, double delta_g); // spherical collapse solver
   // solves for a_virial
   double myscol(double myscolparams[], double acol, double omegacb, double omeganu, double Rthp, double sig1, double sig2, double pars[], int model); // solves for virial quantities and stores them in array myscolparams


   void PrintOutput(realtype t, realtype y1, realtype y2);
   void PrintRootInfo(int root_f1); //

   double funcscol(double xi, void *user_data); //

   scol_error_T error;
};

#endif
