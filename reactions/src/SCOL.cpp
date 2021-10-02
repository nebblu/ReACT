#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_integration.h>
#include "SCOL.h"
#include "SpecialFunctions.h"
#include "Quadrature.h"


#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include <cmath>

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>
#include <cmath>
#include <random>


struct my_f_params { double a; double b; double c; };

static double delta_envc(double x, void * p)
{
  struct my_f_params * params = (struct my_f_params *)p;
  double OM = (params->a);
  double Z = (params->b);
  return 3.0 / 20.0 *
         pow(12*M_PI, 2.0/3.0) *
         (1. - 0.0123 * log10(1. + (1./OM - 1.) / pow(1.+Z,3.)));
}


static double Pzeta(double x, void * p)
{
  struct my_f_params * params = (struct my_f_params *)p;
  double beta = (params->a);
  double varomega = (params->b);
  double delta_envcr = (params->c);
  (void)(params); /* avoid unused parameter warning */
  return - pow(beta,varomega/2.0)/sqrt(M_PI*2.0)  *
         (1.0 + (varomega-1.0)*x / delta_envcr) *
         pow(1.0 - x/delta_envcr, -varomega/2.0-1.0) *
         exp(-(pow(beta,varomega)/2.0) * x * x / pow(1.0 - x/delta_envcr, varomega));
}



// EQ. A4 of  1812.05594 with a = e^t
// derivatives are taken with respect to t (or ln[a])
static int f_modscol(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    realtype y1, y2, ET, ET0, yenv, Fvir, hubble2, dhlnaoh, prefac, term1, term2, term4, IC, Rth, OM, OCB, T1, maxt;
    UserData data;

    data    = (UserData) user_data;
    IC      = data->IC;
    OM      = data->OM;
    OCB     = data->OCB;
    T1      = data->T1;
    Rth     = data->Rth;
    maxt    = data->maxt;

    y1 = Ith(y,1);
    y2 = Ith(y,2);

    ET      = SUNRexp(t - T1);
    ET0     = SUNRexp(T1);

    // make sure we don't exceed y_env spline range - above maximum t we don't care about the solution.
    // Should implement a better fix for this ...
    if (t>maxt){
            yenv = gsl_spline_eval (data->spline, maxt, data->acc);
                }
     else{
            yenv = gsl_spline_eval (data->spline, t, data->acc);
          }

    double myh = (y1 + ET)/ET;
    double myenv = (yenv + ET)/ET;

    Fvir = mymgF(ET*ET0, myh, myenv, Rth, OM, (data->par1),(data->par2),(data->par3), IC, (data->mymodel));
    // (H/H0)^2
    hubble2 = pow2(HAg(ET*ET0, OM, (data->par1), (data->par2), (data->par3),(data->mymodel)));
    // d H/ d ln[a] /H
    dhlnaoh = HA1g(ET*ET0,OM, (data->par1), (data->par2), (data->par3),(data->mymodel))/hubble2;

    prefac  =  SUNRexp(-THREE*t)*OCB; // a^3 Omega_cb

    term1   = -dhlnaoh*y2;
    term2   = (1.+dhlnaoh)*y1;

    term4   = -prefac*(ET+y1)*(
                                (IC+ONE)/(y1/ET + ONE)/(y1/ET + ONE)/(y1/ET + ONE) - ONE
                            ) *(Fvir + ONE)
            /hubble2*HALF;  // RHS of Eq. A4

    Ith(ydot,1) = y2; // y' = y2
    Ith(ydot,2) = term4 + term2 + term1; // y'' = RHS
    return(0);
}



// EQ.A4 of  1812.05594 with F =0 for y_env
static int fscol(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype y1, y2, ET, ET0, prefac, term1, term2, term4, IC, OM, OCB, T1;
  UserData data;

  data    = (UserData) user_data;
  IC      = data->IC;
  OM      = data->OM;
  OCB     = data->OCB;
  T1      = data->T1;
  realtype OL      = 1 - OM;
  y1      = Ith(y,1); y2 = Ith(y,2);

  ET      = SUNRexp(t - T1);
  ET0     = SUNRexp(T1);


  prefac  =  SUNRexp(-THREE*t)*OM;

  term1   =  THREE*prefac*y2*HALF/(prefac+OL);
  term2   = -y1*(prefac-TWO*OL)/(prefac+OL)*HALF;
  term4   = -prefac*(ET+y1)*(
                              (IC+ONE)/(y1/ET + ONE)/(y1/ET + ONE)/(y1/ET + ONE) - ONE
                            )
            /(prefac+OL)*HALF;

  Ith(ydot,1) = y2;
  Ith(ydot,2) = term4 + term2 + term1;
  return(0);
}



static int gscol(realtype t, N_Vector y, realtype *gout, void *user_data)
{
  realtype y1, ET, T1;
  UserData data;

  data    = (UserData) user_data;
  T1      = data->T1;
  y1      = Ith(y,1);
  ET      = SUNRexp(t - T1);
  gout[0] = ET + y1 - ONE;

  return(0);
}




double SCOL::maxP_zeta(double sig2, double dsig2dR, double OM, double Z)
{
  int status,status2;
  int iter = 0, max_iter = 100;
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;
  gsl_function F;
  gsl_function A1; // delta_envc(OM, Z)
  struct my_f_params sig2_params = {sig2, dsig2dR, sqrt(sig2)}; /* sig2, dsig2dR, sig_8*/
  struct my_f_params OM_Z_params = {OM, Z, 1.0/(1.0+Z)}; // OM, Z, a0=1/(1+Z)
  A1.function = &delta_envc;
  A1.params   = &OM_Z_params;
  double gamma = - 8.0 / 3.0 / (sig2_params.a) * (sig2_params.b);
  double d_envcr = GSL_FN_EVAL(&A1,0);
  double a = 0.0;
  double b = d_envcr*0.99;
  double m = 0.0001;
  double varomega = gamma * d_envcr;


  double betatest =pow(5.0/8.0, 3.0/d_envcr) / pow((sig2_params.c), 2.0/varomega);
  struct my_f_params Pzeta_params = {pow(5.0/8.0, 3.0/d_envcr) / pow((sig2_params.c), 2.0/varomega),
                                     varomega,
                                     d_envcr};
  F.function = &Pzeta;
  F.params   = &Pzeta_params;


  T = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc (T);

 gsl_set_error_handler_off();
 status2 = gsl_min_fminimizer_set (s, &F, m, a, b);
 if (status2 ){
   fprintf(stderr, "error: %s\n", gsl_strerror (status2));
   fprintf(stderr, "The offending values of sigma8, sigma8', z, omega_m, varomega, d_envcr, beta, Dl_spt, g_de : %e %e %e %e %e %e %e %e %e  \n", sig2, dsig2dR, Z, OM,varomega, d_envcr, betatest, Dl_spt, g_de );
   error.errorno = status2;
 }
gsl_set_error_handler (NULL);

  do
    {
      iter++;
      status = gsl_min_fminimizer_iterate (s);

      m = gsl_min_fminimizer_x_minimum (s);
      a = gsl_min_fminimizer_x_lower (s);
      b = gsl_min_fminimizer_x_upper (s);

      status
        = gsl_min_test_interval (a, b, 0.001, 0.0);

      if (status == GSL_SUCCESS)
      ;
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_min_fminimizer_free (s);
  return m;
}


// Integrand for growth (EQ. 35 of 1812.05594)
static double xE3(double x, void * p)
{
  struct my_f_params * params = (struct my_f_params *)p;
  double OM = (params->a);
  double OL = 1 - OM;
  return pow(x*sqrt(OM/x/x/x + OL),-3.0);
}



/// Edit function to solve DE with modified poisson.
double SCOL::Delta_Lambda(double OM, double Z)
{
  gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

  double result, error;
  struct my_f_params OM_Z_params = {OM, Z, 1.0/(1.0+Z)}; // OM, Z, a0=1/(1+Z)
  gsl_function F;
  F.function    = &xE3;
  F.params      = &OM_Z_params;

  gsl_integration_qags (&F, 0, OM_Z_params.c, 0, 1e-7, 1000,
                        w, &result, &error);

  gsl_integration_workspace_free (w);

  return 2.5 * OM
          * sqrt(OM/ OM_Z_params.c/ OM_Z_params.c / OM_Z_params.c  + 1-OM)
          / 3e-5 * result;
}


///////////////////////


// Solves for y_env (i.e. F=0 in Eq. A4 of 1812.05594 )
// stores results in xxyy
// parameters: omega_matter, -log(1.0 + myz1) , m1/d1, array to store y/y_i
int SCOL::yenv(double OM_REAL, double OCB_REAL, double XF, double delta_envi, arrays_T xxyy)
{
  realtype reltol, t, tout, T1;
  N_Vector y, abstol;
  SUNMatrix A;
  SUNLinearSolver LS;
  UserData data;
  void *cvode_mem;
  int flag, flagr, iout;
  int rootsfound[2];

      y             = abstol = NULL;
      data          = NULL;
      A             = NULL;
      LS            = NULL;
      cvode_mem     = NULL;

      y             = N_VNew_Serial(NEQ);
      abstol        = N_VNew_Serial(NEQ);
      reltol        = RTOL;
      Ith(abstol,1) = ATOL1;
      Ith(abstol,2) = ATOL2;

      data          = (UserData) malloc(sizeof *data);  /* Allocate data memory */
      data->IC      = delta_envi;  /* this parameter changes the function everytime, and has to be in the loop*/
      data->OM      = OM_REAL;
      data->OCB     = OCB_REAL;

      T1 = T00;

      data->T1      = T1;

      Ith(y,1)      = Y1;
      Ith(y,2)      = -data->IC/THREE; /* this is the initial condition for y' */


      cvode_mem     = CVodeCreate(CV_ADAMS);
      flag          = CVodeInit(cvode_mem, fscol, T1, y);
      flag          = CVodeSVtolerances(cvode_mem, reltol, abstol);
      flag          = CVodeSetUserData(cvode_mem, data);
      flag          = CVodeRootInit(cvode_mem, 1, gscol);
      A             = SUNDenseMatrix(NEQ, NEQ);
      LS            = SUNDenseLinearSolver(y, A);
      flag          = CVDlsSetLinearSolver(cvode_mem, LS, A);

      iout = 0;  tout = 20;

      (*xxyy).xx[iout] = T1;
      (*xxyy).yy[iout] = 0.0;


      while(1)
      {
        flag = CVode(cvode_mem, tout, y, &t, CV_ONE_STEP);
//         printf ("%.7f %.7f \n", t, Ith(y,1));
        iout += 1;
        (*xxyy).count = iout+1;
        (*xxyy).xx[iout] = t;
        (*xxyy).yy[iout] = Ith(y,1); // storing y_env
        if (flag == CV_ROOT_RETURN)
            {
              flagr = CVodeGetRootInfo(cvode_mem, rootsfound);
              // if (check_flagscol(&flagr, "CVodeGetRootInfo", 1)) free(xxyy);
              if (check_flagscol(&flagr, "CVodeGetRootInfo", 1)) return(1);
              goto out;
            }

        if (check_flagscol(&flag, "CVode", 1))
            {
                goto out;
            }
        if (t >= 2.0 ) goto out;
      }

      out:

      /* Free y and abstol vectors */
      N_VDestroy(y);
      N_VDestroy(abstol);

      /* Free integrator memory */
      CVodeFree(&cvode_mem);

      /* Free the linear solver memory */
      SUNLinSolFree(LS);

      /* Free the matrix memory */
      SUNMatDestroy(A);

      free(data);
      // break;
    // free(xxyy);
    // printf("%g\n", data->IC);
    // return xxyy;
      return 0;
}



/*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 */

void SCOL::PrintOutput(realtype t, realtype y1, realtype y2) /*, realtype y3)*/
{
  #if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("%0.9Le     %14.9Le  %14.9Le  \n", t, y1, y2); /* , y3); */
  #elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("%0.9e      %14.9e  %14.9e  \n", t, y1, y2); /* , y3); */
  #else
    printf("%0.9e      %14.9e  %14.9e  \n", t, y1, y2); /* , y3); */
  #endif
    return;
}

void SCOL::PrintRootInfo(int root_f1)
{
  printf("    rootsfound[] = %3d \n", root_f1);

  return;
}

int check_flagscol(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
        funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
          funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
        funcname);
    return(1); }

  return(0);
}



// Solves for delta_i
int SCOL::SphericalCollapse(double *dC, arrays_T3 xxyyzz, UserData data_vec, double TMULT_REAL, double delta_g)
{

    realtype reltol, t, tout, delta_uu, delta_ll, T1;
    N_Vector y, abstol;
    SUNMatrix A;
    SUNLinearSolver LS;
    UserData data;
    void *cvode_mem;
    int flag, flagr, iout, ioutout;
    int count_of_bisec;
    int rootsfound[2];

    delta_uu =  delta_g;
    delta_ll =  4.*delta_g;

// 32 bisections should be enough for most cases, but increase if needed
    int noofbi = 32;
    int noofbi2 = noofbi-1;
  for (count_of_bisec=0; count_of_bisec < noofbi; ++count_of_bisec)
  {
      y             = abstol = NULL;
      data          = NULL;
      A             = NULL;
      LS            = NULL;
      cvode_mem     = NULL;


      y             = N_VNew_Serial(NEQ);
      abstol        = N_VNew_Serial(NEQ);
      reltol        = RTOL;
      Ith(abstol,1) = ATOL1;
      Ith(abstol,2) = ATOL2;

      data          = (UserData) malloc(sizeof *data);  /* Allocate data memory */
      if(check_flagscol((void *)data, "malloc", 2)) { return 1 ;}
      data->Rth     = (*data_vec).Rth;
      data->IC      = (delta_uu+delta_ll)*HALF;  /* this parameter changes the function everytime, and has to be in the loop*/
      data->spline  = (*data_vec).spline;
      data->acc     = (*data_vec).acc;
      data->T1      = (*data_vec).T1;
      data->par1    = (*data_vec).par1;
      data->par2    = (*data_vec).par2;
      data->par3    = (*data_vec).par3;
      data->maxt    = (*data_vec).maxt;
      data->mymodel = (*data_vec).mymodel;
      data->OM      = (*data_vec).OM;
      data->OCB     = (*data_vec).OCB;
      T1            = (*data_vec).T1;


      Ith(y,1)      = Y1;
      Ith(y,2)      = -data->IC/THREE; /* this is the initial condition for y' */


      cvode_mem     = CVodeCreate(CV_ADAMS);
      flag          = CVodeInit(cvode_mem, f_modscol, T1, y);
      flag          = CVodeSVtolerances(cvode_mem, reltol, abstol);
      flag          = CVodeSetUserData(cvode_mem, data);
      flag          = CVodeRootInit(cvode_mem, 1, gscol);
      A             = SUNDenseMatrix(NEQ, NEQ);
      LS            = SUNDenseLinearSolver(y, A);
      flag          = CVDlsSetLinearSolver(cvode_mem, LS, A);

      iout = 0;  tout = -10.; ioutout = 0;

      (*xxyyzz).xx[ioutout] = (*data_vec).T1;
      (*xxyyzz).yy[ioutout] = Y1;
      (*xxyyzz).zz[ioutout] = 0.0;

      while(1)
      {
        if (count_of_bisec == noofbi2)
        {
          flag = CVode(cvode_mem, tout, y, &t, CV_ONE_STEP);
          ioutout += 1;
          (*xxyyzz).count = ioutout+1;
          (*xxyyzz).xx[ioutout] = t;
          (*xxyyzz).yy[ioutout] = Ith(y,1);
          (*xxyyzz).zz[ioutout] = Ith(y,2);
        }
        else
        {
          flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
        }

        if (flag == CV_ROOT_RETURN)
            {
              flagr = CVodeGetRootInfo(cvode_mem, rootsfound);
              if (check_flagscol(&flagr, "CVodeGetRootInfo", 1)) { return 1; };

              delta_ll = data->IC;

              goto out;
            }

        if (check_flagscol(&flag, "CVode", 1))
            {
                delta_ll = data->IC;
                goto out;
            }

          if (flag == CV_SUCCESS)
              {
                iout++;
                tout *= TMULT_REAL;
              }

          if (count_of_bisec == noofbi2)
          {
            if (t > 0.) { goto nout;}
          }
          else
          {
            if (iout == NOUT_DC) goto nout;
          }
        }

        nout:
          delta_uu = data->IC;

        out:
        /* Free y and abstol vectors */
        N_VDestroy(abstol);

        /* Free integrator memory */
        CVodeFree(&cvode_mem);

        /* Free the linear solver memory */
        SUNLinSolFree(LS);

        /* Free the matrix memory */
        SUNMatDestroy(A);
        *dC = data->IC;
	free(data);
	N_VDestroy(y);
 }

// if this pops up, we should increase number of bisections.
    if (SUNRabs(delta_uu - delta_ll) >= EPSILON)
    {
      printf("Wrong Accuracy d_uu - d_ll = %g, epsilon = %g \n", delta_uu - delta_ll, EPSILON);
    }

      return  0;
}




/* Solver for a_vir (eq A6 of 1812.05594), delta_avir and delta_acollapse */
//myscolparams[] is storage for the above quantities (0: delta_collapse, 1: a_virial, 2: delta_virial)

// needs omega0, initial top-hat radius , initial density of enviornment, modification parameters and selection of theory
// model = 1 gives GR collapse, 2+ gives modified collapse - see SpecialFunctions.cpp
double SCOL::myscol(double myscolparams[], double acol, double omegacb, double omeganu, double Rthp, double sig1, double sig2, double pars[], int model)
{

        double omega0 = omegacb + omeganu;

      // theory params
        double p1 = pars[0];
        double p2 = pars[1];
        double p3 = pars[2];

        for(int i=0; i<3; i++){
        myscolparams[i] =0.;
        }

        // Used to initialise y_env, y_h for all z in order to spline and find a_vir
       double myz = 1./acol - 1.;
       double XF = -log(1.0 + myz);

       double m = maxP_zeta (sig1, sig2, omegacb, myz);
       // linear growth divided by initial scale factor for total matter poisson (pseudo) or cb poisson (real)
       double d;
       if (model==1 && omega0 == omegacb) {
         d = Delta_Lambda (omega0, myz);
       }
        else{
          d = Dl_spt/3e-5;
         }


       arrays_T xxyy = (arrays_T)malloc( sizeof(struct arrays) );

       yenv (omega0, omegacb, XF , m/d, xxyy);

       // set max scalefactor
       double maximumt = (*xxyy).xx[(*xxyy).count-1];

       gsl_interp_accel *acc = gsl_interp_accel_alloc ();
       gsl_spline *myspline_yenv    = gsl_spline_alloc (gsl_interp_cspline, (*xxyy).count);
       gsl_spline_init (myspline_yenv, (*xxyy).xx, (*xxyy).yy, (*xxyy).count);

      // we progressively solve for final times from af = 0.000045 (Exp[-10]) to acolapse in steps of TMULT, over NOUT_DC steps.
       double TMULT  = pow(-XF/10.,ONE/(NOUT_DC-1.));

       // holders for delta_c and for a, y(a), y'(a)
       arrays_T3 xxyyzz = (arrays_T3)malloc( sizeof(struct arrays3D) );

         double mydelta;

         if (model == 1) {
           //  take standard guess if GR collapse
           mydelta = DELTA1/acol;
         }
         else{
           // If modified gravity is active, set initial condition for SC to 10% higher than y_{env,initial} if we need to use y_env in spherical collapse as in f(R)
           // This guess allows us to solve for extreme cases (e.g. fr0=10^{-4}, m_nu>0.3eV).
           mydelta = m/d*1.1/acol;
         }

         UserData data;
         data          = (UserData) malloc(sizeof(struct usdat));  /* Allocate data memory */
         data->IC      = (DELTA1+DELTA2)*HALF;  /* this parameter changes the function everytime, and has to be in the loop*/
         data->spline  = myspline_yenv;
         data->acc     = acc;
         data->T1      = T00; // = Log[ai]
         data->OM      = omega0;
         data->OCB     = omegacb;
         data->par1    = p1;
         data->par2    = p2;
         data->par3    = p3;
         data->Rth     = Rthp;
         data->mymodel = model;
	       data->maxt    = maximumt;

          // initialse the delta_i (solver spherical collapse differential equation)
        SphericalCollapse (&myscolparams[0], xxyyzz, data, TMULT, mydelta);
       	free(data);

     double ainit = exp(T00);

     arrays_T3 myen = (arrays_T3)malloc( sizeof(struct arrays3D) );
     arrays_T3 myamax = (arrays_T3)malloc( sizeof(struct arrays3D) );

      // populate arrays for total energy and scale factor lower bound used in root finder
     for(int i = 0; i< (*xxyyzz).count; i++ ){

       double ai = exp((*xxyyzz).xx[i]);
       double arat = ainit/ai;

       // reparametrise to y and y' in the appendix A of  1812.05594

       // Halo y
       double myy = ((*xxyyzz).yy[i] + 1./arat)*arat;
       double myp = ((*xxyyzz).zz[i] + 1./arat)*arat - myy;
       // Enviornment y
       double yenvq = gsl_spline_eval(myspline_yenv, (*xxyyzz).xx[i], acc);
       double myyenv = (yenvq+1./arat)*arat;

       // delta as given in Eq.34 of 1812.05594
       double mydelt =(1.+ myscolparams[0])/pow3(myy) -1.;

       // energy contributions

      // potential energy (always the same):
       double prefac = -omegacb*pow2(myy/ainit)/ai;

       // newtonian contribution
       double wn =  prefac*(1.+mydelt);
        // modified gravity, dark energy  and Kinetic energy contributions
       double wphi, weff, ke;
       if(model==1){
         wphi = 0.;
         weff = 2.*(1.-omega0)*pow2(myy/arat);
         ke =  pow2(HA(ai, omega0)*(myy + myp)/arat);
       }
       else{
         wphi = prefac * mymgF(ai, myy, myyenv, Rthp, omegacb, p1, p2, p3, myscolparams[0],model)*mydelt;
         weff = WEFF(ai,omegacb,p1,p2,p3,model)*pow2(myy/arat);
         ke =  pow2(HAg(ai, omega0,p1,p2,p3,model)*(myy + myp)/arat); // might need to change to cb TO CHECK
       }


        // RHS of virial theorem
       (*myen).xx[i] = ai;
       (*myen).yy[i] = 2.*ke + wn + wphi + weff;

        // used in solver for amax just below to get 2nd root of virial theorem (sets lower bound for solve search )
       (*myamax).xx[i] = ai;
       (*myamax).yy[i] = (*xxyyzz).zz[i] + 1./arat;

         }

      // total energy spline
     gsl_interp_accel *acc2 = gsl_interp_accel_alloc ();
     gsl_spline *energytot    = gsl_spline_alloc (gsl_interp_cspline, (*xxyyzz).count);
     gsl_spline_init (energytot, (*myen).xx, (*myen).yy, (*xxyyzz).count);



     // Create spline of y_h(delta_i)
     // y(ln(a)) spline
     gsl_interp_accel *acc3 = gsl_interp_accel_alloc ();
     gsl_spline *myspline_halo    = gsl_spline_alloc (gsl_interp_cspline, (*xxyyzz).count);
     gsl_spline_init (myspline_halo, (*xxyyzz).xx, (*xxyyzz).yy, (*xxyyzz).count);



     // y(ln(a)) - a/ai spline
     gsl_interp_accel *acc4 = gsl_interp_accel_alloc ();
     gsl_spline *myspline_amax    = gsl_spline_alloc (gsl_interp_cspline, (*xxyyzz).count);
     gsl_spline_init (myspline_amax, (*myamax).xx, (*myamax).yy , (*xxyyzz).count);



/// Find scale factor such that y' = -a/ai (i.e., R_tophat' = 0)
double amax;
// set maximum  and minimum scale factor for which y is actually solved.
double alimit1 = ((*myamax).xx[(*xxyyzz).count-1]);
double alimit2 = ((*myamax).xx[0]);


std::mt19937 gen0(2);

std::uniform_real_distribution<> dis0(alimit2, alimit1);

double pos_pt0 = dis0(gen0);
double neg_pt0 = dis0(gen0);

while ( gsl_spline_eval (myspline_amax, pos_pt0, acc4) < 0.){
    pos_pt0 = dis0(gen0);}

while (gsl_spline_eval (myspline_amax, neg_pt0,acc3)  > 0.){
    neg_pt0 = dis0(gen0);}

 const double about_zero_mag0 = 1e-4;
for (;;)
{
    double mid_pt0 = (pos_pt0 + neg_pt0)/2.0;
    double f_mid_pt0 =  gsl_spline_eval (myspline_amax, mid_pt0,acc4) ;


    if (fabs(f_mid_pt0)  < about_zero_mag0){
       // assign a_vir
        amax = mid_pt0;
        break;
          }

    if (f_mid_pt0 >= 0.){
        pos_pt0 = mid_pt0;}
    else{
        neg_pt0 = mid_pt0;}
}

// Now solve virial theorem

    std::mt19937 gen(2);

		std::uniform_real_distribution<> dis(amax, alimit1);

		double pos_pt = dis(gen);
		double neg_pt = dis(gen);

// find a  2T+W = 0
		while ( gsl_spline_eval (energytot, pos_pt,acc2) < 0.){
				pos_pt = dis(gen);}

		while ( gsl_spline_eval (energytot, neg_pt,acc2)  > 0.){
				neg_pt = dis(gen);}

		 const double about_zero_mag = 1e-4;
		for (;;)
		{
				double mid_pt = (pos_pt + neg_pt)/2.0;
				double f_mid_pt =   gsl_spline_eval (energytot, mid_pt, acc2) ;


				if (fabs(f_mid_pt)  < about_zero_mag){
           // assign a_vir
            myscolparams[1] = mid_pt;
            break;
							}

				if (f_mid_pt >= 0.){
						pos_pt = mid_pt;}
				else{
						neg_pt = mid_pt;}
		}

    // delta (a_vir) (non-linear)  [See eq. 34 of 1812.05594]
     myscolparams[2] = (1.+ myscolparams[0])/pow3((gsl_spline_eval(myspline_halo, log(myscolparams[1]), acc3)*ainit/myscolparams[1] + 1.)) - 1.;

    // delta_acol  (linearly extrapolated.)
     myscolparams[0] = d*myscolparams[0];

    // free all memory
    gsl_spline_free (myspline_halo);
    gsl_spline_free (myspline_yenv);
    gsl_spline_free (energytot);
    gsl_spline_free (myspline_amax);
    free(myamax);
    free(xxyy);
    free(xxyyzz);
    free(myen);
    gsl_interp_accel_free (acc);
    gsl_interp_accel_free (acc2);
    gsl_interp_accel_free (acc3);
    gsl_interp_accel_free (acc4);

    return 0;

}
