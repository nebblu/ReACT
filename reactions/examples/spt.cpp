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

#include <Copter/Cosmology.h>
#include <Copter/LinearPS.h>
#include <Copter/SPT.h>
#include <Copter/GrowthFunction.h>
#include <Copter/SpecialFunctions.h>
#include <Copter/NoWigglePS.h>

using namespace std;
vector<vector<double> > allDatamult;

/* Example code to output the 1-loop powerspectrum for modified gravity in real and redshift space*/

int main(int argc, char* argv[]) {

    /* Set parameters */
	// Input name of transfer_function/cosmology file.
    const char* cstr = "transfers/fr_new";

	// Keep it z=0 to keep Copter's Growth @ 1
    real z = 0;

	// Relative error in magnitude integrations (and angular for some files)
    const real epsrel = 1e-3;

	// Number of k-modes you wish to output between kmin and kmax
    int Nk = 500;
    real kmin = 0.0001;
    real kmax = 0.5;

	 //output file name
    const char* output = "fr_new.dat";

    /* Open output file */
    FILE* fp = fopen(output, "w");

    /* Calculate SPT power spectrum */
    Cosmology C(cstr);
    LinearPS P_l(C, z);
    SPT spt(C, P_l, epsrel);
    IOW iow;
    NoWigglePS now(C, z , EisensteinHu);

// output variables
real p1,p2,p3,p4,p5;

// parameter values
double vars[5];

vars[0] = 1.; // scale factor
vars[1] = 0.24; // omega_matter
vars[2] = 1e-14; // mg param 1 (fr0, omega_rc, w0 --- see src/SpecialFunctions.cpp)
vars[3] = 1.; // mg param 2 (unused)
vars[4] = 1.;// mg param 3 (unused)
// model; 1: GR, 2: f(R) Hu Sawicki, n=1  3: DGP normal branch
int model=1;


// normalise the growth
iow.inite(vars[0],vars[1],vars[2],vars[3],vars[4],model);
spt.remp(fl_spt);


// Lagrangian bias params (See 1607.03149 for example) OR Q-bias params (see PRSD_mg in src/SPT.cpp).
double bias[3];
bias[0] = 1.; // b1
bias[1] = 0.; // b2
bias[2] = 0.; // stochasticity

double sigmav = 5.5; // velocity dispersion
double err = 1e-2; //  absolute error in differential equation solver for numerical PT kernels


for(int i =0; i <Nk;  i ++) {

    real k = kmin * exp(i*log(kmax/kmin)/(Nk-1));

/* for more observables check out the SPT.h file in the src directory */

  p1 = spt.PLOOPn2(0, vars, model, k, 1e-3); // linear spectrum
  p2 = spt.PLOOPn2(1, vars, model, k, 1e-3); // 1-loop spectrum
  p3 = spt.PRSD_mg(2,1,bias,vars, model, sigmav,k,err); // TNS monopole
  p4 = spt.PRSD_mg(2,2,bias,vars, model, sigmav,k,err); // TNS quadrupole
  p5 = spt.PRSD_mg(2,3,bias,vars, model, sigmav,k,err); // TNS hexdecapole


     printf("%e %e %e %e %e %e \n", k, p1,p2,p3,p4,p5); // print to terminal
     fprintf(fp,"%e %e %e %e %e %e  \n", k, p1,p2,p3,p4,p5 ); // print to file

  }


	/*close output file*/
    fclose(fp);
    return 0;
}
