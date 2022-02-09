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

/* Example code to output the 1-loop powerspectrum for modified gravity in real and redshift space*/
int main(int argc, char* argv[]) {

    /* Set parameters */
	// Input name of transfer_function/cosmology file.
    const char* cstr = "transfers/dgp";

	// Keep it z=0 to keep Copter's Growth @ 1
    real z = 0;

	// Relative error in magnitude integrations (and angular for some files)
    const real epsrel = 1e-3;

	// Number of k-modes you wish to output between kmin and kmax
    int Nk = 50;
    real kmin = 0.01;
    real kmax = 0.4;

	 //output file name
    const char* output = "DGP_P024_z0.dat";

    /* Open output file */
    FILE* fp = fopen(output, "w");

    /* Calculate SPT power spectrum */
    Cosmology C(cstr);
    LinearPS P_l(C, z);
    SPT spt(C, P_l, epsrel);
    IOW iow;
    NoWigglePS now(C, z , EisensteinHu);

  // output variables
  real p1,p2,p3,p4,p5,p6,p7,p8,p9;
  // Which gravity or dark energy model?
  // 1: GR  2: f(R) 3: DGP 4: quintessence 5: CPL
    int mymodel = 3;

  // chosen redshift
     double myz = 0.;
  // chosen omega_matter (total)
    double omega0 = 0.281;
  // MG parameter (for f(R) this is |fr0| and for nDGP this is Omega_rc)
     double mgpar = 0.25;

    // parameter values
    double vars[7];

    vars[0] = 1./(1.+myz);
    vars[1] =  omega0;
    vars[2] = mgpar;
    vars[3] = 1.; // wa for CPL
    vars[4] = 1.; // unusedd
    vars[5] = 30.; // number of mass bins
    vars[6] = 0.0; // omega_neutrinos

    double bias[3];
    bias[0] = 1.;
    bias[1] = 0.;
    bias[2] = 0.;

    // normalise the input P_L to initial times
    iow.inite(vars[0],vars[1],vars[2],vars[3],vars[4],mymodel);

    // Calculate (linear) velocity dispersion for the 1-loop RSD spectrum (or for TNS this is a free parameter)
    double pars[3];
    pars[0] = spt.sigmav_init(vars, mymodel);

    // choose RSD model :
    // 0: Kaiser
    // 1: TNS with q-bias (see 1507.01592 for example)
    // 2: TNS with Lagrangian bias (incomplete!)
    // 3: SPT 1-loop RSD spectrum

    double rsd_model = 3;

for(int i =0; i <Nk;  i ++) {

     real k = kmin * exp(i*log(kmax/kmin)/(Nk-1));

    /* for more observables check out the SPT.h file in the src directory */
      p1 = spt.PRSD_mg(rsd_model, 1, bias, vars, mymodel, pars, k, epsrel); // P0 SPT
      p2 = spt.PRSD_mg(rsd_model, 2, bias, vars, mymodel, pars, k, epsrel); // P2 SPT
      p3 = spt.PRSD_mg(rsd_model, 3, bias, vars, mymodel, pars, k, epsrel); // P4 SPT

     printf("%e %e %e %e  \n", k, p1,p2,p3); // print to terminal
     fprintf(fp,"%e %e %e %e  \n", k, p1,p2,p3); // print to file

}


	/*close output file*/
    fclose(fp);
    return 0;
}
