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

#include <omp.h>

#include <Copter/HALO.h>
#include <Copter/Cosmology.h>
#include <Copter/LinearPS.h>
#include <Copter/SpecialFunctions.h>
#include <Copter/SPT.h>

using namespace std;

using std::ifstream;
using std::string;
using std::istringstream;


/* Example code to output the halo model powerspectrum for modified gravity */

int main(int argc, char* argv[]) {
  // Which gravity or dark energy model?
  // 1: GR  2: f(R) 3: DGP 4: quintessence 5: CPL
  int mymodel = 2;

  // Modified gravity active?
  bool modg = true;
  // Is the transfer being fed to ReACT of the target cosmology? If false, the transfer should be LCDM at z=0.
  bool mgcamb = false;

	 //output file name
    const char* output = "ps_f5_z0.dat";
    const char* cstr = "transfers/Matteo_fr";
// 0: scale factor, 1: omega_total, 2-4: mg param (1e-10 ~ GR for default mg functions ), 5: number of points in halo-mass loop in scol_init , 30 works well.
double vars[8];

  // chosen redshift
    double myz = 0.;
  // chosen omega_matter (total)
    double omega0 = 0.3072;
  // MG parameter (for f(R) this is |fr0|)
    double mgpar = 1e-5;

    vars[0] = 1./(1.+myz);
    vars[1] =  omega0;
    vars[2] = mgpar;
    vars[3] = 1.; // wa for CPL
    vars[4] = 1.; // unusedd
    vars[5] = 50.; // number of mass bins
    vars[6] = 0.0; // omega_neutrinos

    /* Open output file */
    FILE* fp = fopen(output, "w");

    // Keep it z=0 to keep Copter's Growth @ 1
    real z = 0;
    // Relative error in magnitude integrations
    real epsrel = 1e-3;

    Cosmology C(cstr);
    LinearPS P_l(C, z);
    HALO halo(C, P_l,P_l,P_l,P_l, epsrel);
    SPT spt(C, P_l, epsrel);

    IOW iow;
    real p1,p2,p3,p4,p5,p6;


// initialise wCDM/LCDM lin growth for PS normalisation
iow.initnorm(vars,mymodel);
/// initialise delta_c(M), a_vir(M), delta_avir(M) and v(M)
halo.scol_init(vars,mgcamb,mymodel);
halo.scol_initp(vars,mgcamb);
halo.react_init_nu(vars,mgcamb,modg,mymodel);

//#pragma omp parallel for
int Nk =100;
double kmin = 0.001;
double kmax = 10.;

 for(int i =0; i < Nk;  i ++) {

  real k = kmin * exp(i*log(kmax/kmin)/(Nk-1));

      p1 =  halo.one_halo(k, vars);
      p2 =  halo.one_halop(k, vars);
      p3 =  halo.reaction_nu(k, vars, mgcamb);

     printf("%d %e %e %e %e \n", i, k, p1,p2,p3); // print to terminal
     fprintf(fp,"%e %e %e %e \n", k, p1,p2, p3); // print to file

}

	/*close output file*/
    fclose(fp);
    return 0;
}
