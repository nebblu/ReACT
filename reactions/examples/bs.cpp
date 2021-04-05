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
#include <Copter/BSPT.h>

using namespace std;

/* Example code to output the real space bispectrum for modified gravity */

int main() {

	 //output file name
    const char* output = "check.dat";
    const char* cstr = "transfers/wmap9b";


    //const char* cstr = "transfers/test";

    /* Open output file */
    FILE* fp = fopen(output, "w");

    // Which gravity or dark energy model?
    // 1 for DGP, 2 for Hu-Sawicki,
    int mymodel = 2;

    // Keep it z=0 to keep Copter's Growth @ 1
    real z = 0;
    // Relative error in magnitude integrations
    real epsrel = 1e-3;

    int Nk =50;
    double kmin = 0.01;
    double kmax = 3.;

    Cosmology C(cstr);
    LinearPS P_l(C, z);
    HALO halo(C, P_l, P_l, P_l, epsrel);
    BSPT bspt(C, P_l, epsrel);
    SPT spt(C, P_l, epsrel);

    IOW iow;
    real b1,b2,b3,b4;

    double vars[4];


// 0: scale factor, 1: omega_total, 2: mg param (1e-10 ~ GR for default mg functions ), 3: =1 for DGP, =2 for f(R)
    vars[0] = 1./(1.+0.541);
    vars[1] = 0.279;
    vars[2] = 1e-5;

// initialise GR lin growth for PS normalisation
iow.inite2(vars[0],vars[1],vars[2], 1.,1.,mymodel);
spt.remp(fl_spt);
// initialise fitting function params (n_effective, k_nl, sigma_8)
bspt.mypspline(vars[0],vars[1],1);
double vars2[9];

// redefine inputs (should fix this!)
// 0 = omega0, 1 = mg1, 2=mg2, 3=mg3, 4=a, 5= 1 for DGP, 2 for Hu-Sawicki,
// 6 = integration accuracy, 7,8 = ODE solver absolute and relative accuracy
vars2[0] = vars[1];

vars2[1] = vars[2];
vars2[2] = 1.;
vars2[3] = 1.;

vars2[4] = vars[0];

vars2[5] = (double)mymodel;

vars2[6] = 1e-2;
vars2[7] = 1e-3;
vars2[8] = 1e-3;


 for(int i =0; i < Nk;  i ++) {

  real k =  kmin * exp(i*log(kmax/kmin)/(Nk-1));

// bispectra for equilateral config

 double mu = -0.5; // cosine of angle between k and k1
 double k1 = k;

 b1 = bspt.Btreen(vars2,k,k1,mu); // numerically computed tree level
 b2 = bspt.Bloopn(vars2,k,k1,mu); //numerically computed
 b3 = bspt.Bloop(1,k,k1,mu); //EdS approximated GR spectrum
 b4 = bspt.Bfit(k,k1,mu); // Gil-Marin/Namikawa et al fitting formula (GR and DGP only)


 printf("%e %e %e %e %e %e \n",k,b1,b2,b3,b4,b1/b2); // print to terminal
 fprintf(fp,"%e %e %e %e %e \n",k,b1,b2,b3,b4); // print to terminal


}

	/*close output file*/
    fclose(fp);
    return 0;
}
