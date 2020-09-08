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


/* Example code to output the halo model powerspectrum for modified gravity */

int main(int argc, char* argv[]) {

	 //output file name
    const char* output = "ps_f5_z0.dat";

    const char* cstr = "transfers/dgp";

// 0: scale factor, 1: omega_total, 2-4: mg param (1e-10 ~ GR for default mg functions ), 5: number of points in halo-mass loop in scol_init , 30 works well.
double vars[6];

    vars[0] = 1.;
    vars[1] = 0.3072;
    //vars[2] = 1e-5;
    vars[2] = 0.1;
    vars[3] = 1.;
    vars[4] = 1.;
    vars[5] = 30.;
    int model = 3;

    /* Open output file */
    FILE* fp = fopen(output, "w");

    // Keep it z=0 to keep Copter's Growth @ 1
    real z = 0;
    // Relative error in magnitude integrations
    real epsrel = 1e-2;

    Cosmology C(cstr);
    LinearPS P_l(C, z);
    HALO halo(C, P_l, epsrel);
    SPT spt(C, P_l, epsrel);

    IOW iow;
    real p1,p2,p3,p4,p5,p6;


// initialise wCDM/LCDM lin growth for PS normalisation
iow.initnorm(vars);
/// initialise delta_c(M), a_vir(M), delta_avir(M) and v(M)
halo.scol_init(vars,model);
halo.scol_initp(vars,model);
halo.react_init(vars,model);

//#pragma omp parallel for
int Nk =100;
double kmin = 0.0001;
double kmax = 10.;

 for(int i =0; i < Nk;  i ++) {

  real k =  kmin * exp(i*log(kmax/kmin)/(Nk-1));

      p1 =  halo.one_halo(k, vars);
      p2 =  halo.one_halop(k, vars);

      p3 =  halo.reaction(k, vars);

     printf("%d %e %e %e %e \n", i, k, p1,p2,p3); // print to terminal
     fprintf(fp,"%e %e %e %e \n", k, p1,p2, p3); // print to file

}

	/*close output file*/
    fclose(fp);
    return 0;
}
