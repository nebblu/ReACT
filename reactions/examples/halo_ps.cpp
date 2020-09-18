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

vector<vector<double> > mypk;

/* Example code to output the halo model powerspectrum for modified gravity */

int main(int argc, char* argv[]) {

	 //output file name
    const char* output = "ps_f5_z0.dat";

    const char* cstr = "transfers/Matteo_fr";

// 0: scale factor, 1: omega_total, 2-4: mg param (1e-10 ~ GR for default mg functions ), 5: number of points in halo-mass loop in scol_init , 30 works well.
double vars[7];

    vars[0] = 1./(1.+1.);
    vars[1] = 0.3072;
    vars[2] = 1e-5;
    vars[3] = 1.;
    vars[4] = 1.;
    vars[5] = 30.;
    vars[6] = 0.0;

    /* Open output file */
    FILE* fp = fopen(output, "w");

    // Keep it z=0 to keep Copter's Growth @ 1
    real z = 0;
    // Relative error in magnitude integrations
    real epsrel = 1e-2;

    Cosmology C(cstr);
    LinearPS P_l(C, z);
    HALO halo(C, P_l,P_l,P_l,P_l, epsrel);
    SPT spt(C, P_l, epsrel);

    IOW iow;
    real p1,p2,p3,p4,p5,p6;


// initialise wCDM/LCDM lin growth for PS normalisation
iow.initnorm(vars);
/// initialise delta_c(M), a_vir(M), delta_avir(M) and v(M)
halo.scol_init(vars,false);
halo.scol_initp(vars,false);
halo.react_init(vars);

// load in k-binning from sims
ifstream fin2("validate/matteo_data/reactions/HM_reaction_standard_HMF_F5_z1.dat");


// Load in the data
string line2;
    while (getline(fin2, line2)) {      // for each line
            vector<double> lineData;           // create a new row
            double val;
            istringstream lineStream(line2);
            while (lineStream >> val) {          // for each value in line
                    lineData.push_back(val);           // add to the current row
            }
            mypk.push_back(lineData);         // add row to allData
    }

    int Nk =  mypk.size();


//#pragma omp parallel for
//int Nk =100;
double kmin = 0.0001;
double kmax = 10.;

 for(int i =0; i < Nk;  i ++) {

  real k =  mypk[i][0];//kmin * exp(i*log(kmax/kmin)/(Nk-1));

      p1 =  halo.one_halo(k, vars);
      p2 =  halo.one_halop(k, vars);

      p3 =  halo.reaction(k, vars)/mypk[i][1];

     printf("%d %e %e %e %e \n", i, k, p1,p2,p3); // print to terminal
     fprintf(fp,"%e %e %e %e \n", k, p1,p2, p3); // print to file

}

	/*close output file*/
    fclose(fp);
    return 0;
}
