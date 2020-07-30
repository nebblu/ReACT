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

//using namespace std;
using std::ifstream;
using std::string;
using std::istringstream;

vector<vector<double> > mytrans;
vector<vector<double> > mypk;


/* Example code to output the reaction and halo spectra for mg + neutrinos */
int main(int argc, char* argv[]) {


// Load transfer function at z from MGCAMB with all species at some redshift
ifstream fin("validate/mgcamb/transfer_out_z1.dat");

// Load in the transfer data
string line;
    while (getline(fin, line)) {      // for each line
            vector<double> lineData;           // create a new row
            double val;
            istringstream lineStream(line);
            while (lineStream >> val) {          // for each value in line
                    lineData.push_back(val);           // add to the current row
            }
            mytrans.push_back(lineData);         // add row to allData
    }

// Populate arrays and normalise to 1 at large scales for Cosmology class input
int Nkr = mytrans.size();
int* Nkt = &Nkr;
array Tm(*Nkt);
array Tcb(*Nkt);
array Tnu(*Nkt);
array ki(*Nkt);
          for(int i = 0; i< Nkr; i++){
                  ki[i] = mytrans[i][0];
                  Tm[i] = mytrans[i][6]/mytrans[0][6];
                  Tcb[i] = mytrans[i][7]/mytrans[0][7];
                  Tnu[i] = mytrans[i][5]/mytrans[0][5]; // should be adjusted if this column is all 0s ....
              }


// integration error
real epsrel = 1e-2;

// // // Specify params
double h  = 0.6898;
double n_s = 0.969; // spectral index
double Omega_m =  0.281474; // total matter fraction
double Omega_b  = 0.0473; //  baryon fraction
double Omega_nu = 0.00902645;  // neutrino fraction

// z=0
// double sigma_8m = 0.7296; // IMP
// double sigma_8cb = 0.7434; // IMP
// double sigma_8nu = 0.3271; // IMP

// z=1
double sigma_8m = 0.4507; // IMP
double sigma_8cb = 0.4602; // IMP
double sigma_8nu =  0.1730; // IMP


double modg = 1e-15;  //  modified gravity param
double massb = 30.; // number of mass bins between 5<Log10[M]<18
// // input transfer redshift
double myz = 1.;


// if false, Omega_nu should be 0 --- create flag.
bool mgcamb = true; // mgcamb transfer input or camb transfer input -- maybe rename to massnu because it's more correctly massive neutrinos or not

// LCDM transfer at z=0 for Copter growth --- current ReACT implementation (for mgcamb=false)
// const char* cstr = "tests/nu0/test";
// Cosmology C(cstr);
// LinearPS P_l(C, 0.);
// HALO halo(C, P_l, P_l, P_l, epsrel);

// Load cosmology classes
Cosmology Cm(h, n_s, Omega_m, Omega_b, sigma_8m, ki, Tm);
Cosmology Ccb(h, n_s, Omega_m, Omega_b, sigma_8cb, ki, Tcb);
Cosmology Cnu(h, n_s, Omega_m, Omega_b, sigma_8nu, ki, Tnu);

// Get linear P(k) from input transfer
LinearPS P_l(Cm, 0.);
LinearPS P_cb(Ccb, 0.);
LinearPS P_nu(Cnu, 0.);

// Load halo class witth all linear P(k)
HALO halo(Cm, P_l, P_cb, P_nu, epsrel);
// Special functions class
IOW iow;

//SPT spt(Cm, P_l, epsrel);

// store params for passing into React functions
double vars[7];
    vars[0] = 1./(myz+1.); //  scale factor
    vars[1] = Omega_m;
    vars[2] = modg; //  modified gravity param
    vars[3] = 1.;  // extra
    vars[4] = 1.; // extra
    vars[5] = massb; // number of mass bins between 5<Log10[M]<18
    vars[6] = Omega_nu;

// initialise spherical collapse quantities and reaction quantities
halo.initialise(vars,mgcamb);

//halo.react_init(vars);
// Output section
double kmin = 1e-3;
double kmax = 0.5;



// load Matteo's data pk for testing
//ifstream fin2("validate/matteo_data/1-halo/mnu_0.4_P1h_pseudo_z0.dat");
ifstream fin2("validate/matteo_data/reactions/mnu_0.4_reaction_z1.dat");

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


double p1,p2,p3;
 //output file name
 const char* output = "mytest.dat";
 /* Open output file */
 FILE* fp = fopen(output, "w");

 for(int i =0; i < Nk;  i ++) {
      real k = mypk[i][0];
      // quantity to validate
      double react_q = halo.reaction_nu(k,vars);
      p1 = react_q/mypk[i][1];

//if(fabs(p1-1.)>0.001){
     printf("%d %e %e \n", i, k, p1); // print to terminal
     fprintf(fp,"%e %e \n", k, p1); // print to file
//}
}

	/*close output file*/
    fclose(fp);
    return 0;
}
