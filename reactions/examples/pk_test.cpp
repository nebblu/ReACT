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



// Load transfer function at z from MGCAMB with all species
ifstream fin("tests/nu0/f4/test_transfer_z05.dat");

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
                  Tnu[i] = 1e10;//mytrans[i][5]/(mytrans[0][5]+100);
              }


// integration error
real epsrel = 1e-2;

// // // Specify params
double h  = 0.7;
double n_s = 0.96; // spectral index
double Omega_m =  0.3; // total matter fraction
double Omega_b  = 0.0462; //  baryon fraction
double Omega_nu = 0.000;  // neutrino fraction
double sigma_8m = 0.6699; // needs to be specified at each redshift ---- ISSUE? Or we specify pk at each z ...
double modg = 1e-4;  //  modified gravity param
double massb = 30.; // number of mass bins between 5<Log10[M]<18
// // input transfer redshift
double myz = 0.5;

bool mgcamb = true; // mgcamb transfer input or camb transfer input

// LCDM transfer at z=0 for Copter growth --- current ReACT implementation (for mgcamb=false below)
// const char* cstr = "tests/nu0/test";
// Cosmology C(cstr);
// LinearPS P_l(C, 0.);
// HALO halo(C, P_l, P_l, P_l, epsrel);

// Load cosmology classes
Cosmology Cm(h, n_s, Omega_m, Omega_b, sigma_8m, ki, Tm);
Cosmology Ccb(h, n_s, Omega_m, Omega_b, sigma_8m, ki, Tcb);
Cosmology Cnu(h, n_s, Omega_m, Omega_b, sigma_8m, ki, Tnu);

// Get linear P(k) from input transfer
LinearPS P_l(Cm, 0.);
LinearPS P_cb(Ccb, 0.);
LinearPS P_nu(Cnu, 0.);

// Load halo class witth all linear P(k)
HALO halo(Cm, P_l, P_cb, P_nu, epsrel);


// Special functions class
IOW iow;


// load MGCAMB pk for testing
ifstream fin2("testing.dat");

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
double kmin = 1e-2;
double kmax = 10.;
int Nk = 100;

double p1,p2,p3;
 //output file name
 const char* output = "testing3.dat";
 /* Open output file */
 FILE* fp = fopen(output, "w");

 for(int i =0; i < mypk.size();  i ++) {
      real k =  mypk[i][0];
      p1 =  halo.reaction_nu(k, vars)/mypk[i][2];


if(fabs(p2-1.)>0.001){
  //if(fabs(p2/p3-1.)>0.01){
     printf("%d %e %e \n", i, k, p1); // print to terminal
     fprintf(fp,"%e %e \n", k, p1); // print to file
}

}

	/*close output file*/
    fclose(fp);

// Check pure P_l input - output
// int Nk = mypk.size();
// int* Nk2 = &Nk;
//
// array kii(*Nk2);
// array Tnu(*Nk2);
//
// for(int i = 0; i< Nk; i++){
//         kii[i] = mypk[i][0];
//         Tnu[i] = mypk[i][1];
//     }
//
// Cosmology Cnu(h, n_s, Omega_nu, Omega_b, sigma_8nu, kii, Tnu);
//
// myLinearPS P_nu(Cnu, myz);
//
// for(int i = 0; i<Nk; i++){
//     double k = mypk[i][0];
//     double pkcamb = mypk[i][1];
//     double pkcop = P_nu(k);
//     printf("%e %e %e %e \n", k, pkcamb, pkcop, pkcamb/pkcop);
// }

    return 0;
}
