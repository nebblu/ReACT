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
vector<vector<double> > mytransl;
vector<vector<double> > mypk;


/* Example code to output the reaction and halo spectra for mg + neutrinos */
int main(int argc, char* argv[]) {

// 1: GR  2: f(R) 3: DGP 4: quintessence 5: CPL
int mymodel = 5;

// target redshift
double myz = 1.;
double Omega_nu = 0.0053;  // neutrino fraction mv = 0.24
// CPL parameters
double w0 = -0.9;
double wa = 0.1;


// Modified gravity active?
bool modg = false;

// Is the transfer being fed to ReACT of the target cosmology? If false, the transfer should be LCDM at z=0.
bool mgcamb = true;

// Load transfer function at z from MGCAMB with all species at some redshift
ifstream fin("transfers/wcdm_nu024_transfer_1.dat");

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


ifstream finlcdm("transfers/lcdm_nu024_transfer_1.dat");
// Load in the transfer data
string linelcdm;
    while (getline(finlcdm, linelcdm)) {      // for each line
            vector<double> lineData;           // create a new row
            double val;
            istringstream lineStream(linelcdm);
            while (lineStream >> val) {          // for each value in line
                    lineData.push_back(val);           // add to the current row
            }
            mytransl.push_back(lineData);         // add row to allData
    }

// Populate arrays and normalise to 1 at large scales for Cosmology class input
int Nkr = mytrans.size();
int Nkrl = mytransl.size();

int* Nkt = &Nkr;
int* Nktl = &Nkrl;


array Tm(*Nkt);
array Tcb(*Nkt);
array Tcbl(*Nktl);
array Tnu(*Nkt);
array ki(*Nkt);
array kil(*Nktl);

// integration error
real epsrel = 1e-3;

/*Specify params*/
double h  = 0.7; // hubble constant
double n_s = 0.972; // spectral index
double Omega_m = 0.2793; // total matter fraction
double Omega_b  = 0.0463; //  baryon fraction

double pscale = 0.002;
double As = 2.392e-09;

// number of mass bins between 5<Log10[M]<20
double massb = 50.;

// store params for passing into React functions
double vars[8];
    vars[0] = 1./(myz+1.); //  scale factor
    vars[1] = Omega_m;
    vars[2] = w0; //  modified gravity param or w0 (see SpecialFunctions.cpp)
    vars[3] = wa;  // extra, in the CPL case it is wa
    vars[4] = 1.; // extra
    vars[5] = massb; // number of mass bins between 5<Log10[M]<18
    vars[6] = Omega_nu;


IOW iow;

// Store modified transfers
  for(int i = 0; i< Nkr; i++){
          ki[i] = mytrans[i][0];
           Tm[i] =  mytrans[i][6]; // total
           Tcb[i] = mytrans[i][7]; // CB
           Tnu[i] = mytrans[i][5]; // neutrinos
      }
// store LCDM CB transfer
  for(int i = 0; i< Nkrl; i++){
          kil[i] = mytransl[i][0];
          Tcbl[i] = mytransl[i][7];
      }


// Load cosmology classes
Cosmology Cm(h, n_s, Omega_m, Omega_b, As, pscale, ki, Tm);
Cosmology Ccb(h, n_s, Omega_m, Omega_b, As, pscale, ki, Tcb);
Cosmology Cnu(h, n_s, Omega_m, Omega_b, As, pscale, ki, Tnu);
Cosmology Ccbl(h, n_s, Omega_m, Omega_b, As, pscale, kil, Tcbl);


// Get linear P(k) from input transfer
LinearPS_as P_l(Cm, 0.);
LinearPS_as P_cb(Ccb, 0.);
LinearPS_as P_nu(Cnu, 0.);
LinearPS_as P_cbl(Ccbl, 0.);


// Load halo class witth all linear P(k)
HALO halo(Cm, P_l, P_cb, P_nu, P_cbl, epsrel);
SPT spt(Cm, P_cb, epsrel);


//initialise spherical collapse quantities and reaction quantities
halo.initialise(vars,mgcamb,modg,mymodel);
halo.phinit_pseudo(vars,mgcamb);


ifstream fin2("benchmark_data/nu24_wcdm_z1.dat");

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


int Nk = mypk.size();
double p1,p2,p3;

 for(int i =0; i <Nk;  i ++) {

      real k = mypk[i][0];

      p1 = P_l(k)/mypk[i][1];
      p2 = halo.reaction_nu(k,vars,mgcamb)/mypk[i][2]; ;
      p3 = halo.PHALO_pseudo(k,mgcamb)/mypk[i][3];

      if(fabs(p1-1.)>0.01 || fabs(p2-1.)>0.01 || fabs(p3-1.)>0.01){
      printf("%e %e %e %e  \n", k, p1,p2, p3);
      std::cout << "Test failed: check output";
      return 0.;
      }
}

  std::cout << "Test passed";

    return 0;
}
