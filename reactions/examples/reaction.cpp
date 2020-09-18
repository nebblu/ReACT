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


// Load transfer function at z from MGCAMB with all species at some redshift
ifstream fin("dustgrain/0.0ev/F6_transfer_out_z0.dat");
//ifstream fin("validate/matteo_data/transfers/F5_matteo_transfer_out_z1.dat");

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

ifstream finlcdm("dustgrain/0.0ev/GR_transfer_out_z0.dat");
//ifstream finlcdm("validate/matteo_data/transfers/GR_matteo_transfer_out_z1.dat");

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

// Store modified transfers
  for(int i = 0; i< Nkr; i++){
          ki[i] = mytrans[i][0];
          Tm[i] = mytrans[i][6];
          Tcb[i] = mytrans[i][7];
          Tnu[i] = mytrans[i][5]; // should be adjusted if this column is all 0s ....
      }
// store LCDM CB transfer
  for(int i = 0; i< Nkr; i++){
          kil[i] = mytransl[i][0];
          Tcbl[i] = mytransl[i][7];
      }



// integration error
real epsrel = 1e-2;

/*Specify params*/

/* Dustgrain */
double h  = 0.6731;
double n_s = 0.9658; // spectral index
double Omega_m = 0.31345; // total matter fraction
double Omega_b  = 0.0481; //  baryon fraction
double Omega_nu = 0.0;  // neutrino fraction
double pscale = 0.05;
double As = 2.199e-9;
double modg = 1e-6;  //  modified gravity param
double massb = 50.; // number of mass bins between 5<Log10[M]<18

/* Matteo */
// double h  = 0.68;
// double n_s = 0.9645; // spectral index
// double Omega_m = 0.3072; // total matter fraction
// double Omega_b  = 0.0481; //  baryon fraction
// double Omega_nu = 0.0;  // neutrino fraction
// double pscale = 0.05;
// double As = 2.085e-9;
// double modg = 1e-5;  //  modified gravity param
// double massb = 36.; // number of mass bins between 5<Log10[M]<18
//

// // input transfer redshift
double myz = 0.;

// if false, Omega_nu should be 0 --- create flag.
bool mgcamb = true; // mgcamb transfer input or camb transfer input -- maybe rename to massnu because it's more correctly massive neutrinos or not

// Load cosmology classes
Cosmology Cm(h, n_s, Omega_m, Omega_b, As, pscale, ki, Tm);
Cosmology Ccb(h, n_s, Omega_m, Omega_b, As, pscale, ki, Tcb);
Cosmology Ccbl(h, n_s, Omega_m, Omega_b, As, pscale, ki, Tcbl);
Cosmology Cnu(h, n_s, Omega_m, Omega_b, As, pscale, ki, Tnu);


// Get linear P(k) from input transfer
LinearPS_as P_l(Cm, 0.);
LinearPS_as P_cb(Ccb, 0.);
LinearPS_as P_nu(Cnu, 0.);
LinearPS_as P_cbl(Ccbl, 0.);

// Load halo class witth all linear P(k)
HALO halo(Cm, P_l, P_cb, P_nu, P_cbl, epsrel);

IOW iow;
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
halo.phinit_pseudo(vars,mgcamb);

// Output section
double kmin = 1e-3;
double kmax = 0.5;

// load in k-binning from sims
ifstream fin2("dustgrain/LCDM_power_all_033.txt");
//ifstream fin2("validate/matteo_data/reactions/HM_reaction_standard_HMF_F5_z1.dat");


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
 const char* output = "F6_0ev_z0.dat";
 /* Open output file */
 FILE* fp = fopen(output, "w");

 for(int i =0; i <Nk;  i ++) {

      real k = mypk[i][0];
       p1 = halo.reaction(k,vars);
       p2 = halo.PHALO_pseudo(k,mgcamb);
       p3 = p2*p1;

     printf("%d %e %e %e %e \n", i, k, p1,p2,p3); // print to terminal
     fprintf(fp,"%e %e %e %e \n", k, p1,p2,p3); // print to file
}

	/*close output file*/
    fclose(fp);
    return 0;
}
