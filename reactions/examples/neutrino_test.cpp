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


extern "C" {


int main(int argc, char* argv[]) {

// input transfer file from MGCAMB
ifstream fin("tests/test_transfer_z0.5.dat");

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


int Nkr = mytrans.size();
int* Nkt = &Nkr;
array Tm(*Nkt);
array Tcb(*Nkt);
// array Tnu(*Nkt);
array ki(*Nkt);

//double Tm[Nkr],Tcb[Nkr],Tnu[Nkr],ki[Nkr];
// Populate arrays and normalise to 1 at large scales
          for(int i = 0; i< Nkr; i++){
                  ki[i] = mytrans[i][0];
                  Tm[i] = mytrans[i][6]/mytrans[0][6];
                  Tcb[i] = mytrans[i][7]/mytrans[0][7];
                //  Tnu[i] = mytrans[i][5]/mytrans[0][5];
              }

// Relative error in magnitude integrations
real epsrel = 1e-2;

// Specify params
double h  = 0.68;
double n_s = 0.96; // IMP
double Omega_m =  0.3072; // total matter fraction
double Omega_b  = 0.048; //  baryon fraction
double Omega_nu = 0.0006;  // neutrino fraction
double sigma_8m = 0.8971; // IMP
double sigma_8cb = 0.8989; // IMP
double sigma_8nu = 0.0375; // IMP
double modg = 1e-5;  //  modified gravity param
double massb = 30.; // number of mass bins between 5<Log10[M]<18

// desired output redshift
double myz = 0.;

Cosmology Cm(h, n_s, Omega_m, Omega_b, sigma_8m, ki, Tm);
Cosmology Ccb(h, n_s, Omega_m, Omega_b, sigma_8cb, ki, Tcb);
// Cosmology Cnu(h, n_s, Omega_nu, Omega_b, sigma_8nu, ki, Tnu);

LinearPS P_l(Cm, myz);
LinearPS P_cb(Ccb, myz);
// myLinearPS P_nu(Cnu, myz);

ifstream fin2("tests/test_nupower_z05.dat");

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
int* Nk2 = &Nk;

array kii(*Nk2);
array Tnu(*Nk2);

for(int i = 0; i< Nk; i++){
        kii[i] = mypk[i][0];
        Tnu[i] = mypk[i][1];
    }

Cosmology Cnu(h, n_s, Omega_nu, Omega_b, sigma_8nu, kii, Tnu);

myLinearPS P_nu(Cnu, myz);

for(int i = 0; i<Nk; i++){
    double k = mypk[i][0];
    double pkcamb = mypk[i][1];
    double pkcop = P_nu(k);
    printf("%e %e %e %e \n", k, pkcamb, pkcop, pkcamb/pkcop);
}


// HALO halo_mu(C, P_m, P_cb, P_nu, epsrel)
// IOW iow;

// double vars[7];
//
//     vars[0] = 1./(myz+1.); //  scale factor
//     vars[1] = Omega_m; // total matter fraction
//     vars[2] = Omega_nu; // neutrino fraction
//     vars[3] = modg; //  modified gravity param
//     vars[4] = 1.;  // extra
//     vars[5] = 1.; // extra
//     vars[6] = massb; // number of mass bins between 5<Log10[M]<18
//     vars[7] = mypar;
// initialise all  relevant growth factors
// iow.initnorm_nu(vars);
// /// initialise delta_c(M), a_vir(M), delta_avir(M) and v(M)
// halo_mu.scol_initnu(vars);
// halo_mu.scol_initpnu(vars);
// halo_mu.react_initnu(vars);

// double kmin = 1e-2;
// double kmax = 10.;
// int Nk = 100;

//output file name
//  const char* output = "ps_f5nu_z0.dat";
//  /* Open output file */
//  FILE* fp = fopen(output, "w");
//
//  for(int i =0; i < Nk;  i ++) {
//       // real k =  kmin * exp(i*log(kmax/kmin)/(Nk-1));
//       // p1 =  halo_mu.one_halonu(k, vars);
//       // p2 =  halo_mu.one_halopnu(k, vars);
//       // p3 =  halo_mu.reactionnu(k, vars);
//
//      printf("%d %e %e %e %e \n", i, k, p1,p2,p3); // print to terminal
//      fprintf(fp,"%e %e %e %e \n", k, p1,p2, p3); // print to file
//
// }

	/*close output file*/
//    fclose(fp);
    return 0;
}

}
