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
#include <Copter/cosmo_react.h>

using namespace std;
vector<vector<double> > allDatamultz,allDatamultk,allDatamultpofk;


/* Code to test the reaction libraries */

int main() {


    ifstream fin1("z.txt");
    string line1;
    while (getline(fin1, line1)) {      // for each line
            vector<double> lineData;           // create a new row
            double val;
            istringstream lineStream(line1);
            while (lineStream >> val) {          // for each value in line
                    lineData.push_back(val);           // add to the current row
            }
            allDatamultz.push_back(lineData);         // add row to allData
    }

    ifstream fin2("k.txt");
    string line2;
    while (getline(fin2, line2)) {      // for each line
            vector<double> lineData;           // create a new row
            double val;
            istringstream lineStream(line2);
            while (lineStream >> val) {          // for each value in line
                    lineData.push_back(val);           // add to the current row
            }
            allDatamultk.push_back(lineData);         // add row to allData
    }


    ifstream fin3("pofk.txt");
    string line3;
    while (getline(fin3, line3)) {      // for each line
            vector<double> lineData;           // create a new row
            double val;
            istringstream lineStream(line3);
            while (lineStream >> val) {          // for each value in line
                    lineData.push_back(val);           // add to the current row
            }
            allDatamultpofk.push_back(lineData);         // add row to allData
    }


	 //output file name
    const char* output = "check.dat";

    /* Open output file */
    FILE* fp = fopen(output, "w");


 int Nk = allDatamultk.size();
 int Nz = allDatamultz.size();


double fr0 = -5;
double cosmology[5] = {0.68,0.9645,0.308489,0.0475779,0.815088}; // specify the cosmology
int mass_loop =30;



double kvals[Nk],pofk[Nk],zvals[Nz];

for(int i=0; i< Nz ; i++){
  zvals[i] = allDatamultz[Nz-i-1][0];
}

for(int i=0; i< Nk ; i++){
  kvals[i] = allDatamultk[i][0];
  pofk[i] =  allDatamultpofk[i][0];
}


// output array
double myoutput_react[(Nk-1)*(Nz)+Nz-1];
double myoutput_pl[(Nk-1)*(Nz)+Nz-1];

bool myb = 0;
// compute reactions
int initit = compute_reaction(&Nk, pofk, &Nk, kvals, &Nz, zvals, &myb,
                        &cosmology[0], &cosmology[1], &cosmology[2], &cosmology[3], &cosmology[4],
                        &fr0, &mass_loop,
                        &Nk, &Nz, myoutput_react,
                        &Nk,  &Nz,  myoutput_pl);


// output
for(int j =0; j< Nz; j++){
  double myz = zvals[j];
for(int i =0; i < Nk;  i ++) {
  double k = kvals[i];
  double b1 = myoutput_react[i*(Nz)+j];
  double b2 = myoutput_pl[i*(Nz)+j];

fprintf(fp,"%e %e %e %e\n",myz,k,b1,b2); // print to terminal
}
}

	/*close output file*/
    fclose(fp);
    return 0;
}
