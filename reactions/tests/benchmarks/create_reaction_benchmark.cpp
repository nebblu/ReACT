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
vector<vector<double> > allDatamult;

int main() {


    ifstream fin("transfer.dat");
    string line;
    while (getline(fin, line)) {      // for each line
        vector<double> lineData;           // create a new row
        double val;
        istringstream lineStream(line);
        while (lineStream >> val) {          // for each value in line
            lineData.push_back(val);           // add to the current row
        }
        allDatamult.push_back(lineData);         // add row to allData
    }


    int Nk = allDatamult.size(); // transfer function number of k-values

   int Nz = 11; // number of redshifts
   double zvals[11] = {2.,1.8,1.6,1.4,1.2,1.0,0.8,0.6,0.4,0.2,0.};


 // specify the redshifts
    double cosmology[5] = {0.68,0.9645,0.308489,0.0475779,0.815088}; // specify the cosmology
    double fR0 = 1e-5;
    double transfer[Nk], kvals[Nk]; // declar transfer function array and the k-values it's calculated at
    // initialise those arrays
    int mass_loop = 30;

    for(int i=0; i< Nk ; i++){
        kvals[i] = allDatamult[i][0];
        transfer[i] = allDatamult[i][1];
    }

    bool is_tranfer = true;

    // // output array
    // double myoutput[Nk*Nz];

    // compute reactions
    // int initit = compute_reaction(&Nk, transfer, &Nk, kvals, &Nz, zvals,
    //             &is_tranfer,
    //             &cosmology[0], &cosmology[1], &cosmology[2], &cosmology[3], &cosmology[4],
    //             &fR0,
    //             &Nk, &Nz, myoutput);


    // output array
    double myoutput_react[(Nk-1)*(Nz)+Nz-1];
    double myoutput_pl[(Nk-1)*(Nz)+Nz-1];

     int initit = compute_reaction(&Nk, transfer, &Nk, kvals, &Nz, zvals, &is_tranfer,
                          &cosmology[0], &cosmology[1], &cosmology[2], &cosmology[3], &cosmology[4],
                          &fR0, &mass_loop,
                          &Nk, &Nz, myoutput_react,
                          &Nk,  &Nz,  myoutput_pl);

    // output
    //output file name
    const char* output = "reaction.dat";
    /* Open output file */
    FILE* fp = fopen(output, "w");
    for(int j =0; j< Nz; j++){
            for(int i =0; i < Nk;  i ++) {
            double b1 = myoutput_react[i*(Nz)+(Nz-1-j)];
            fprintf(fp, "%e ", b1);
            }
        fprintf(fp, "\n");
    }
    fclose(fp);

    fp = fopen("reaction_k_h.dat", "w");
    for(int i =0; i < Nk;  i ++) {
        fprintf(fp,"%e\n",kvals[i]); // print to terminal
    }
    fclose(fp);

    fp = fopen("reaction_z.dat", "w");
    for(int i =0; i < Nz;  i ++) {
        fprintf(fp,"%e\n",zvals[Nz-1-i]); // print to terminal
    }
    fclose(fp);



    return 0;
}
