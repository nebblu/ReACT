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

extern "C" {
    extern const int ERROR_MESSAGE_LEN = 512;
    char error_message[ERROR_MESSAGE_LEN];

    void react_error_numerical(double val1 = 0, double val2  = 0,double  val3 = 0, double  val4 = 0, double val5 = 0) {
        const int truncate = ERROR_MESSAGE_LEN-10;

        #pragma omp critical
        {
            snprintf(error_message, truncate, "Reaction error: %e %e %e %e %e ", val1, val2, val3, val4, val5);
            fprintf(stderr, " %e %e %e %e %e \n", val1, val2,val3,val4,val5);
        }
    }

    void react_error(const char* msg) {
        const int truncate = ERROR_MESSAGE_LEN-10;

        #pragma omp critical
        {
            snprintf(error_message, truncate, "Reaction error: %s", msg);
            fprintf(stderr, "%s\n", error_message);
        }
    }

    int test_func(int* N, int* M, double* array)
    {
        std::cout << "N: " << *N << " M: " << *M << "\n";
        for(int i=0; i<*N; i++)
        {
            for(int j=0; j<*M; j++)
            {
                std::cout<<array[i*(*M)+j] << " ";
            }
            std::cout<<"\n";
        }
        return 42;
    }


    /* Option 1: Run wrapper as normally, i.e. if input transfer/pofk is LCDM at z=0 */
    int compute_reaction(int* N_k_pk, double* powerspectrum,
                         int* N_k, double* kvals,
                         int* N_z, double* zvals,
                         bool* is_transfer,
                         double* h, double* n_s, double* Omega_m, double* Omega_b, double* sigma_8,
                         double* mg1, double* mg2, double* mg3, int* mass_loop, int* model,
                         int* N_k_react, int* N_z_react, double* output_react,
                         int* N_k_pl, int* N_z_pl, double* output_pl,
                         double* modsig8,
                         int* verbose)
    {
        int status = 0;
        if(*N_k != *N_k_pk || *N_k != *N_k_react || *N_k != *N_k_pl){
            react_error("Inconsistency in k array sizes");
            return 1;
        }
        if(*N_z != *N_z_react || *N_z != *N_z_pl){
            react_error("Inconsistency in z array sizes");
            return 1;
        }

        if(kvals[0] > 1e-4) {
            react_error("k must at least go down to 1e-4");
            return 1;
        }
        if(zvals[0] < zvals[*N_z-1]) {
            react_error("z must be in decreasing order.");
            return 1;
        }
        if(zvals[0] > 2.5) {
            react_error("z must be smaller than 2.5");
            return 1;
        }


        if(*verbose > 0) {
            std::cout<<"z: " << zvals[0] << " -> " << zvals[*N_z-1] << "\n";
            std::cout<<"k: " << kvals[0] << " -> " << kvals[*N_k-1] << "\n";
            std::cout<<"h: " << *h << "\n";
            std::cout<<"n_s: " << *n_s << "\n";
            std::cout<<"Omega_m: " << *Omega_m << "\n";
            std::cout<<"Omega_b: " << *Omega_b << "\n";
            std::cout<<"sigma_8: " << *sigma_8 << "\n";
            std::cout<<"mg1: " << *mg1 << "\n";
            std::cout<<"mg2: " << *mg2 << "\n";
            std::cout<<"mg3: " << *mg3 << "\n";
            std::cout<<"mass loop: " << *mass_loop << "\n";
            std::cout<<"model: "  << *model << " (1:GR, 2:f(R), 3:DGP, 4:quintessence, 5: CPL)\n";

        }


        array Ti(*N_k);
        array ki(*N_k, kvals);

        double T0;
        if(!(*is_transfer)) {
            T0 = sqrt(powerspectrum[0] / pow(ki[0], *n_s));
        }

        for(int i = 0; i<*N_k; i++){
            if(*is_transfer) {
                Ti[i] = powerspectrum[i];
            }
            else {
                // Copter takes transfer function as input - Easiest fix is to just convert input PS to T
                Ti[i] = sqrt(powerspectrum[i] / pow(ki[i], *n_s));
                // Must normalise T to 1. at large scales.
                Ti[i] /= T0;
            }
        }

        if(*verbose > 0) {
            std::cout<<"T(k): " << Ti[0] << ", " << Ti[1] << " -> " << Ti[*N_k-1] << "\n";
            std::cout<<"is_transfer: " << *is_transfer << "\n";
        }


        // Relative error in magnitude integrations
        real epsrel = 1e-2;

        Cosmology C(*h, *n_s, *Omega_m, *Omega_b, *sigma_8, ki, Ti);
        LinearPS P_l(C, 0.0);
        HALO halo(C, P_l, P_l, P_l, P_l, epsrel);
        SPT spt(C, P_l, epsrel);
        IOW iow;

        /* array set up for 1-loop spt splining over z and reaction splining over k */
        int loop_nk = 60;
        double ploopr[*N_z],ploopp[*N_z],react_tab[loop_nk],k2val[loop_nk];
        /* declare splines */
        Spline mysr,mysp,myreact;
        /*declare variable array*/
        double vars[8];
        vars[0] = 1.;
        vars[1] = *Omega_m;
        vars[2] = *mg1;
        vars[3] = *mg2;
        vars[4] = *mg3;
        vars[5] = *mass_loop;
        vars[6] = 0.0; // massive neutrino mass sum, not implemented in pyreact yet!

        bool modcamb = false;
        int mod = *model;
        // initialise power spectrum normalisation before running 1-loop computations
        iow.initnorm(vars,mod);

        // initialise splines over redshift for real and pseudo 1-loop spectra @ k0 = 0.06h/Mpc
        double k0 = 0.06;
        for(int i=0; i<*N_z; i++) {
            ploopr[i] = 0.;
            ploopp[i] = 0.;
        }

        // Switch on 1-loop modified gravity correction if DGP or f(R)
        bool modg;
        if (mod == 2 || mod == 3) {
          // 1-loop computations at all redshifts @ k0
          spt.ploop_init(ploopr,ploopp, zvals , *N_z, vars, mod, k0);
          modg = true;
        }
        else{
          modg = false;
        }

        double myscalef[*N_z];
        for(int i = 0; i<*N_z ; i++){
            myscalef[i]= 1./(1.+zvals[i]);
        }

        mysr = CubicSpline(*N_z, myscalef, ploopr);
        mysp = CubicSpline(*N_z, myscalef, ploopp);


        // main loop
        for(int j = 0; j<*N_z; j++) {
            vars[0] = 1./(1.+zvals[j]);

            // initialise wCDM/LCDM lin growth for PS normalisation
            iow.initnorm(vars,mod);
            // Spherical collapse stuff
            /// initialise delta_c(M), a_vir(M), delta_avir(M) and v(M)
            status = halo.scol_init(vars,modcamb,mod);
            status |= halo.scol_initp(vars,modcamb);

            // store modified sigma 8 at z=0
            if(j == *N_z-1) {
                modsig8[0] = vars[7];
            }

            if(status != 0) {
                react_error("Failed to compute spherical collapse.");
                return 1;
            }

            // initialise k_star and mathcal{E}
            halo.react_init_nu_multiz(vars,mysr,mysp,modcamb,modg);  // kstar and mathcal E for reaction for multiple redshifts

            // reaction
            //#pragma omp parallel for
            for(int i =0; i < loop_nk;  i ++) {
                // 1.01 is put so that spline covers full range of input k values
                k2val[i] = kvals[0] * exp(i*log(kvals[*N_k-1]*1.01/kvals[0])/(loop_nk-1.));
                if(k2val[i]<0.01){
                    react_tab[i] = 1.;
                }
                else {
                    react_tab[i] = halo.reaction_nu(k2val[i], vars, modcamb);
                }
            }

            myreact = CubicSpline(loop_nk, k2val, react_tab);

            for(int i =0; i < *N_k;  i ++) {
                output_react[i*(*N_z)+j] =  myreact(kvals[i]);
                output_pl[i*(*N_z)+j] = pow2(halo.Lin_Grow(kvals[i]))*P_l(kvals[i]);
                if(*verbose > 1) {
                    printf(" %e %e %e \n",zvals[j], kvals[i],halo.plinear_cosmosis(kvals[i]));
                }
            }
        }

        return 0;
    }


    /* Option 2: Run wrapper assuming input is transfer/pofk array in k and z for all species in target cosmology */
    int compute_reaction_nu(int* N_pk_m, double* torpk_m,
                         int* N_k, double* kvals,
                         int* N_z, double* zvals,
                         int* N_pk_cb, double* torpk_cb,
                         int* N_pk_lcdm, double* torpk_lcdm,
                         int* N_k_lcdm, double* kvals_lcdm,
                         bool* is_transfer,
                         double* h, double* n_s, double* Omega_m, double* Omega_b, double* Omega_nu,
                         double* As, double* pscale,
                         double* mg1, double* mg2, double* mg3, int* mass_loop, int* model,
                         int* N_k_react, int* N_z_react, double* output_react,
                         int* N_k_pl, int* N_z_pl, double* output_pl,
                         double* modsig8,
                         int* verbose)
    {

        int status = 0;

        if(*N_k != *N_k_react || *N_k != *N_k_pl || *N_pk_m/(*N_z) != *N_k || *N_pk_cb/(*N_z) != *N_k ){
            react_error("Inconsistency in k array sizes");
            return 1;
        }
        if(*N_z != *N_z_react || *N_z != *N_z_pl || *N_z != *N_pk_m/(*N_k)){
            react_error("Inconsistency in z array sizes");
            return 1;
        }

        if(kvals[0] > 1e-4) {
            react_error("k must at least go down to 1e-4");
            return 1;
        }
        if(zvals[0] < zvals[*N_z-1]) {
            react_error("z must be in decreasing order.");
            return 1;
        }
        if(zvals[0] > 2.5) {
            react_error("z must be smaller than 2.5");
            return 1;
        }


        if(*verbose > 0) {
            std::cout<<"z: " << zvals[0] << " -> " << zvals[*N_z-1] << "\n";
            std::cout<<"k: " << kvals[0] << " -> " << kvals[*N_k-1] << "\n";
            std::cout<<"h: " << *h << "\n";
            std::cout<<"n_s: " << *n_s << "\n";
            std::cout<<"Omega_m: " << *Omega_m << "\n";
            std::cout<<"Omega_b: " << *Omega_b << "\n";
            std::cout<<"Omega_nu: " << *Omega_nu << "\n";
            std::cout<<"As: " << *As << "\n";
            std::cout<<"pivot scale: " << *pscale << "\n";
            std::cout<<"mg1: " << *mg1 << "\n";
            std::cout<<"mg2: " << *mg2 << "\n";
            std::cout<<"mg3: " << *mg3 << "\n";
            std::cout<<"mass loop: " << *mass_loop << "\n";
            std::cout<<"model: "  << *model << " (1:GR, 2:f(R), 3:DGP, 4:quintessence, 5: CPL)\n";
        }


/* Switch on 1-loop modified gravity correction if DGP or f(R) */
        bool modg;
        /* Which model ?*/
        int mod = *model;

        if (mod == 2 || mod == 3) {
          modg = true;
        }
        else{
          modg = false;
        }

        // Relative error in magnitude integrations
        real epsrel = 1e-2;

/* Initialise arrays for various species */
        array Tm(*N_k);
        array Tcb(*N_k);
        array Tnu(*N_k);
        array ki(*N_k, kvals);

        /* Perhaps we have a different binning for the LCDM run?*/
        array Tlcdm(*N_k_lcdm);
        array kilcdm(*N_k_lcdm, kvals_lcdm);

        /*declare variable array*/
        double vars[8];
        vars[0] = 1.;
        vars[1] = *Omega_m;
        vars[2] = *mg1;
        vars[3] = *mg2;
        vars[4] = *mg3;
        vars[5] = *mass_loop;
        vars[6] = *Omega_nu;

        bool modcamb = true;


        /* Perform some checks */
        if ( *pscale<1e-3 || *pscale>0.1 || *As<1.5e-9 || *As>2.5e-9 ) {
          react_error("Values of As and pivot scale are not compatible -  set 1e-3<pscale<0.1  and 1.5e-9<As<2.5e-9");
        }

        if (!*is_transfer) {
          react_error("You must supply the input transfer functions if using a modified linear input");
        }


        // Special function class
        IOW iow;
        // initialise power spectrum normalisation before running
        iow.initnorm(vars,mod);



/* Option 2: Run wrapper assuming input is transfer/pofk array in k and z for all species in target cosmology */
/* If transfer or pofk inputs are at the redshifts required and already modified we need to initialise classes for each redshift */


          /* PART 1: Initialise vector of classes for each z*/

          // Create vector of cosmology structure
          vector<Cosmology> mycosmo;
          mycosmo.reserve((*N_z)*4);
          // Create vector of classes
          vector<LinearPS_as> mylps;
          mylps.reserve((*N_z)*4);

          vector<HALO> myhalo;
          myhalo.reserve(*N_z);

          vector<SPT> myspt;
          myspt.reserve(*N_z);


          // Load classes for each redshift
          for(int i=0; i<*N_z; i++){

            // Load transfers for real cosmology
              for(int j = 0; j<*N_k; j++){
                  if(*is_transfer) {
                      Tm[j] = torpk_m[(*N_z-i-1)*(*N_k)+j];
                      Tcb[j] = torpk_cb[(*N_z-i-1)*(*N_k)+j];
                      if (*Omega_nu<1e-10) {
                        Tnu[j] = Tm[j];
                      }
                      else{
                      Tnu[j] = (*Omega_m*Tm[j] - (*Omega_m-*Omega_nu)*Tcb[j])/(*Omega_nu);
                      }
                  }
                  else {
                      // Copter takes transfer function as input - Easiest fix is to just convert input PS to T
                      Tm[j] = sqrt(torpk_m[(*N_z-i-1)*(*N_k)+j] / (2*M_PI*M_PI * (*As) * ki[j] * pow(*h,4) * pow(ki[j]*(*h)/(*pscale), *n_s-1.)));
                      Tcb[j] = sqrt(torpk_cb[(*N_z-i-1)*(*N_k)+j] / (2*M_PI*M_PI * (*As) * ki[j] * pow(*h,4) * pow(ki[j]*(*h)/(*pscale), *n_s-1.)));
                      if (*Omega_nu<1e-10) {
                        Tnu[j] = Tm[j];
                      }
                      else{
                      Tnu[j] = (*Omega_m*Tm[j] - (*Omega_m-*Omega_nu)*Tcb[j])/(*Omega_nu);
                        }
                    }
                  }

            // same for LCDM
                  for(int j = 0; j<*N_k_lcdm; j++){
                      if(*is_transfer) {
                        Tlcdm[j] = torpk_lcdm[(*N_z-i-1)*(*N_k_lcdm)+j];
                      }
                      else {
                          // Copter takes transfer function as input - Easiest fix is to just convert input PS to T
                        Tlcdm[j] = sqrt(torpk_lcdm[(*N_z-i-1)*(*N_k_lcdm)+j] / (2*M_PI*M_PI * (*As) * ki[j] * pow(*h,4) * pow(ki[j]*(*h)/(*pscale), *n_s-1.)));
                        }
                      }

                  // Load cosmology classes
                  mycosmo.push_back(Cosmology(*h, *n_s, *Omega_m, *Omega_b, *As, *pscale, ki, Tm));
                  mycosmo.push_back(Cosmology(*h, *n_s, *Omega_m, *Omega_b, *As, *pscale, ki, Tcb));
                  mycosmo.push_back(Cosmology(*h, *n_s, *Omega_m, *Omega_b, *As, *pscale, ki, Tnu));
                  mycosmo.push_back(Cosmology(*h, *n_s, *Omega_m, *Omega_b, *As, *pscale, kilcdm, Tlcdm));

                  mylps.push_back(LinearPS_as(mycosmo[4*i], 0.));
                  mylps.push_back(LinearPS_as(mycosmo[4*i+1], 0.));
                  mylps.push_back(LinearPS_as(mycosmo[4*i+2], 0.));
                  mylps.push_back(LinearPS_as(mycosmo[4*i+3], 0.));

                  // Load halo class with all linear P(k)
                  myhalo.push_back(HALO(mycosmo[4*i], mylps[4*i] , mylps[4*i+1], mylps[4*i+2], mylps[4*i+3], epsrel));

                  // Load Standard Perturbation Theory class
                  myspt.push_back(SPT(mycosmo[4*i+1], mylps[4*i+1], epsrel));
          }


        /* PART 2:  1-loop Perturbation theory part to get k_star */

        /* array set up for 1-loop spt splining over z */
        double ploopr[*N_z],ploopp[*N_z];

        // initialise splines over redshift for real and pseudo 1-loop spectra @ k0 = 0.06h/Mpc
        double k0 = 0.06;
        for(int i=0; i<*N_z; i++) {
            ploopr[i] = 0.;
            ploopp[i] = 0.;
        }


      if (modg) {
        /* Initialise arrays for P_L(k0)  for 1-loop real and pseudo computations */
        double pkz0[*N_z], pkz0p[*N_z];
        /* Similiarly for the integrated P_L: P_L (|k-p|) and P_L(p). We integrate with 50 absiccae in cos(theta) = k.p/(|k||p|) and 1600 points in |p| - See src/SPT.cpp for more details */
        /* Do not change number of abscissae - this is assumed in SPT.cpp's ploop_init_nu function */
        double pkz1[*N_z][50], pkz2[*N_z][1600], pkz1p[*N_z][50], pkz2p[*N_z][1600];

        // Populate power spectrum arrays
        double KMAX = QMAXp/k0; // max r = p/k (QMAXp = 30 - see SpecialFunctions.h)
        double KMIN = QMINp/k0; // min r = p/k (QMINp = 1e-4 - ...)

        double k2val,xv,karg;

        for(int i=0; i<*N_z; i++){
            pkz0[i] = mylps[i*4+1](k0);
            pkz0p[i] = mylps[i*4](k0);
          for(int j=0; j<50; j++){
              k2val = k0*KMIN * exp(j*log(KMAX/KMIN)/(50.-1.));
              pkz1[i][j] = mylps[i*4+1](k2val);
              pkz1p[i][j] = mylps[i*4](k2val);
                for(int jj=0; jj<32; jj++){
                  xv = x32[jj];
                  karg = sqrt(k2val*k2val+k0*k0-2.*k2val*k0*xv);
                  pkz2[i][jj*50+j] = mylps[i*4+1](karg);
                  pkz2p[i][jj*50+j] = mylps[i*4](karg);
            }
          }
        }

        /* Initialise 1-loop arrays as function of z @ k0 */
        myspt[0].ploop_init_nu(pkz0,pkz1,pkz2,pkz0p,pkz1p,pkz2p,ploopr,ploopp, zvals, *N_z, vars, mod, k0);

      }

        double myscalef[*N_z];
        for(int i = 0; i<*N_z ; i++){
            myscalef[i]= 1./(1.+zvals[i]);
        }

        /* Spline resulting arrays as functions of scale factor */
        Spline mysr = CubicSpline(*N_z, myscalef, ploopr);
        Spline mysp = CubicSpline(*N_z, myscalef, ploopp);


        /* PART 3: main loop to initialise reaction */
        int loop_nk = 60;
        double react_tab[loop_nk],k2val[loop_nk];

        for(int j = 0; j<*N_z; j++) {
            vars[0] = 1./(1.+zvals[j]);

            // initialise wCDM/LCDM lin growth for PS normalisation
            iow.initnorm(vars,mod);


            // Spherical collapse stuff
            /// initialise delta_c(M), a_vir(M), delta_avir(M) and v(M)
            status = myhalo[j].scol_init(vars,modcamb,mod); // real spherical collapse quantities
            status |= myhalo[j].scol_initp(vars,modcamb); // pseudo spherical collapse quantities


            // store modified sigma 8 at z=0
            if(j == *N_z-1) {
                modsig8[0] = vars[7];
            }

            // initialise k_star and mathcal{E}
            myhalo[j].react_init_nu_multiz(vars,mysr,mysp,modcamb,modg);


            // reaction
            //#pragma omp parallel for
            for(int i =0; i < loop_nk;  i ++) {
                // 1.01 is put so that spline covers full range of input k values
                k2val[i] = kvals[0] * exp(i*log(kvals[*N_k-1]*1.01/kvals[0])/(loop_nk-1.));
                if(k2val[i]<0.01){
                    react_tab[i] = 1.;
                }
                else {
                    react_tab[i] = myhalo[j].reaction_nu(k2val[i], vars, modcamb);
                }
            }


            Spline myreact = CubicSpline(loop_nk, k2val, react_tab);

        /* PART 4: Store results to output arrays */
            for(int i =0; i < *N_k;  i ++) {
                output_react[i*(*N_z)+j] =  myreact(kvals[i]);
                output_pl[i*(*N_z)+j] = mylps[4*j](kvals[i]);
                if(*verbose > 1) {
                    printf(" %e %e %e \n",zvals[j], kvals[i], mylps[4*j](kvals[i]));
                }
            }
        }

        return 0;
  }
}
