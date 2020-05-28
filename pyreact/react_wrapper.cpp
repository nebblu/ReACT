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

    int compute_reaction(int* N_k_pk, double* powerspectrum,
                         int* N_k, double* kvals,
                         int* N_z, double* zvals,
                         bool* is_transfer,
                         double* h, double* n_s, double* Omega_m, double* Omega_b, double* sigma_8,
                         double* mg1, int* mass_loop,
                         int* N_k_react, int* N_z_react, double* output_react,
                         int* N_k_pl, int* N_z_pl, double* output_pl,
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
            std::cout<<"mass loop: " << *mass_loop << "\n";
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
        HALO halo(C, P_l, epsrel);
        SPT spt(C, P_l, epsrel);
        IOW iow;

        /* array set up for 1-loop spt splining over z and reaction splining over k */
        int loop_nk = 60;
        double ploopr[*N_z],ploopp[*N_z],react_tab[loop_nk],k2val[loop_nk];
        /* declare splines */
        Spline mysr,mysp,myreact;
        /*declare variable array*/
        double vars[6];
        vars[0] = 1.;
        vars[1] = *Omega_m;
        vars[2] = *mg1;//pow(10.0,*mg1); // edit this for new mg param
        vars[3] = 1.;
        vars[4] = 1.;
        vars[5] = *mass_loop;
        // initialise power spectrum normalisation before running 1-loop computations
        iow.initnorm(vars);

        // initialise splines over redshift for real and pseudo 1-loop spectra @ k0 = 0.06h/Mpc
        double k0 = 0.06;
        for(int i=0; i<*N_z; i++) {
            ploopr[i] = 0.;
            ploopp[i] = 0.;
        }

        // 1-loop computations at all redshifts @ k0
        spt.ploop_init(ploopr,ploopp, zvals , *N_z, vars, k0);

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
            iow.initnorm(vars);
            // Spherical collapse stuff
            /// initialise delta_c(M), a_vir(M), delta_avir(M) and v(M)
            status = halo.scol_init(vars);
            status |= halo.scol_initp(vars);

            if(status != 0) {
                react_error("Failed to compute spherical collapse.");
                return 1;
            }

              // initialise k_star and mathcal{E}
            halo.react_init2(vars,mysr,mysp);

            // reaction
            //#pragma omp parallel for
            for(int i =0; i < loop_nk;  i ++) {
                // 1.01 is put so that spline covers full range of input k values
                k2val[i] = kvals[0] * exp(i*log(kvals[*N_k-1]*1.01/kvals[0])/(loop_nk-1.));
                if(k2val[i]<0.01){
                    react_tab[i] = 1.;
                }
                else {
                    react_tab[i] = halo.reaction(k2val[i], vars);
                }
            }

            myreact = CubicSpline(loop_nk, k2val, react_tab);

            for(int i =0; i < *N_k;  i ++) {
                output_react[i*(*N_z)+j] =  myreact(kvals[i]);
                output_pl[i*(*N_z)+j] = halo.plinear_cosmosis(kvals[i]);
                if(*verbose > 1) {
                    printf(" %e %e %e \n",zvals[j], kvals[i],halo.plinear_cosmosis(kvals[i]));
                }
            }
        }

        return 0;
    }
}
