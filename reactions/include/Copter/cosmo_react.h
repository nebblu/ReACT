#ifndef COSMO_REACT_NEW_H
#define COSMO_REACT_NEW_H

extern "C" {

//int test_func(int* N, int* M, double array[*N][*M]);

// //int myreaction(int Nk, int Nz, double powerspectrum[][100], int Nk2, double kvals[], int Nz2, const double zvals[], double cosmology[], double fR0, int No1, int No2, double output[][100]);
// int compute_reaction(int* N_k_pk, int* N_z_pk, double powerspectrum[*N_k_pk][*N_z_pk],
//                      int* N_k, double kvals[*N_k],
//                      int* N_z, double zvals[*N_z],
//                      double* h, double* n_s, double* Omega_m, double* Omega_b, double* sigma_8,
//                      double* fR0,
//                      int* N_k_react, int* N_z_react, double output[*N_k_react][*N_z_react]);
//
//

 int test_func(int* N, int* M, double* array);

 // int compute_reaction(int* N_k_pk, double* powerspectrum,
 //                         int* N_k, double* kvals,
 //                         int* N_z, double* zvals,
 //                         bool* is_transfer,
 //                         double* h, double* n_s, double* Omega_m, double* Omega_b, double* sigma_8,
 //                         double* fR0,
 //                         int* N_k_react, int* N_z_react, double* output);
 //

 int compute_reaction(int* N_k_pk, double* powerspectrum,
                      int* N_k, double* kvals,
                      int* N_z, double* zvals,
                      bool* is_transfer,
                      double* h, double* n_s, double* Omega_m, double* Omega_b, double* sigma_8,
                      double* fR0, int* mass_loop,
                      int* N_k_react, int* N_z_react, double* output_react,
                      int* N_k_pl, int* N_z_pl,  double* output_pl);



// specify transfer function at z=0 (transfer) at k-values given by kvalst in an array of size N_k_t
// Output at scales given in array kvalsout (size N_k_out) and at redshifts given in zvals array
int compute_reaction_new(int* N_k_t, int* N_z, double* transfer, double* kvalst, double* zvals,
                         double* h, double* n_s, double* Omega_m, double* Omega_b, double* sigma_8,
                         double* fR0,
                         int* N_k_out, double* kvalsout, double* output);
}

#endif
