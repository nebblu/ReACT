#ifndef CORRELATION_FUNCTION_H
#define CORRELATION_FUNCTION_H
#include "Common.h"
#include "array.h"
#include "SPT.h"
#include "Cosmology.h"


/* Compute the quantity
 *   \xi_l^m(r) = \int_0^\infty \frac{dk}{2\pi^2} k^m j_l(kr) P(k)
 * using Simpson's rule.  The parameters Nk, kmin, and kmax determine the
 * accuracy of the integration.  Note that Nk must be even. */
void ComputeXiLM(int l, int m, const PowerSpectrum& P,
                 int Nr, const double r[], double xi[],
                 int Nk = 32768, double kmin = 0., double kmax = 30.);


class CorrelationFunction {
public:
    CorrelationFunction(const Cosmology& C, const PowerSpectrum& P, real kmin = 1e-3, real kmax = 30.);

    real Evaluate(real r) const;
    real operator()(real r) const { return Evaluate(r); }

    array EvaluateMany(const array& r) const;
    array operator()(const array& r) const { return EvaluateMany(r); }

    const PowerSpectrum& GetPowerSpectrum() const { return P; }

    /* Gaussian Streaming Model */
    // calls all initialization functions needed for GSM
    void gsm_init(int a, int b, double s8, double scalef, double omega0, double smo, double mg1, double mg2, double mg3)const;

    //Correlation function : Linear, 1-loop and FT of TNS multipoles
    //a =1 : Linear(GR)
    //a =2 : Regpt 1-loop(GR)
    //a =3 : linear(MG)
    //a =4 : Regpt 1-loop(MG)
    //a =5 : TNS_0
    //a =6 : TNS_2
    //a =7 : TNS_4
      real xi_real(double r, int a) const ;

    // Precompute 1-loop PS in KMIN_e, KMAX_e range ( extra parameters are for scale dep numerical PS)
    // a  chooses which PS
    // b  chooses scale dep or not (b=0 scale indep, b=1 scale dep)
    void loop_init(int a, int b, double scalef, double smo, double omega0, double p1, double p2, double p3 )const;

    // Precompute integrands for A terms with bias b as functions of k
    void Aterms_init(int a, int b, double scalef, double omega0, double p1, double p2, double p3)const;
    // Precompute sig_12 terms as functions of r
    void sig_init(int a, int loop)const;

      // r is separation
      // b is galaxy bias
      // a  chooses scale dep or not (a=1 scale indep, a=2 scale dep)

      /* Mean infall velocity */
    	real v12(int a, double b, double r) const;
    	/* Mean linear infall velocity */
    	real v12_L(int a, double b, double r) const;

      // These terms require spline initialization done by sig_init (see above)
      // sig_init chooses numerical or analytic forms
      // u is the square of the cosine between separation vector and LINE OF SIGHT
    /*Mean velocity dispersion */
    	real sig12(double b, double r, double u) const;
    /*Mean linear velocity dispersion */
    	real sig12_L(double r, double u) const;

      // gsm integral functions
    real gsm_integrand_inner(int a, double b, double r_sigma, double r_pi, double y) const;
    void gsm_xi_init(int a, double b,  double s, double mu_s)const;

    // 2D GSM
    // s is separation and mu_s is cosine of angle between LOS and s(vector)
    // a  chooses scale dep or not (a=1 scale indep, a=2 scale dep)
    // b is galaxy bias
    real gsm_xi(int a, double b, double s, double  mu_s ) const;
    // multipoles of GSM
    // order : 0: monopole, 2: quadrupole, 4: hexdecapole
    real gsm_multi(int a, int order, double b, double s ) const;

    // Same for LSM
    real lsm_integrand_inner(double b, double r_sigma, double r_pi, double y) const;
    void lsm_xi_init(double b,  double s, double mu_s)const;
    real lsm_xi(double b, double s, double  mu_s ) const;
    real lsm_multi(int order, double b, double s ) const;


protected:
    const PowerSpectrum& P;
    const Cosmology& C;
    real kmin, kmax;
};

#endif // CORRELATION_FUNCTION_H
