
#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/bind.hpp>
using boost::cref;

#include "Cosmology.h"
#include "PowerSpectrum.h"
#include "Quadrature.h"
#include "CDE.h"
#include "BSPTN.h"
#include "SPT.h"
#include "Spline.h"
#include "SpecialFunctions.h"
#include "LinearPS.h"


#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include <math.h>       /* pow */

/////////////////CLUSTERING DARK ENERGY (AND MG RESUMMED SPECTRUM) LIBRARY ///////////////////
// Load Linear PS
CDE::CDE(const Cosmology& C, const PowerSpectrum& P_l, double epsrel_)
: C(C), P_l(P_l)
{
    epsrel = epsrel;
}


/* Tree level cross bispectrum for  isosceles or equilateral configurations */
// vars specifies omega_m, cs^2, w , p3 and scale factor
double CDE::Bcde(int a, double vars[], double k, double p, double x) const {
IOW iow;
iow.initn_cde_btree(vars[4], k, p, x, vars[0], vars[1], vars[2],vars[3]);
double d =  sqrt(pow2(k) + pow2(p) + 2*p*k*x);
double omega_m = vars[0]/pow3(vars[4]);
double A = -3.*(1.+vars[2]);
double omegaf = pow(vars[4],A);
double omega_q= (1.-vars[0])*omegaf;
double cs2 = vars[1];
double norm = 2./pow4(dnorm_spt);
double terms[3];
if(d < 1e-6)
    return 0;
else
    switch(a) {
      // theta_m density_m density_m
        case 1:
						return    norm*(P_l(k)*P_l(p)*G1_nk*F1p_nk[0]*F2A_nk[0] +  P_l(d)*P_l(k)*G1_nk*F2D_nk[0]*F1kmp_nk[0]
                          + P_l(p)*P_l(d)*G2_nk[0]*F1p_nk[0]*F1kmp_nk[0]);
            break;
      // theta_q density_m density_m
        case 2:
            return    norm*(P_l(k)*P_l(p)*G1q_nk*F1p_nk[0]*F2A_nk[0] +  P_l(d)*P_l(k)*G1q_nk*F2D_nk[0]*F1kmp_nk[0]
                      + P_l(p)*P_l(d)*G2q_nk[0]*F1p_nk[0]*F1kmp_nk[0]);
            break;
      // theta_m \delta_q \delta_q
        case 3:
            return    norm*(P_l(k)*P_l(p)*G1_nk*F1pq_nk[0]*F2Aq_nk[0] +  P_l(d)*P_l(k)*G1_nk*F2Dq_nk[0]*F1kmpq_nk[0]
                      + P_l(p)*P_l(d)*G2_nk[0]*F1pq_nk[0]*F1kmpq_nk[0]);
            break;
       // theta_q \delta_q \delta_q
        case 4:
            return    norm*(P_l(k)*P_l(p)*G1q_nk*F1p_nk[0]*F2A_nk[0] +  P_l(d)*P_l(k)*G1q_nk*F2D_nk[0]*F1kmp_nk[0]
                      + P_l(p)*P_l(d)*G2q_nk[0]*F1p_nk[0]*F1kmp_nk[0]);
           break;
      // theta_m \delta_m \delta_q
        case 5:
            return    norm*(P_l(k)*P_l(p)*G1_nk*F1p_nk[0]*F2Aq_nk[0] +  P_l(d)*P_l(k)*G1_nk*F2D_nk[0]*F1kmpq_nk[0]
                      + P_l(p)*P_l(d)*G2_nk[0]*F1p_nk[0]*F1kmpq_nk[0]);
            break;
     // theta_q \delta_m \delta_q
        case 6:
            return    norm*(P_l(k)*P_l(p)*G1q_nk*F1p_nk[0]*F2Aq_nk[0] +  P_l(d)*P_l(k)*G1q_nk*F2D_nk[0]*F1kmpq_nk[0]
                      + P_l(p)*P_l(d)*G2q_nk[0]*F1p_nk[0]*F1kmpq_nk[0]);
           break;
     // theta_m \delta_t \delta_t [\delta_t = \omega_m*(\delta_m + \omega_q/\omega_m*(1+cs^2)*\delta_q) ]
        case 7:
        // theta_m density_m density_m
  				 terms[0]= norm*(P_l(k)*P_l(p)*G1_nk*F1p_nk[0]*F2A_nk[0] +  P_l(d)*P_l(k)*G1_nk*F2D_nk[0]*F1kmp_nk[0]
                         + P_l(p)*P_l(d)*G2_nk[0]*F1p_nk[0]*F1kmp_nk[0]);

        // theta_m \delta_q \delta_q
           terms[1]= norm*(P_l(k)*P_l(p)*G1_nk*F1pq_nk[0]*F2Aq_nk[0] +  P_l(d)*P_l(k)*G1_nk*F2Dq_nk[0]*F1kmpq_nk[0]
                     + P_l(p)*P_l(d)*G2_nk[0]*F1pq_nk[0]*F1kmpq_nk[0]);

        // theta_m \delta_m \delta_q
           terms[2]= norm*(P_l(k)*P_l(p)*G1_nk*F1p_nk[0]*F2Aq_nk[0] +  P_l(d)*P_l(k)*G1_nk*F2D_nk[0]*F1kmpq_nk[0]
                     + P_l(p)*P_l(d)*G2_nk[0]*F1p_nk[0]*F1kmpq_nk[0]);

          return pow2(omega_m)*(terms[0] + 2.*omega_q/omega_m*(1.+cs2)*terms[2] + pow2(omega_q/omega_m*(1.+cs2))*terms[1]);
           break;

           // theta_q \delta_t \delta_t [\delta_t = \omega_m*(\delta_m + \omega_q/\omega_m*(1+cs^2)*\delta_q) ]
       case 8:
          // theta_q density_m density_m
            terms[0]= norm*(P_l(k)*P_l(p)*G1q_nk*F1p_nk[0]*F2A_nk[0] +  P_l(d)*P_l(k)*G1q_nk*F2D_nk[0]*F1kmp_nk[0]
                      + P_l(p)*P_l(d)*G2q_nk[0]*F1p_nk[0]*F1kmp_nk[0]);

          // theta_q \delta_m \delta_q
            terms[1]=  norm*(P_l(k)*P_l(p)*G1q_nk*F1p_nk[0]*F2Aq_nk[0] +  P_l(d)*P_l(k)*G1q_nk*F2D_nk[0]*F1kmpq_nk[0]
                      + P_l(p)*P_l(d)*G2q_nk[0]*F1p_nk[0]*F1kmpq_nk[0]);

          // theta_q \delta_q \delta_q
            terms[2]=  norm*(P_l(k)*P_l(p)*G1q_nk*F1p_nk[0]*F2A_nk[0] +  P_l(d)*P_l(k)*G1q_nk*F2D_nk[0]*F1kmp_nk[0]
                      + P_l(p)*P_l(d)*G2q_nk[0]*F1p_nk[0]*F1kmp_nk[0]);

            return pow2(omega_m)*(terms[0] + 2.*omega_q/omega_m*(1.+cs2)*terms[2] + pow2(omega_q/omega_m*(1.+cs2))*terms[1]);
            break;
        case 9:
          // d_m d_m d_m
          return    norm*(P_l(k)*P_l(p)*F1_nk*F1p_nk[0]*F2A_nk[0] +  P_l(d)*P_l(k)*F1_nk*F2D_nk[0]*F1kmp_nk[0]
                        + P_l(p)*P_l(d)*F2_nk[0]*F1p_nk[0]*F1kmp_nk[0]);
          break;
      default:
            warning("SPT: invalid indices, a = %d\n", a);
            return 0;
    }
}

// 1-loop resummed spectrum

// Spline linear growth factors and no wiggle spectra
Spline F1_spline, G1_spline, EHUNWcde, RENWcde;
double ANWcde,sigmav;


double KMIN_cde = 1e-4;
double KMAX_cde = 50.;
void CDE::F1_init(int a, double scalef, double omega0, double par1, double par2, double par3) const{
  IOW iow;
  vector<double> kval_table, F1_table, G1_table;
  int n3 = 300;
  for(int i = 0; i<n3; i++){
  double k = KMIN_cde*exp(i*log(KMAX_cde/(KMIN_cde))/(n3-1.));
  iow.initn_lin(a, scalef, k, omega0, par1, par2, par3);
  double ling1=F1_nk;
  double ling2=G1_nk;

  kval_table.push_back(k);
  F1_table.push_back(ling1);
  G1_table.push_back(ling2);

  }
  F1_spline = LinearSpline(kval_table,F1_table);
  G1_spline = LinearSpline(kval_table,G1_table);
}

static double ANWcde_integrand(double q, double k){
return pow2(F1_spline(k))*RENWcde(k)*(1.-j0(k*q))*pow2(q);
}

static double nowig_integrand(const PowerSpectrum& P_l, double k, double q){
double lambda = 0.25*pow(k/0.05,0.04);
return P_l(q)/EHUNWcde(q)*exp(-pow2(log(k/q)/lambda)/2.)/q;
}


////////////////
// FOR RSD ONLY

static double cde_exp(const PowerSpectrum& P_l, double q){
  return  P_l(q)*pow2(F1_spline(q))/(6.*pow2(M_PI));
}

static void sigv_init(const PowerSpectrum& P_l){
    sigmav = Integrate(bind(cde_exp,cref(P_l), _1), KMIN_cde, KMAX_cde, 1e-3);
  }
//////////////


// Initialize the NW power spectrum, A^{nw,0l} and sigma_v

// run this before running other kernel initializations (initn or initn_cde)
// a chooses cde(1) or mg(0)
void CDE::initcde(int a, double scalef, double omega0, double par1, double par2, double par3) const{
NoWigglePS now(C, 0. , EisensteinHu);
IOW iow;

// linear spectrum initialization in range [1e-4,50]
F1_init(a, scalef, omega0, par1, par2, par3);

// sigma_v iinitialization
//sigv_init(cref(P_l));

// EHUNWcde initialization
int n3 = 500;
vector<double> kval_table, ehu_table, now_table;;
for(int i = 0; i<n3; i++ ){
double  k = KMIN_cde*exp(i*log(KMAX_cde/KMIN_cde)/(n3-1.));
double  ehups = now.Evaluate(k);
  kval_table.push_back(k);
  ehu_table.push_back(ehups);
    }

EHUNWcde = LinearSpline(kval_table,ehu_table);

// Renormalized linear spectrum and ANWcde initialization
const double qmin = 10.;
const double qmax = 300.;
double c[2] = {qmin,KMIN_cde};
double d[2] = {qmax,KMAX_cde};

for(int i = 0; i<n3; i++ ){
  double k = kval_table[i];

  double lambda = 0.25*pow(k/0.05,0.04);
  double nowps = EHUNWcde(k)/sqrt(2*M_PI*pow2(lambda))*Integrate<ExpSub>(bind(nowig_integrand,cref(P_l),k,_1),KMIN_cde,KMAX_cde,epsrel);
  now_table.push_back(nowps);
    }

  RENWcde = LinearSpline(kval_table,now_table);
  ANWcde = 1./((pow3(qmax)-pow3(qmin))*pow2(M_PI))*Integrate<2>(bind(ANWcde_integrand,_1,_2),c,d,1e-3,1e-3);

}


// NUMERICAL KERNEL RESUMMED 1-LOOP TERMS
//Integrating over angle - P22 numerical:
static double F2F(int y, int a,  const PowerSpectrum& P_l, double k, double r){
	double temp_ps;
	double myresult=0.;
	switch (a) {
			case 1:
			for( int i = 0; i < n1; i++ )
				{
        double d = 1+ r*r - 2*r*x128[i];
		 		if(d < 1e-5){
			 	temp_ps=0.;
		 		}
		  	else {
		  	temp_ps = RENWcde(k*sqrt(d)) * pow2(F2_nk[i*n2 + y]);
				}
        myresult += w128[i] * temp_ps;
				}
        return 2. * r * r * myresult;
        break;
			case 2:
      myresult = 0.;
			for( int i = 0; i < n1; i++ )
				{
        double d = 1+ r*r - 2*r*x128[i];
				if(d < 1e-5){
				temp_ps=0.;
				}
				else {
				temp_ps = RENWcde(k*sqrt(d)) * G2_nk[i*n2 + y]*F2_nk[i*n2 + y];
				}
        myresult += w128[i] * temp_ps;
				}
        return  2. * r * r * myresult;
        break;
			case 3:
			for( int i = 0; i < n1; i++ )
				{
				double d = 1+ r*r - 2*r*x128[i];
				if(d < 1e-5){
				temp_ps=0.;
				}
				else {
				temp_ps = RENWcde(k*sqrt(d)) * pow2(G2_nk[i*n2 + y]);
				}
        myresult += w128[i] * temp_ps;
				}
        return  2. * r * r * myresult;
        break;
			}
	}

  //Integrating over angle -P13 numerical:
  static double F3F(int y, int a,  const PowerSpectrum& P_l, double k, double r){
    double myresult = 0.;
    // Set the limits of angular integrations
  	switch (a) {
  			case 1:
  			 for( int i = 0; i < n1; i++ )
  				{
  				myresult += w128[i] * F1_nk * F3_nk[i*n2 + y];
  				}
        return   6. * r * r * myresult;
        break;
  			case 2:
  			 for( int i = 0; i < n1; i++ )
  				{
  					myresult += w128[i] * (G1_nk*F3_nk[i*n2 + y] + F1_nk*G3_nk[i*n2 + y]);
  				}
        return   3. * r * r * myresult;
        break;
  			case 3:
  				for( int i = 0; i < n1; i++ )
  				{
  					myresult += w128[i] * G1_nk * G3_nk[i*n2 + y];
  				}
          return    6. * r * r * myresult;
          break;

  		}
    }

/*Selection function for numerical kernel magnitude */

inline int num_selec_mag(double kmin, double kmax, double y){
      double YMIN =QMINp/kmax;
      double YMAX = QMAXp/kmin;
    	int mag_int;
  //  	mag_int = (y-YMIN)*(n2*1.-1)/(YMAX-YMIN); // linear
    //mag_int = sqrt((y-QMINp/kmax)*kmin/QMAXp)*(n2-1); //quadratic
    mag_int = (int)round((n2-1)*log(y/YMIN)/log(YMAX/YMIN)); //exponential
    		return mag_int;
}

/* P^{(22)} */
// a = 1 : delta delta
// a = 2 : delta theta
// a = 3 : theta theta
double CDE::P22nres(double kmin, double kmax,  int a, double k) const {
    	double KMAX = QMAXp/k;
    	double KMIN = QMINp/k;
      int y1 = num_selec_mag(kmin, kmax, KMIN);
      int y2 = num_selec_mag(kmin, kmax, KMAX);
      double y[n2];
      double integrand[n2];
      for (int i = y1; i<=y2; i++){
      y[i] = QMINp/kmax * exp(i*log(QMAXp*kmax/(QMINp*kmin))/(n2*1.-1.)); // exponential sampling
  //    y[i] = i*1./(n2-1.)*(QMAXp/kmin-QMINp/kmax)+QMINp/kmax; // linear sampling
      integrand[i] = RENWcde(k*y[i]) * F2F(i, a, cref(P_l), k, y[i]);
}
double res = 0.;
  for( int i = y1+1; i <= y2; ++ i ){
res += 0.5 * (y[i] - y[i-1])*(integrand[i] + integrand[i-1]);
}
return  k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt) * res;
}

double CDE::P13nres(double kmin, double kmax,  int a, double k) const {
  double KMAX = QMAXp/k;
  double KMIN = QMINp/k;
  int y1 = num_selec_mag(kmin, kmax, KMIN);
  int y2 = num_selec_mag(kmin, kmax, KMAX);
  double y[n2];
  double integrand[n2];
  for (int i = y1; i<=y2; i++){
  y[i] = QMINp/kmax * exp(i*log(QMAXp*kmax/(QMINp*kmin))/(n2*1.-1.));
//  y[i] = i*1./(n2-1.)*(QMAXp/kmin-QMINp/kmax)+QMINp/kmax; // linear sampling
  integrand[i] = RENWcde(k*y[i]) * F3F(i, a, cref(P_l), k, y[i]);
  }
double res = 0.;
  for( int i = y1+1; i <= y2; ++ i ){
res += 0.5 * (y[i] - y[i-1])*(integrand[i] + integrand[i-1]);
}
  return  k*k*k/(4*M_PI*M_PI)/pow4(dnorm_spt) * RENWcde(k) * res;
}


double CDE::PNWloopn(double kmin, double kmax, double k, int a ) const{
  switch(a) {
      case 1:
          return   pow2(F1_spline(k)/dnorm_spt)*RENWcde(k) + P13nres(kmin, kmax, a, k) + P22nres(kmin, kmax, a, k);
          break;
      case 2:
          return  F1_spline(k)*G1_spline(k)/pow2(dnorm_spt)*RENWcde(k) + P13nres(kmin, kmax, a, k) + P22nres(kmin, kmax, a, k);
          break;
      case 3:
          return  pow2(G1_spline(k)/dnorm_spt)*RENWcde(k) + P13nres(kmin, kmax, a, k) + P22nres(kmin, kmax, a, k);
          break;
      default:
          warning("CDE: invalid indices, a = %d\n", a);
          return 0;
  }
}

/* Resummed 1-loop spectra a la Fonseca */

double CDE::Presumn(double kmin, double kmax, double k, int a) const{
  SPT spt(C,P_l,epsrel);
  switch(a){
  case 1:
   return   PNWloopn(kmin,kmax,k,1) + exp(-0.5*pow2(k)*ANWcde/pow2(dnorm_spt))*(spt.PLOOPn(kmin,kmax,1,k)-PNWloopn(kmin,kmax,k,1) + 0.5*pow2(k)*ANWcde/pow2(dnorm_spt)*pow2(F1_spline(k)/dnorm_spt)*(P_l(k) - RENWcde(k)));
      break;
  case 2:
   return PNWloopn(kmin,kmax,k,2) + exp(-0.5*pow2(k)*ANWcde/pow2(dnorm_spt))*(spt.PLOOPn(kmin,kmax,2,k)-PNWloopn(kmin,kmax,k,2) + 0.5*pow2(k)*ANWcde/pow2(dnorm_spt)*G1_spline(k)*F1_spline(k)/pow2(dnorm_spt)*(P_l(k) - RENWcde(k)));
      break;
  case 3:
      return PNWloopn(kmin,kmax,k,3) + exp(-0.5*pow2(k)*ANWcde/pow2(dnorm_spt))*(spt.PLOOPn(kmin,kmax,3,k)-PNWloopn(kmin,kmax,k,3) + 0.5*pow2(k)*ANWcde/pow2(dnorm_spt)*pow2(G1_spline(k)/dnorm_spt)*(P_l(k) - RENWcde(k)));
      break;
  case 4:
      return pow2(F1_spline(k)/dnorm_spt)*RENWcde(k);
      //IR CHECK
  case 5:
      return P13nres(kmin, kmax, 1, k) + P22nres(kmin, kmax, 1, k) + exp(-0.5*pow2(k)*ANWcde/pow2(dnorm_spt))*(spt.P22n(kmin,kmax,1,k) + spt.P13n(kmin,kmax,1,k) -(P13nres(kmin, kmax, 1, k) + P22nres(kmin, kmax, 1, k)));
  default:
          warning("CDE: invalid indices, a = %d\n", a);
          return 0;
}
}
