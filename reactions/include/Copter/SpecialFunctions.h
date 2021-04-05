#ifndef SPECIALFUNCTIONS_H
#define SPECIALFUNCTIONS_H

#include <vector>

// GL quadrature for 256 abscissae - used for numerical kernel angular integration
// abscissae
const double x128[256] = {-9.999561e-01,-9.997684e-01,-9.994309e-01,-9.989435e-01,-9.983063e-01,-9.975193e-01,-9.965826e-01,-9.954965e-01,-9.942610e-01,-9.928763e-01,-9.913428e-01,-9.896605e-01,-9.878297e-01,-9.858508e-01,-9.837240e-01,-9.814496e-01,-9.790280e-01,-9.764595e-01,-9.737446e-01,-9.708836e-01,-9.678769e-01,-9.647251e-01,-9.614285e-01,-9.579877e-01,-9.544032e-01,-9.506755e-01,-9.468052e-01,-9.427929e-01,-9.386392e-01,-9.343446e-01,-9.299099e-01,-9.253357e-01,-9.206227e-01,-9.157716e-01,-9.107831e-01,-9.056580e-01,-9.003970e-01,-8.950010e-01,-8.894707e-01,-8.838069e-01,-8.780106e-01,-8.720826e-01,-8.660238e-01,-8.598350e-01,-8.535173e-01,-8.470715e-01,-8.404987e-01,-8.337997e-01,-8.269757e-01,-8.200277e-01,-8.129566e-01,-8.057636e-01,-7.984497e-01,-7.910160e-01,-7.834637e-01,-7.757938e-01,-7.680076e-01,-7.601062e-01,-7.520907e-01,-7.439624e-01,-7.357225e-01,-7.273723e-01,-7.189129e-01,-7.103457e-01,-7.016719e-01,-6.928929e-01,-6.840099e-01,-6.750243e-01,-6.659375e-01,-6.567508e-01,-6.474655e-01,-6.380831e-01,-6.286050e-01,-6.190327e-01,-6.093674e-01,-5.996107e-01,-5.897641e-01,-5.798290e-01,-5.698070e-01,-5.596994e-01,-5.495079e-01,-5.392340e-01,-5.288792e-01,-5.184450e-01,-5.079331e-01,-4.973450e-01,-4.866822e-01,-4.759465e-01,-4.651394e-01,-4.542624e-01,-4.433174e-01,-4.323058e-01,-4.212294e-01,-4.100898e-01,-3.988887e-01,-3.876278e-01,-3.763087e-01,-3.649331e-01,-3.535028e-01,-3.420195e-01,-3.304849e-01,-3.189007e-01,-3.072686e-01,-2.955905e-01,-2.838680e-01,-2.721029e-01,-2.602971e-01,-2.484521e-01,-2.365699e-01,-2.246523e-01,-2.127009e-01,-2.007176e-01,-1.887042e-01,-1.766625e-01,-1.645943e-01,-1.525014e-01,-1.403856e-01,-1.282488e-01,-1.160927e-01,-1.039192e-01,-9.173013e-02,-7.952729e-02,-6.731252e-02,-5.508766e-02,-4.285453e-02,-3.061497e-02,-1.837082e-02,-6.123912e-03, 0.0061239123751895295011702,0.0183708184788136651179263,0.0306149687799790293662786,0.0428545265363790983812423,0.0550876556946339841045614,0.0673125211657164002422903,0.0795272891002329659032271,0.0917301271635195520311456,0.1039192048105094036391969,0.1160926935603328049407349,0.1282487672706070947420496,0.1403856024113758859130249,0.1525013783386563953746068,0.1645942775675538498292845,0.1766624860449019974037218,0.1887041934213888264615036,0.2007175933231266700680007,0.2127008836226259579370402,0.2246522667091319671478783,0.2365699497582840184775084,0.2484521450010566668332427,0.2602970699919425419785609,0.2721029478763366095052447,0.2838680076570817417997658,0.2955904844601356145637868,0.3072686197993190762586103,0.3189006618401062756316834,0.3304848656624169762291870,0.3420194935223716364807297,0.3535028151129699895377902,0.3649331078236540185334649,0.3763086569987163902830557,0.3876277561945155836379846,0.3988887074354591277134632,0.4100898214687165500064336,0.4212294180176238249768124,0.4323058260337413099534411,0.4433173839475273572169258,0.4542624399175899987744552,0.4651393520784793136455705,0.4759464887869833063907375,0.4866822288668903501036214,0.4973449618521814771195124,0.5079330882286160362319249,0.5184450196736744762216617,0.5288791792948222619514764,0.5392340018660591811279362,0.5495079340627185570424269,0.5596994346944811451369074,0.5698069749365687590576675,0.5798290385590829449218317,0.5897641221544543007857861,0.5996107353629683217303882,0.6093674010963339395223108,0.6190326557592612194309676,0.6286050494690149754322099,0.6380831462729113686686886,0.6474655243637248626170162,0.6567507762929732218875002,0.6659375091820485599064084,0.6750243449311627638559187,0.6840099204260759531248771,0.6928928877425769601053416,0.7016719143486851594060835,0.7103456833045433133945663,0.7189128934599714483726399,0.7273722596496521265868944,0.7357225128859178346203729,0.7439624005491115684556831,0.7520906865754920595875297,0.7601061516426554549419068,0.7680075933524456359758906,0.7757938264113257391320526,0.7834636828081838207506702,0.7910160119895459945467075,0.7984496810321707587825429,0.8057635748129986232573891,0.8129565961764315431364104,0.8200276660989170674034781,0.8269757238508125142890929,0.8337997271555048943484439,0.8404986523457627138950680,0.8470714945172962071870724,0.8535172676795029650730355,0.8598350049033763506961731,0.8660237584665545192975154,0.8720825999954882891300459,0.8780106206047065439864349,0.8838069310331582848598262,0.8894706617776108888286766,0.8950009632230845774412228,0.9003970057703035447716200,0.9056579799601446470826819,0.9107830965950650118909072,0.9157715868574903845266696,0.9206227024251464955050471,0.9253357155833162028727303,0.9299099193340056411802456,0.9343446275020030942924765,0.9386391748378148049819261,0.9427929171174624431830761,0.9468052312391274813720517,0.9506755153166282763638521,0.9544031887697162417644479,0.9579876924111781293657904,0.9614284885307321440064075,0.9647250609757064309326123,0.9678769152284894549090038,0.9708835784807430293209233,0.9737445997043704052660786,0.9764595497192341556210107,0.9790280212576220388242380,0.9814496290254644057693031,0.9837240097603154961666861,0.9858508222861259564792451,0.9878297475648606089164877,0.9896604887450652183192437,0.9913427712075830869221885,0.9928763426088221171435338,0.9942609729224096649628775,0.9954964544810963565926471,0.9965826020233815404305044,0.9975192527567208275634088,0.9983062664730064440555005,0.9989435258434088565550263,0.9994309374662614082408542,0.9997684374092631861048786,0.9999560500189922307348012};
// weights
const double w128[256] = {1.127890e-04,2.625349e-04,4.124633e-04,5.623490e-04,7.121542e-04,8.618537e-04,1.011424e-03,1.160844e-03,1.310089e-03,1.459137e-03,1.607967e-03,1.756556e-03,1.904881e-03,2.052920e-03,2.200652e-03,2.348053e-03,2.495102e-03,2.641777e-03,2.788055e-03,2.933916e-03,3.079336e-03,3.224294e-03,3.368769e-03,3.512738e-03,3.656180e-03,3.799074e-03,3.941398e-03,4.083130e-03,4.224250e-03,4.364737e-03,4.504569e-03,4.643725e-03,4.782184e-03,4.919926e-03,5.056930e-03,5.193175e-03,5.328642e-03,5.463309e-03,5.597156e-03,5.730164e-03,5.862312e-03,5.993581e-03,6.123951e-03,6.253402e-03,6.381915e-03,6.509470e-03,6.636050e-03,6.761633e-03,6.886203e-03,7.009739e-03,7.132224e-03,7.253639e-03,7.373966e-03,7.493186e-03,7.611283e-03,7.728238e-03,7.844033e-03,7.958652e-03,8.072077e-03,8.184291e-03,8.295278e-03,8.405020e-03,8.513501e-03,8.620705e-03,8.726616e-03,8.831218e-03,8.934495e-03,9.036432e-03,9.137013e-03,9.236223e-03,9.334048e-03,9.430473e-03,9.525483e-03,9.619065e-03,9.711203e-03,9.801885e-03,9.891096e-03,9.978823e-03,1.006505e-02,1.014977e-02,1.023297e-02,1.031464e-02,1.039475e-02,1.047331e-02,1.055029e-02,1.062570e-02,1.069950e-02,1.077171e-02,1.084230e-02,1.091126e-02,1.097858e-02,1.104426e-02,1.110828e-02,1.117063e-02,1.123131e-02,1.129031e-02,1.134761e-02,1.140321e-02,1.145709e-02,1.150926e-02,1.155970e-02,1.160841e-02,1.165538e-02,1.170060e-02,1.174406e-02,1.178576e-02,1.182570e-02,1.186386e-02,1.190024e-02,1.193483e-02,1.196764e-02,1.199865e-02,1.202785e-02,1.205526e-02,1.208086e-02,1.210464e-02,1.212661e-02,1.214676e-02,1.216509e-02,1.218159e-02,1.219626e-02,1.220911e-02,1.222012e-02,1.222930e-02,1.223665e-02,1.224216e-02,1.224583e-02,1.224767e-02,0.0122476716402897559040703,0.0122458343697479201424639,0.0122421601042728007697281,0.0122366493950401581092426,0.0122293030687102789041463,0.0122201222273039691917087,0.0122091082480372404075141,0.0121962627831147135181810,0.0121815877594817721740476,0.0121650853785355020613073,0.0121467581157944598155598,0.0121266087205273210347185,0.0121046402153404630977578,0.0120808558957245446559752,0.0120552593295601498143471,0.0120278543565825711612675,0.0119986450878058119345367,0.0119676359049058937290073,0.0119348314595635622558732,0.0119002366727664897542872,0.0118638567340710787319046,0.0118256971008239777711607,0.0117857634973434261816901,0.0117440619140605503053767,0.0117005986066207402881898,0.0116553800949452421212989,0.0116084131622531057220847,0.0115597048540436357726687,0.0115092624770394979585864,0.0114570935980906391523344,0.0114032060430391859648471,0.0113476078955454919416257,0.0112903074958755095083676,0.0112313134396496685726568,0.0111706345765534494627109,0.0111082800090098436304608,0.0110442590908139012635176,0.0109785814257295706379882,0.0109112568660490397007968,0.0108422955111147959952935,0.0107717077058046266366536,0.0106995040389797856030482,0.0106256953418965611339617,0.0105502926865814815175336,0.0104733073841704030035696,0.0103947509832117289971017,0.0103146352679340150682607,0.0102329722564782196569549,0.0101497741990948656546341,0.0100650535763063833094610,0.0099788230970349101247339,0.0098910956966958286026307,0.0098018845352573278254988,0.0097112029952662799642497,0.0096190646798407278571622,0.0095254834106292848118297,0.0094304732257377527473528,0.0093340483776232697124660,0.0092362233309563026873787,0.0091370127604508064020005,0.0090364315486628736802278,0.0089344947837582075484084,0.0088312177572487500253183,0.0087266159616988071403366,0.0086207050884010143053688,0.0085135010250224906938384,0.0084050198532215357561803,0.0082952778462352254251714,0.0081842914664382699356198,0.0080720773628734995009470,0.0079586523687543483536132,0.0078440334989397118668103,0.0077282379473815556311102,0.0076112830845456594616187,0.0074931864548058833585998,0.0073739657738123464375724,0.0072536389258339137838291,0.0071322239610753900716724,0.0070097390929698226212344,0.0068862026954463203467133,0.0067616333001737987809279,0.0066360495937810650445900,0.0065094704150536602678099,0.0063819147521078805703752,0.0062534017395424012720636,0.0061239506555679325423891,0.0059935809191153382211277,0.0058623120869226530606616,0.0057301638506014371773844,0.0055971560336829100775514,0.0054633085886443102775705,0.0053286415939159303170811,0.0051931752508692809303288,0.0050569298807868423875578,0.0049199259218138656695588,0.0047821839258926913729317,0.0046437245556800603139791,0.0045045685814478970686418,0.0043647368779680566815684,0.0042242504213815362723565,0.0040831302860526684085998,0.0039413976414088336277290,0.0037990737487662579981170,0.0036561799581425021693892,0.0035127377050563073309711,0.0033687685073155510120191,0.0032242939617941981570107,0.0030793357411993375832054,0.0029339155908297166460123,0.0027880553253277068805748,0.0026417768254274905641208,0.0024951020347037068508395,0.0023480529563273120170065,0.0022006516498399104996849,0.0020529202279661431745488,0.0019048808534997184044191,0.0017565557363307299936069,0.0016079671307493272424499,0.0014591373333107332010884,0.0013100886819025044578317,0.0011608435575677247239706,0.0010114243932084404526058,0.0008618537014200890378141,0.0007121541634733206669090,0.0005623489540314098028152,0.0004124632544261763284322,0.0002625349442964459062875,0.0001127890178222721755125};

// GL quadrature for 64 (proper) abscissae - used for numerical kernel angular integration
const double x32[32] = {-0.997264, -0.985612, -0.964762, -0.934906, -0.896321, -0.849368,
-0.794484, -0.732182, -0.663044, -0.587716, -0.5069, -0.421351,
-0.331869, -0.239287, -0.144472, -0.0483077, 0.0483077, 0.144472,
0.239287, 0.331869, 0.421351, 0.5069, 0.587716, 0.663044, 0.732182,
0.794484, 0.849368, 0.896321, 0.934906, 0.964762, 0.985612, 0.997264};
const double w32[32] = {0.00701862, 0.0162744, 0.0253919, 0.0342739, 0.0428358, 0.0509981,
0.0586841, 0.0658222, 0.0723458, 0.0781939, 0.0833119, 0.0876521,
0.0911739, 0.0938444, 0.0956387, 0.0965401, 0.0965401, 0.0956387,
0.0938444, 0.0911739, 0.0876521, 0.0833119, 0.0781939, 0.0723458,
0.0658222, 0.0586841, 0.0509981, 0.0428358, 0.0342739, 0.0253919,
0.0162744, 0.00701862};


// angular integration limits
const real XMAX = x128[255];
const real XMIN = x128[0];

// PS integration magnitude limits
const real QMINp = 1e-4;
const real QMAXp = 30.;


//Kernel array declation

extern double F1_nk;
extern double G1_nk;

extern double F1q_nk;
extern double G1q_nk;

extern double F1kmp_nk;
extern double G1kmp_nk;

extern double F1p_nk;
extern double G1p_nk;

extern double F1p_nkp;
extern double G1p_nkp;

extern double F2A_nk;
extern double G2A_nk;

extern double F2B_nk;
extern double G2B_nk;

extern double F2C_nk;
extern double G2C_nk;

extern double F2D_nk;
extern double G2D_nk;

extern double F2_nk;
extern double G2_nk;

extern double F2_nkp;
extern double G2_nkp;

extern double F3_nk;
extern double G3_nk;

extern double F3_nkp;
extern double G3_nkp;


  ///// EVOLUTION FACTORS ////////

  // LCDM
extern double Dl_spt;
  // DGP
extern double D_spt; //
extern double F_spt; //
extern double Cx_spt;
extern double I_spt;
extern double J_spt;
extern double K_spt;
extern double KL_spt;
extern double F3_spt; //


  //derivatives of evolution factors
  //LCDM
 extern double fl_spt;
  //DGP
 extern double fdgp_spt;
 extern double Dd_spt;
 extern double Fd_spt;
 extern double Cdx_spt;
 extern double Id_spt;
 extern double Jd_spt;
 extern double Kd_spt;
 extern double KLd_spt;
 extern double F3d_spt;

 // 4th order evolution factors for DGP
 extern double evol4[34];
 extern double devol4[34];

//Hubble at chosen scale factor
 extern double H_spt;

// Normalisation of linear power spectrum (Dl @ a=1)
 extern double dnorm_spt;


// correction to evolving DE virial concentration (see HALO.cpp, cvirial)
extern double g_de;

class IOW {
public:

  	// Initialisation of evolution factors for analytic kernels as well as linear power spectrum normalisation (dnorm_spt)
  	// A is the scale factor
  	// omega0 is the matter density param
  	// par1,par2,par3 are the mg parameters (generally only par1 is used)
    // model selects MG or DE model (1 = GR, MG: 2 = f(R), 3 = DGP, DE models: 4 = quintessence, 5 = CPL, 6 = HyP)
  	void inite(double A, double omega0, double par1, double par2, double par3, int model );
    // Evolution factors including those for f4/g4 DGP KERNELS (slower than inite)
    void inite2(double A, double omega0, double par1, double par2, double par3, int model );

// numerical kernel initialisation functions
    void initn_lin(double A, double k, double omega0, double par1, double par2, double par3, int model ); // numerical linear growth function solver
    void initnorm(double vars[], int model ); // wCDM and LCDM power spectrum normalisation  - roughly the same as above function but with no MG and with dark energy. It solves at a=1 to normalise general spectrum as well solves with LCDM at given a
    void initn2(double A, double k[], double x[], double kargs[], double omega0, double par1, double par2, double par3, double omeganu, int model ); // 1-loop kernel solver
    void initn2_pseudo(double A, double k[], double x[], double kargs[], double omega0, double par1, double par2, double par3, double omeganu, int model ); // 1-loop pseudo spectrum kernel solver
    void initn3(double redshifts[], int noz, double k[], double x[], double kargs[], double omega0, double par1, double par2, double par3, double mykernelarray[][20], int model ); // 1-loop kernel (real and pseudo) solver for array of redshifts -- used in reaction for lensing
    void initn_rsd(double A, double k[], double x[], double kargs[], double omega0, double par1, double par2, double par3, int model ); // 1-loop TNS solver (also solves for ABC correction term kernels)

};


// hubble functions ---  lcdm
double HA(double a, double omega0);
double HA1(double a, double omega0);
double HA2(double a, double omega0);

// hubble functions ---   general
double HAg(double a, double omega0, double p1, double p2,  double p3, int model );
double HA1g(double a, double omega0, double p1, double p2,  double p3, int model );
double HA2g(double a, double omega0, double p1, double p2,  double p3, int model );
double HA2g2(double a, double omega0, double p1, double p2,  double p3, int model );
// Spherical collapse contribution of DE
double WEFF(double a, double omega0, double p1, double p2, double p3, int model );

// Modified gravity functions
// 1-loop PT (1606.02520)
double mu(double a, double k0, double omega0, double p1, double p2, double p3, int model  );
double gamma2(double a, double omega0, double k0, double k1, double k2, double u1, double p1, double p2, double p3 , int model );
double gamma3(double a, double omega0, double k0, double k1, double k2, double k3, double u1,double u2, double u3, double p1, double p2, double p3, int model );
// spherical collapse (1812.05594 - appendix)
double mymgF(double a, double yh, double yenv, double Rth, double omega0, double p1, double p2, double p3, double delta_initial, int model);


// standard PT kernels
double alpha(double k1, double k2, double u1);
double alphas(double k1, double k2, double u1);
double beta1(double k1, double k2, double u1);


//// EINSTEIN DE SITTER THEORY KERNEL FUNCTIONS /////

/*EDS Kernels for DGP and GR*/

//2nd order LCDM
double F2eds(double k1, double k2, double u1);
double G2eds(double k1, double k2, double u1);

// 2nd order DGP
double F2ndgp(double k1, double k2, double u1);
double G2ndgp(double k1, double k2, double u1);

//3rd order LCDM
double F3eds(double k1, double k2, double k3, double u1, double u2, double u3);
double G3eds(double k1, double k2, double k3, double u1, double u2, double u3);
// optimised input
double F3edsb(double k1, double k2, double k3, double k4, double k5, double k6, double x1, double x2, double x3);
double G3edsb(double k1, double k2, double k3, double k4, double k5, double k6, double x1, double x2, double x3);

//3rd order DGP
double F3ndgp(double k1, double k2, double k3, double x1, double x2, double x3);
double G3ndgp(double k1, double k2, double k3, double x1, double x2, double x3);
// optimised input
double F3bdgp(double k[],double x[]);
double G3bdgp(double k[],double x[]);

//4th order LCDM
double F4eds(double k[], double x[]);
double G4eds(double k[], double x[]);
// optimised input
double F4edsb(double k1, double k2, double k3, double k4, double x1, double x2,double x3, double x4, double x5);
double G4edsb(double k1, double k2, double k3, double k4, double x1, double x2,double x3, double x4, double x5);

// 4th order DGP
double F4dgp(double k[], double x[]);
double G4dgp(double k[], double x[]);

// F2 fitting functions as in Gil-Marin and Scoccimaro Bispectrum papers
double F2fit(double vars[], double k1, double k2, double x);
double F2fitsc(double vars[], double k1, double k2, double x);

/* multipole factors : Example for TNS: (2l+1)/2 * integ dmu  e^(-(k*f*sigma_v*mu)^2) * mu^(2n) * P_l(mu) */
// a =1 : monopole
// a =2 : quadrupole
// a =3 : hexdecapole
// e selects multipole factor (1-3: Resummed, 4-5: Matsubara, 6:TNS)
double factL(real k, real sigma_v, real F0, real anw, int n, int a, int e );

/* Fingers of god term in exponential form */
// a =1 : LCDM
// a =2 : nDGP
// a =3 : numerical (arbitrary model)
double DFOG(real k, real u,  double sigma_v, int a);

/* Minimum and maximum functions */
double Min(real a, real b);
double Max(real a, real b);

/* Cylindrical Bessel functions */
double BesselJ0(double x);
double BesselJ1(double x);
double BesselJn(int n, double x);

// Legendre polynomials
double legendre(int order, double x );

/* Spherical Bessel functions */
double SphericalBesselJ0(double x);
double SphericalBesselJ1(double x);
double SphericalBesselJ2(double x);
double SphericalBesselJ3(double x);
double SphericalBesselJ4(double x);


/* Gamma function:
 *   \Gamma(a) = \int_0^\infty t^{a-1} e^{-t} dt */
double Gamma(double a);

/* Natural logarithm of gamma function:
 *   \log \Gamma(a) */
double LogGamma(double a);

/* Unnormalized lower incomplete gamma function:
 *  \gamma(a,x) = \int_0^x t^{a-1} e^{-t} dt */
double LowerGamma(double a, double x);

/* Unnormalized upper incomplete gamma function:
 *  \Gamma(a,x) = \int_x^\infty t^{a-1} e^{-t} dt */
double UpperGamma(double a, double x);

class integral
{
    ////////// Con-/destructor & initializer//////////
public:
    integral(  );
    ~integral(  );
    void clear(  );
    void clear( const unsigned & size );

    ////////// Data buffer //////////
private:
    std::vector<double> x_buf, y_buf;

    ////////// Read data //////////
public:
    void read( const double & x, const double & y );

    ////////// Feed back //////////
public:
    double result(  );

    ////////// Gauss-Legendre //////////
public:
    void   gl_clear (  );
    double gl_xi    ( const int & i );
    void   gl_read  ( const int & i, const double & kernel );
    double gl_result(  );
    static const int gl_num;
public:
    double gl_intg_res;
    static const double x[  ];
    static const double w[  ];
};

#endif // SPECIALFUNCTIONS_H
