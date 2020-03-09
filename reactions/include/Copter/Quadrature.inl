#include <cassert>
#include <cfloat>
#include <cstdlib>
#include <map>

#include "Common.h"

const int GK_LIMIT = 8192;

const int GM_MINPTS = 1000;
const int GM_MAXPTS = 10000000;


/* GKWorkspace: Workspace for adaptive Gauss-Kronrod algorithm */
struct GKWorkspace {
    const int limit;
    int size;
    double* alist;
    double* blist;
    double* rlist;
    double* elist;
    int* order;

    /* Initialize workspace.  limit is the maximum number of subintervals allowed. */
    GKWorkspace(int limit = GK_LIMIT);
    ~GKWorkspace();

    void GetMaxInterval(double& a, double& b, double& r, double& e);
    void Sort();
    void Update(double a1, double b1, double area1, double error1,
                double a2, double b2, double area2, double error2);
    void SumResults(double& result, double& abserr);

    static void ReportError(int error_code);

    /* General Gauss-Kronrod rule */
    static double GaussKronrod(const int n, const double wg[], const double wgk[],
                               double fc, double fv1[], double fv2[], double half,
                               double* abserr, double* resabs, double* resasc);

    /* 15-point Gauss-Kronrod rule with error estimate */
    template<typename Function>
    static double GaussKronrod15(Function f, double a, double b, double* abserr, double* resabs, double* resasc);

    /* Integrate f(x) from a to b */
    template<typename Function>
    double Integrate(Function f, double a, double b, double epsrel, double epsabs, double* abserr, int* neval);

    static const double wg15[4];
    static const double wgk15[8];
    static const double xgk15[8];
};


/***** GKWorkspace *****/

template<typename Function>
double GKWorkspace::GaussKronrod15(Function f, double a, double b,
                                    double* abserr, double* resabs, double* resasc)
{
    const int n = 8;
    double fv1[n], fv2[n];
    const double c = 0.5*(a + b);         // center of the interval
    const double fc = f(c);               // f evaluated at c
    const double half = 0.5*(b - a);      // half the interval
    double abscissa;
    for(int j = 0; j < (n-1)/2; j++) {
        abscissa = half*xgk15[2*j+1];
        fv1[2*j+1] = f(c - abscissa);
        fv2[2*j+1] = f(c + abscissa);
    }
    for(int j = 0; j < n/2; j++) {
        abscissa = half*xgk15[2*j];
        fv1[2*j] = f(c - abscissa);
        fv2[2*j] = f(c + abscissa);
    }
    return GaussKronrod(n, wg15, wgk15, fc, fv1, fv2, half, abserr, resabs, resasc);
}

/* Adaptive integration using GaussKronrod15 for fundamental intervals */
template<typename Function>
double GKWorkspace::Integrate(Function f, double a, double b, double epsrel, double epsabs, double* pabserr, int* pneval) {
    int error_code = 0;
    int neval = 0;
    double result = 0, abserr = 0, resabs = 0, resasc = 0;

    /* Make sure precision request is reasonable */
    if(epsabs <= 0 && (epsrel < 50*DBL_EPSILON || epsrel < 0.5e-28)) {
        ReportError(1);
        return 0;
    }

    /* Integrate from left to right, fixing the sign at the end */
    double sign = +1;
    if(a > b) {
        double tmp = a;
        a = b;
        b = tmp;
        sign = -1;
    }

    /* Perform first integration */
    result = GaussKronrod15(f, a, b, &abserr, &resabs, &resasc);
    neval += 15;

    /* Initialize workspace */
    size = 1;
    alist[0] = a;
    blist[0] = b;
    rlist[0] = result;
    elist[0] = abserr;
    order[0] = 0;

    double area = result;
    double errsum = abserr;
    int iteration = 1;

    /* Adaptively subdivide interval until convergence is achieved */
    double a1, b1, a2, b2;
    double a_i, b_i, r_i, e_i;
    double area1, area2, area12;
    double error1, error2, error12;
    double resasc1, resasc2, resabs1, resabs2;
    double tolerance;
    do {
        area1 = area2 = area12 = error1 = error2 = error12 = 0; 

        /* Bisect the interval with the largest error estimate */
        GetMaxInterval(a_i, b_i, r_i, e_i);
        a1 = a_i;
        b1 = 0.5*(a_i + b_i);
        a2 = b1;
        b2 = b_i;

        /* Integrate over the two subintervals */
        area1 = GaussKronrod15(f, a1, b1, &error1, &resabs1, &resasc1);
        area2 = GaussKronrod15(f, a2, b2, &error2, &resabs2, &resasc2);
        neval += 30;
        area12 = area1 + area2;
        error12 = error1 + error2;
        errsum += error12 - e_i;
        area += area12 - r_i;

        tolerance = fmax(epsabs, epsrel*fabs(area));

        if(errsum > tolerance) {
            /* Check for bad behavior at a point of the integration range */
            double tmp = (1 + 100*DBL_EPSILON) * (fabs(a2) + 1000*DBL_MIN);
            if(fabs(a1) <= tmp && fabs(b2) <= tmp)
                error_code = 3;

            /* Check if we've reached the maximum number of subintervals */
            if(iteration == limit)
                error_code = 4;
        }

        Update(a1, b1, area1, error1, a2, b2, area2, error2);
        iteration++;
    }
    while(iteration < limit && error_code == 0 && errsum > tolerance);

    /* Re-sum results to minimize round-off error */
    SumResults(result, abserr);

    if(abserr > tolerance)
        ReportError(error_code);

    if(pabserr)
        *pabserr = abserr;
    if(pneval)
        *pneval = neval;
    return sign*result;
}


/* GMWorkspace: workspace for n-dimensional Genz-Malik algorithm */
struct GMWorkspace {
    int n;              // dimension of integral
    void* wrkstr;       // storage space for region objects

    GMWorkspace(int n);
    ~GMWorkspace();
};


/* Singleton class for managing workspaces.  The rationale here is that we
 * might need to call Integrate() many many times, and we don't want to have to
 * allocate and free workspace memory each time.  So instead we keep a list
 * of already allocated workspaces, and a flag to indicate whether or not it's
 * currently in use.  When a new workspace is requested, we first check for
 * an existing unused workspace.  Only if no workspaces are available do we
 * allocate memory for a new one.  This class handles all of this logic. */
struct WorkspaceManager {
    typedef std::map<GKWorkspace*, int> GKWorkspaceList;
    typedef std::map<GMWorkspace*, int> GMWorkspaceList;

    /* List of instantiated workspaces.  The key is a pointer to a workspace
     * object; the value is either 1 if the workspace is currently in use, or
     * 0 otherwise. */
    GKWorkspaceList gk_workspaces;
    GMWorkspaceList gm_workspaces;

    /* Get or release a GKWorkspace for 1-dimensional integration. */
    GKWorkspace* get_gk_workspace();
    void release_gk_workspace(GKWorkspace* workspace);

    /* Get or release a GMWorkspace for n-dimensional integration. */
    GMWorkspace* get_gm_workspace(int n);
    void release_gm_workspace(GMWorkspace* workspace);

    WorkspaceManager();
    ~WorkspaceManager();
};

/* Global WorkspaceManager object, constructed at startup */
extern WorkspaceManager workspace_manager;


/***** Integrate *****/

template<typename Function>
double Integrate(Function f, double a, double b, double epsrel, double epsabs, double* abserr, int* neval) {
    GKWorkspace* workspace = workspace_manager.get_gk_workspace();
    double result = workspace->Integrate(f, a, b, epsrel, epsabs, abserr, neval);
    workspace_manager.release_gk_workspace(workspace);
    return result;
}


/***** Integrate<Sub> *****/

template<typename Function, typename Sub>
struct SubFunc {
    Function f;
    Sub sub;
    SubFunc(Function f_, Sub sub_ = Sub()) : f(f_), sub(sub_) {}
    double operator()(double u) { return f(sub.x(u)) * sub.dxdu(u); }
};

template<typename Sub, typename Function>
double Integrate(Function f, double a, double b, double epsrel, double epsabs, double* abserr, int* neval, Sub sub) {
    GKWorkspace* workspace = workspace_manager.get_gk_workspace();
    SubFunc<Function, Sub> F(f, sub);
    double result = workspace->Integrate(F, sub.u(a), sub.u(b), epsrel, epsabs, abserr, neval);
    workspace_manager.release_gk_workspace(workspace);
    return result;
}


/***** Integrate<n> *****/

template<int n>
struct region {
    double center[n];
    double width[n];
    double val;
    double err;
    int divaxn;
};

template<typename Function, int n>
struct apply {
    double integrand;
    apply(Function f, double* z) { integrand = f(z); }
    operator double() const { return integrand; }
};

template<typename Function>
struct apply<Function, 1> {
    double integrand;
    apply(Function f, double* z) { integrand = f(z[0]); }
    operator double() const { return integrand; }
};

template<typename Function>
struct apply<Function, 2> {
    double integrand;
    apply(Function f, double* z) { integrand = f(z[0], z[1]); }
    operator double() const { return integrand; }
};

template<typename Function>
struct apply<Function, 3> {
    double integrand;
    apply(Function f, double* z) { integrand = f(z[0], z[1], z[2]); }
    operator double() const { return integrand; }
};

template<typename Function>
struct apply<Function, 4> {
    double integrand;
    apply(Function f, double* z) { integrand = f(z[0], z[1], z[2], z[3]); }
    operator double() const { return integrand; }
};

template<typename Function>
struct apply<Function, 5> {
    double integrand;
    apply(Function f, double* z) { integrand = f(z[0], z[1], z[2], z[3], z[4]); }
    operator double() const { return integrand; }
};

template<typename Function>
struct apply<Function, 6> {
    double integrand;
    apply(Function f, double* z) { integrand = f(z[0], z[1], z[2], z[3], z[4], z[5]); }
    operator double() const { return integrand; }
};

template<int n, typename Function>
double Integrate(Function f, double* a, double* b, double epsrel, double epsabs, double* pabserr, int* pneval) {
    assert(n >= 1 && n <= 15);
    const int two_to_the_n = (1 << n);
    const int rulcls = two_to_the_n + 2*n*n + 2*n + 1;  // number of function calls per basic rule
    assert(GM_MAXPTS >= rulcls);

    /** Initialization **/
    int ifail = 3;
    double finest = 0.;
    double abserr = 0.;
    int funcls = 0;
    int divflg = 1;
    int subrgn = 1;
    int sbrgns = 1;
    int subtmp;
    int maxrgns = 1 + (1 + GM_MAXPTS/rulcls)/2;

    double center[n], width[n], widthl[n], z[n];
    for(int j = 0; j < n; j++) {
        center[j] = 0.5*(a[j] + b[j]);
        width[j] = 0.5*(b[j] - a[j]);
    }

    /** Basic rule initialization **/
    const double lambda2 = sqrt(9./70.);
    const double lambda4 = sqrt(9./10.);
    const double lambda5 = sqrt(9./19.);
    const double wt1 = (12824. - 9120.*n + 400.*n*n)/19683.;
    const double wt2 = 980./6561.;
    const double wt3 = (1820. - 400.*n)/19683.;
    const double wt4 = 200./19683.;
    const double wt5 = 6859./19683./two_to_the_n;
    const double wtp1 = (729. - 950.*n + 50.*n*n)/729.;
    const double wtp2 = 245./486.;
    const double wtp3 = (265. - 100.*n)/1458.;
    const double wtp4 = 25./729.;
    const double ratio = pow2(lambda2/lambda4);

    GMWorkspace* workspace = workspace_manager.get_gm_workspace(n);
    region<n>* wrkstr = (region<n>*)workspace->wrkstr;

    double sum1, sum2, sum3, sum4, sum5, difmax, f1, f2, f3, f4, df1, df2, dif;
    int divaxn = 0, divaxo = 0;
    while(true) {
        /** Begin basic rule **/
        double rgnvol = (double)two_to_the_n;
        for(int j = 0; j < n; j++) {
            rgnvol *= width[j];
            z[j] = center[j];
        }
        sum1 = apply<Function, n>(f, z);

        /* Compute symmetric sums of f(lambda2,0,...,0) and f(lambda4,0,...,0), and
         * maximum fourth difference */
        sum2 = sum3 = difmax = 0.;
        for(int j = 0; j < n; j++) {
            z[j] = center[j] - lambda2*width[j];
            f1 = apply<Function, n>(f, z);
            z[j] = center[j] + lambda2*width[j];
            f2 = apply<Function, n>(f, z);
            widthl[j] = lambda4*width[j];
            z[j] = center[j] - widthl[j];
            f3 = apply<Function, n>(f, z);
            z[j] = center[j] + widthl[j];
            f4 = apply<Function, n>(f, z);
            sum2 += f1 + f2;
            sum3 += f3 + f4;
            df1 = f1 + f2 - 2*sum1;
            df2 = f3 + f4 - 2*sum1;
            dif = fabs(df1 - ratio*df2);
            if(dif >= difmax) {
                difmax = dif;
                divaxn = j;
            }
            z[j] = center[j];
        }

        /* Compute symmetric sum of f(lambda4,lambda4,0,...,0) */
        sum4 = 0.;
        for(int j = 1; j < n; j++) {
            for(int k = j; k < n; k++) {
                for(int l = 1; l <= 2; l++) {
                    widthl[j-1] = -widthl[j-1];
                    z[j-1] = center[j-1] + widthl[j-1];
                    for(int m = 1; m <= 2; m++) {
                        widthl[k] = -widthl[k];
                        z[k] = center[k] + widthl[k];
                        sum4 += apply<Function, n>(f, z);
                    }
                }
                z[k] = center[k];
            }
            z[j-1] = center[j-1];
        }

        /* Compute symmetric sum of f(lambda5,lambda5,...,lambda5) */
        sum5 = 0.;
        for(int j = 0; j < n; j++) {
            widthl[j] = -lambda5*width[j];
            z[j] = center[j] + widthl[j];
        }
        for(int j = 0; j != n; ) {
            sum5 += apply<Function, n>(f, z);
            for(j = 0; j < n; j++) {
                widthl[j] = -widthl[j];
                z[j] = center[j] + widthl[j];
                if(widthl[j] > 0)
                    break;
            }
        }

        /* Compute fifth and seventh degree rules and error */
        double rgncmp = rgnvol*(wtp1*sum1 + wtp2*sum2 + wtp3*sum3 + wtp4*sum4);
        double rgnval = rgnvol*(wt1*sum1 + wt2*sum2 + wt3*sum3 + wt4*sum4 + wt5*sum5);
        double rgnerr = fabs(rgnval - rgncmp);

        /** (End basic rule) **/

        finest += rgnval;
        abserr += rgnerr;
        funcls += rulcls;

        /** Place results of basic rule into partially ordered list according to
         ** subregion error **/
        if(divflg == 0) {
            /* When divflg == 0 start at top of list and move down list tree to
             * find corect position for results from first half of recently
             * divided subregion */
            while(true) {
                subtmp = 2*subrgn;
                if(subtmp > sbrgns)
                    break;
                if(subtmp < sbrgns && wrkstr[subtmp].err < wrkstr[subtmp+1].err)
                    subtmp++;
                if(rgnerr >= wrkstr[subtmp].err)
                    break;

                wrkstr[subrgn] = wrkstr[subtmp];
                subrgn = subtmp;
            }
        }
        else {
            /* When divflg == 1 start at bottom right branch and move up list tree
             * to find correct position for results from second half of recently
             * divided subregion */
            while(true) {
                subtmp = subrgn/2;
                if(subtmp == 0 || rgnerr < wrkstr[subtmp].err)
                    break;
                wrkstr[subrgn] = wrkstr[subtmp];
                subrgn = subtmp;
            }
        }

        /** Store results of basic rule in correct position in list **/
        for(int j = 0; j < n; j++) {
            wrkstr[subrgn].center[j] = center[j];
            wrkstr[subrgn].width[j] = width[j];
        }
        wrkstr[subrgn].val = rgnval;
        wrkstr[subrgn].err = rgnerr;
        wrkstr[subrgn].divaxn = divaxn;

        if(divflg == 0) {
            /* When divflg == 0 prepare for second application of basic rule */
            center[divaxo] += 2*width[divaxo];
            sbrgns++;
            subrgn = sbrgns;
            divflg = 1;
            continue;
        }
        /** (End ordering and storage of basic rule results) **/

        /** Make checks for possible termination of routine **/
        double relerr = abserr/fabs(finest);
        if(sbrgns + 1 > maxrgns)
            ifail = 2;
        if(funcls + 2*rulcls > GM_MAXPTS)
            ifail = 1;
        if((relerr < epsrel || abserr < epsabs) && funcls >= GM_MINPTS)
            ifail = 0;
        if(ifail < 3)
            break;

        /** Prepare to use basic rule on each half subregion with largest
         ** error **/
        divflg = 0;
        subrgn = 1;
        abserr -= wrkstr[subrgn].err;
        finest -= wrkstr[subrgn].val;
        divaxo = wrkstr[subrgn].divaxn;
        for(int j = 0; j < n; j++) {
            center[j] = wrkstr[subrgn].center[j];
            width[j] = wrkstr[subrgn].width[j];
        }
        width[divaxo] *= 0.5;
        center[divaxo] -= width[divaxo];
    }

    workspace_manager.release_gm_workspace(workspace);

    if(ifail == 1)
        warning("Integrate: did not converge after %d function evaluations\n", funcls);
    else if(ifail == 2)
        warning("Integrate: not enough storage space (this should never happen)\n");

    if(pabserr)
        *pabserr = abserr;
    if(pneval)
        *pneval = funcls;

    return finest;
}
