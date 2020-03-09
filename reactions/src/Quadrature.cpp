#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <map>

#include "Common.h"
#include "Quadrature.h"


/******************************************************************************
 * DiscreteIntegrate
 ******************************************************************************/

double DiscreteIntegrate(int n, double* f, double h) {
    /* Use Simpson's rule when we have an even number of intervals (odd number of points) */
    if((n % 2) == 1) {
        double S = f[0] + f[n-1];
        for(int i = 1; i < n-1; i++)
            S += 2*(1 + (i%2))*f[i];
        return S*h/3;
    }
    /* Otherwise use Hollingsworth and Hunter's 3rd-order formula */
    else if(n == 2)
        return h/2 * (f[0] + f[1]);
    else if(n == 4)
        return h/8 * (3*f[0] + 9*f[1] + 9*f[2] + 3*f[3]);
    else if(n == 6)
        return h/24 * (9*f[0] + 28*f[1] + 23*f[2] + 23*f[3] + 28*f[4] + 9*f[5]);
    else {
        double S = ( 9*(f[0] + f[n-1]) + 28*(f[1] + f[n-2]) + 23*(f[2] + f[n-3]) )/24;
        for(int i = 3; i < n-3; i++)
            S += f[i];
        return S*h;
    }
}


/******************************************************************************
 * GKWorkspace
 ******************************************************************************/

GKWorkspace::GKWorkspace(int n) : limit(n) {
    assert(n > 0);
    size = 0;
    alist = (double*)malloc(n*sizeof(double));
    blist = (double*)malloc(n*sizeof(double));
    rlist = (double*)malloc(n*sizeof(double));
    elist = (double*)malloc(n*sizeof(double));
    order = (int*)malloc(n*sizeof(int));
}

GKWorkspace::~GKWorkspace() {
    free(alist);
    free(blist);
    free(rlist);
    free(elist);
    free(order);
}

void GKWorkspace::GetMaxInterval(double& a, double& b, double& r, double& e) {
    const int i_max = order[0];
    a = alist[i_max];
    b = blist[i_max];
    r = rlist[i_max];
    e = elist[i_max];
}

void GKWorkspace::Update(double a1, double b1, double area1, double error1,
                          double a2, double b2, double area2, double error2)
{
    const int i_max = order[0];
    const int i_new = size;

    /* We replace the interval i_max with the larger (meaning the one with the
     * larger error) of the two newly created subintervals.  The smaller
     * interval is appended at the end of the list.  This speeds up the
     * sort routine somewhat. */

    /* Append the newly-created intervals to the list */
    if(error2 > error1) {
        alist[i_max] = a2;      // blist[i_max] is already == b2
        rlist[i_max] = area2;
        elist[i_max] = error2;

        alist[i_new] = a1;
        blist[i_new] = b1;
        rlist[i_new] = area1;
        elist[i_new] = error1;
    }
    else {
        blist[i_max] = b1;      // alist[i_max] is already == a1
        rlist[i_max] = area1;
        elist[i_max] = error1;

        alist[i_new] = a2;
        blist[i_new] = b2;
        rlist[i_new] = area2;
        elist[i_new] = error2;
    }
    size++;

    Sort();
}

void GKWorkspace::Sort() {
    const int last = size - 1;

    /* Check whether the list contains more than two error estimates */
    if(last < 2) {
        order[0] = 0;
        order[1] = 1;
        return;
    }

    int i_max = order[0];
    double errmax = elist[i_max];

    /* Compute the number of elements in the list to be maintained in
     * descending order.  This number depends on the number of subdivisions
     * still allowed. */
    int top = (last < limit/2 + 2) ? last : limit - last + 1;

    /* Insert errmax by traversing the list top-down */
    int j = 1;
    while(j < top && errmax < elist[order[j]]) {
        order[j-1] = order[j];
        j++;
    }
    order[j-1] = i_max;

    /* Insert errmin by traversing the list bottom-up */
    double errmin = elist[last];
    int k = top - 1;
    while(k > j-2 && errmin >= elist[order[k]]) {
        order[k+1] = order[k];
        k--;
    }
    order[k+1] = last;
}

void GKWorkspace::SumResults(double& result, double& abserr) {
    result = abserr = 0;
    for(int i = 0; i < size; i++) {
        result += rlist[i];
        abserr += elist[i];
    }
}

void GKWorkspace::ReportError(int error_code) {
    switch(error_code) {
    case 0:
        break;
    case 1:
        warning("Integrate: tolerance cannot be achieved\n");
        break;
    case 2:
        warning("Integrate: roundoff error prevents tolerance from being achieved\n");
        break;
    case 3:
        warning("Integrate: bad integrand behavior found in the integration interval\n");
        break;
    case 4:
        warning("Integrate: maximum number of subdivisions reached\n");
        break;
    default:
        warning("Integrate: could not integrate function\n");
        break;
    }
}

/* Weights for 15-point Gauss-Kronrod rule */
const double GKWorkspace::wg15[4] = {
    0.1294849661688697e+00,
    0.2797053914892767e+00,
    0.3818300505051189e+00,
    0.4179591836734694e+00
};
const double GKWorkspace::wgk15[8] = {
    0.2293532201052922e-01,
    0.6309209262997855e-01,
    0.1047900103222502e+00,
    0.1406532597155259e+00,
    0.1690047266392679e+00,
    0.1903505780647854e+00,
    0.2044329400752989e+00,
    0.2094821410847278e+00
};
const double GKWorkspace::xgk15[8] = {
    0.9914553711208126e+00,
    0.9491079123427585e+00,
    0.8648644233597691e+00,
    0.7415311855993944e+00,
    0.5860872354676911e+00,
    0.4058451513773972e+00,
    0.2077849550078985e+00,
    0.0e+00
};

/* General Gauss-Kronrod integration */
double GKWorkspace::GaussKronrod(const int n, const double wg[], const double wgk[],
                                 double fc, double fv1[], double fv2[], double half,
                                 double* abserr, double* resabs, double* resasc)
{
    double result_gauss = (n % 2 == 0) ? fc * wg[n/2-1] : 0;
    double result_kronrod = fc * wgk[n-1];
    double result_abs = fabs(result_kronrod);
    double result_asc = 0;
    double mean = 0, err = 0;

    double f1 = 0, f2 = 0;
    for(int j = 0; j < (n-1)/2; j++) {
        f1 = fv1[2*j+1];
        f2 = fv2[2*j+1];
        result_gauss += wg[j] * (f1 + f2);
        result_kronrod += wgk[2*j+1] * (f1 + f2);
        result_abs += wgk[2*j+1] * (fabs(f1) + fabs(f2));
    }
    for(int j = 0; j < n/2; j++) {
        f1 = fv1[2*j];
        f2 = fv2[2*j];
        result_kronrod += wgk[2*j] * (f1 + f2);
        result_abs += wgk[2*j] * (fabs(f1) + fabs(f2));
    }

    mean = 0.5 * result_kronrod;

    result_asc = wgk[n-1] * fabs(fc - mean);
    for(int j = 0; j < n-1; j++)
        result_asc += wgk[j] * (fabs(fv1[j]-mean) + fabs(fv2[j]-mean));

    /* Scale by the width of the integration region */
    err = fabs(result_kronrod - result_gauss) * half;
    result_kronrod *= half;
    result_abs *= fabs(half);
    result_asc *= fabs(half);

    if(abserr) {
        if(result_asc != 0 && err != 0) {
            double scale = pow((200*err/result_asc), 1.5);
            if(scale < 1)
                err = result_asc * scale;
            else
                err = result_asc;
        }
        if(result_abs > DBL_MIN/(50*DBL_EPSILON))
            err = fmax(err, 50*DBL_EPSILON*result_abs);
        *abserr = err;
    }
    if(resabs)
        *resabs = result_abs;
    if(resasc)
        *resasc = result_asc;
    return result_kronrod;
}


/******************************************************************************
 * GKWorkspace
 ******************************************************************************/

GMWorkspace::GMWorkspace(int n_) : n(n_) {
    wrkstr = malloc(GM_MAXPTS*((2*n+2)*sizeof(double) + sizeof(int)));
}

GMWorkspace::~GMWorkspace() {
    free(wrkstr);
}


/******************************************************************************
 * WorkspaceManager
 ******************************************************************************/

GKWorkspace* WorkspaceManager::get_gk_workspace() {
    GKWorkspace* ws = NULL;

    /* Guard against race condition where multiple OpenMP parallel regions
     * request a workspace simultaneously. */
    #pragma omp critical
    {
        /* Look for a pre-existing workspace */
        for(GKWorkspaceList::iterator it = gk_workspaces.begin(); it != gk_workspaces.end(); it++) {
            if(it->second == 0) {
                it->second = 1;
                ws = it->first;
                break;
            }
        }
        if(ws == NULL) {
            /* Otherwise create a new one */
            ws = new GKWorkspace();
            gk_workspaces[ws] = 1;
        }
    }

    return ws;
}

void WorkspaceManager::release_gk_workspace(GKWorkspace* ws) {
    if(gk_workspaces.find(ws) != gk_workspaces.end())
        gk_workspaces[ws] = 0;
}

GMWorkspace* WorkspaceManager::get_gm_workspace(int n) {
    GMWorkspace* ws = NULL;

    /* Guard against race condition where multiple OpenMP parallel regions
     * request a workspace simultaneously. */
    #pragma omp critical
    {
        /* Look for a pre-existing workspace */
        for(GMWorkspaceList::iterator it = gm_workspaces.begin(); it != gm_workspaces.end(); it++) {
            if(it->first->n == n && it->second == 0) {
                it->second = 1;
                ws = it->first;
                break;
            }
        }
        if(ws == NULL) {
            /* Otherwise create a new one */
            ws = new GMWorkspace(n);
            gm_workspaces[ws] = 1;
        }
    }

    return ws;
}

void WorkspaceManager::release_gm_workspace(GMWorkspace* ws) {
    if(gm_workspaces.find(ws) != gm_workspaces.end())
        gm_workspaces[ws] = 0;
}

WorkspaceManager::WorkspaceManager() {
}

WorkspaceManager::~WorkspaceManager() {
    /* Free allocated workspaces */
    for(GKWorkspaceList::iterator it = gk_workspaces.begin(); it != gk_workspaces.end(); it++)
        delete it->first;
    for(GMWorkspaceList::iterator it = gm_workspaces.begin(); it != gm_workspaces.end(); it++)
        delete it->first;
}

/* Global WorkspaceManager singleton */
WorkspaceManager workspace_manager;
