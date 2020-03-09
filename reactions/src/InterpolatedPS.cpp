#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <cstdio>

#include "Datafile.h"
#include "InterpolatedPS.h"
#include "array.h"

static Spline make_spline(bool logscale, int N, const real* k, const real* p) {
    if(logscale) {
        array logk(N);
        for(int i = 0; i < N; i++)
            logk[i] = log(k[i]);
        return CubicSpline(N, logk, p);
    }
    else
        return CubicSpline(N, k, p);
}

/* Make an heuristic guess about whether the points are linearly or
 * logarithmically spaced (favor linear spacing for ambiguous cases) */
static bool is_logscale(int N, const real* k) {
    real kmin = k[0];
    real kmax = k[N-1];
    real kmid = k[N/2];
    real kmid_lin = (kmax + kmin)/2;
    real kmid_log = sqrt(kmin*kmax);
    bool logscale = 10*(fabs(kmid - kmid_log) < fabs(kmid - kmid_lin));
    return logscale;
}


InterpolatedPS::InterpolatedPS() {
}

InterpolatedPS::InterpolatedPS(const Cosmology& C_, const char* filename, int kcol, int pcol) {
    C = &C_;
    n = C->n;
    Datafile F(filename);
    array k = F.GetColumn(kcol);
    array p = F.GetColumn(pcol);
    kmin = k.front();
    kmax = k.back();
    pmin = p.front();
    pmax = p.back();
    logscale = is_logscale(k.size(), k.data());
    interp = make_spline(logscale, k.size(), k.data(), p.data());
}

InterpolatedPS::InterpolatedPS(const Cosmology& C_, const array& k, const array& p) {
    C = &C_;
    n = C->n;
    kmin = k.front();
    kmax = k.back();
    pmin = p.front();
    pmax = p.back();
    logscale = is_logscale(k.size(), k.data());
    interp = make_spline(logscale, k.size(), k.data(), p.data());
}

InterpolatedPS::InterpolatedPS(const Cosmology& C_, int N, const real* k, const real* p) {
    C = &C_;
    n = C->n;
    kmin = k[0];
    kmax = k[N-1];
    pmin = p[0];
    pmax = p[N-1];
    logscale = is_logscale(N, k);
    interp = make_spline(logscale, N, k, p);
}

InterpolatedPS::InterpolatedPS(const Cosmology& C_, const array& k, const array& p, bool logscale_) {
    C = &C_;
    n = C->n;
    kmin = k.front();
    kmax = k.back();
    pmin = p.front();
    pmax = p.back();
    logscale = logscale_;
    interp = make_spline(logscale, k.size(), k.data(), p.data());
}

InterpolatedPS::InterpolatedPS(const Cosmology& C_, int N, const real* k, const real* p, bool logscale_) {
    C = &C_;
    n = C->n;
    kmin = k[0];
    kmax = k[N-1];
    pmin = p[0];
    pmax = p[N-1];
    logscale = logscale_;
    interp = make_spline(logscale, N, k, p);
}

real InterpolatedPS::Evaluate(real k) const {
    if(k <= kmin)
        return pmin * pow(k/kmin, n);   // interpolate to 0
    else if(k >= kmax)
        return pmax * exp(-100*pow2((k-kmax)/kmax)); // cut off small-scale power smoothly
    else
        return logscale ? interp(log(k)) : interp(k);
}

const Cosmology& InterpolatedPS::GetCosmology() const {
    return *C;
}
