#include <cmath>

template<typename Function>
vector<real> RungeKutta4(real t0, real t1, real y0, Function dydt, int N) {
    /* Initialize output array */
    vector<real> Y(N+1);
    Y[0] = y0;

    real t;
    real dt = (t1 - t0)/N;
    real k1, k2, k3, k4;

    /* Evolve N steps of Runge-Kutta 4 */
    for(int i = 0; i < N; i++) {
        t = t0 + i*dt;
        k1 = dydt(t, Y[i]);
        k2 = dydt(t + 0.5*dt, Y[i] + 0.5*dt*k1);
        k3 = dydt(t + 0.5*dt, Y[i] + 0.5*dt*k2);
        k4 = dydt(t + dt, Y[i] + dt*k3);
        Y[i+1] = Y[i] + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    }

    return Y;
}

template<typename Function>
vector<vector<real> > RungeKutta4(real t0, real t1, const vector<real>& y0, Function dydt, int N) {
    int m = (int)y0.size();

    /* Initialize output array */
    vector<vector<real> > Y(m, vector<real>(N+1));
    for(int j = 0; j < m; j++)
        Y[j][0] = y0[j];

    real t;
    real dt = (t1 - t0)/N;
    vector<real> k1, k2, k3, k4;
    vector<real> y(m);

    /* Evolve N steps of Runge-Kutta 4 */
    for(int i = 0; i < N; i++) {
        t = t0 + i*dt;
        for(int j = 0; j < m; j++)
            y[j] = Y[j][i];
        k1 = dydt(t, y);
        for(int j = 0; j < m; j++)
            y[j] = Y[j][i] + 0.5*dt*k1[j];
        k2 = dydt(t + 0.5*dt, y);
        for(int j = 0; j < m; j++)
            y[j] = Y[j][i] + 0.5*dt*k2[j];
        k3 = dydt(t + 0.5*dt, y);
        for(int j = 0; j < m; j++)
            y[j] = Y[j][i] + dt*k3[j];
        k4 = dydt(t + dt, y);
        for(int j = 0; j < m; j++)
            Y[j][i+1] = Y[j][i] + (dt/6)*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
    }

    return Y;
}

namespace DormandPrince {
    /* Butcher tableau values */
    const real a1 = 1/5., a2 = 3/10., a3 = 4/5., a4 = 8/9., a5 = 1., a6 = 1.;
    const real b11 = 1/5.,
               b21 = 3/40.,       b22 = 9/40.,
               b31 = 44/45.,      b32 = -56/15.,      b33 = 32/9.,
               b41 = 19372/6561., b42 = -25360/2187., b43 = 64448/6561., b44 = -212/729.,
               b51 = 9017/3168.,  b52 = -355/33.,     b53 = 46732/5247., b54 = 49/176.,  b55 = -5103/18656.,
               b61 = 35/384.,     b62 = 0.,           b63 = 500/1113.,   b64 = 125/192., b65 = -2187/6784., b66 = 11/84.;
    const real c1 = 35/384., c2 = 0., c3 = 500/1113., c4 = 125/192., c5 = -2187/6784., c6 = 11/84., c7 = 0.;
    const real d1 = 5179/57600., d2 = 0., d3 = 7571/16695., d4 = 393/640., d5 = -92097/339200., d6 = 187/2100., d7 = 1/40.;
    const real e1 = c1-d1, e2 = c2-d2, e3 = c3-d3, e4 = c4-d4, e5 = c5-d5, e6 = c6-d6, e7 = c7-d7;
}

template<typename Function>
int RKDP(real t0, real t1, real y0, Function dydt, vector<real>& t, vector<real>& y, real epsrel, real h0) {
    using namespace DormandPrince;

    /* Set initial conditions */
    t.resize(1);
    y.resize(1);
    t[0] = t0;
    y[0] = y0;

    const int maxsteps = 1000000;
    const real hmin = 1e-7*(t1 - t0);

    real k1, k2, k3, k4, k5, k6, k7;
    real ytmp, yerr, ytol, r;
    real h = (h0 > 0) ? h0 : (t1 - t0);
    int n = 0;          // step number
    int neval = 1;      // number of derivative evaluations
    k1 = dydt(t[0], y[0]);
    while(t[n] < t1) {
        if(h < hmin) {
            warning("RKDP: step size too small: h = %g (after %d steps)\n", h, n+1);
            return -1;
        }
        if(n > maxsteps) {
            warning("RKDP: too many steps: n = %d\n", n);
            return -1;
        }

        k2 = dydt(t[n] + h*a1, y[n] + h*b11*k1);
        k3 = dydt(t[n] + h*a2, y[n] + h*(b21*k1 + b22*k2));
        k4 = dydt(t[n] + h*a3, y[n] + h*(b31*k1 + b32*k2 + b33*k3));
        k5 = dydt(t[n] + h*a4, y[n] + h*(b41*k1 + b42*k2 + b43*k3 + b44*k4));
        k6 = dydt(t[n] + h*a5, y[n] + h*(b51*k1 + b52*k2 + b53*k3 + b54*k4 + b55*k5));
        ytmp = y[n] + h*(b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5 + b66*k6);
        k7 = dydt(t[n] + h*a6, ytmp);
        neval += 6;

        yerr = h*(e1*k1 + e2*k2 + e3*k3 + e4*k4 + e5*k5 + e6*k6 + e7*k7);
        ytol = epsrel * (fabs(y[n]) + h*fabs(k1) + 1e-20);
        r = fabs(yerr/ytol);
        if(r > 1) {
            /* Error too large: try smaller h */
            h *= fmax(0.01, 0.8*pow(r, -0.25));
        }
        else {
            /* Error within tolerance: update result vectors and h */
            t.push_back(t[n] + h);
            y.push_back(ytmp);
            n++;

            k1 = k7;
            h *= fmin(10, 0.8*pow(r, -0.2));
            if(t[n] + h > t1) { // adjust h if the next step would take us past t1
                h = t1 - t[n];
                if(h <= hmin) { // if the new h is very small (e.g. because we did this on the last step) then break
                    t[n] = t1;
                    break;
                }
            }
        }
    }

    return n+1;
}

template<typename Function>
int RKDP(real t0, real t1, const vector<real>& y0, Function dydt, vector<real>& y, real epsrel, real h0) {
    using namespace DormandPrince;
    int m = (int)y0.size();
    if(m > (int)y.size())
        y.resize(m);

    /* Set initial conditions */
    real t = t0;
    for(int j = 0; j < m; j++)
        y[j] = y0[j];

    const int maxsteps = 1000000;
    const real hmin = 1e-7*(t1 - t0);

    vector<real> k1(m), k2(m), k3(m), k4(m), k5(m), k6(m), k7(m);
    vector<real> ytmp = y0;
    real yerr, ytol, r;
    real h = (h0 > 0) ? h0 : (t1 - t0);
    int n = 0;          // step number
    int neval = 1;      // number of derivative evaluations
    k1 = dydt(t0, ytmp);
    while(t < t1) {
        if(h < hmin) {
            warning("RKDP: step size too small: h = %g (after %d steps)\n", h, n);
            return -1;
        }
        if(n > maxsteps) {
            warning("RKDP: too many steps: n = %d\n", n);
            return -1;
        }

        for(int j = 0; j < m; j++)
            ytmp[j] = y[j] + h*b11*k1[j];
        k2 = dydt(t + h*a1, ytmp);
        for(int j = 0; j < m; j++)
            ytmp[j] = y[j] + h*(b21*k1[j] + b22*k2[j]);
        k3 = dydt(t + h*a2, ytmp);
        for(int j = 0; j < m; j++)
            ytmp[j] = y[j] + h*(b31*k1[j] + b32*k2[j] + b33*k3[j]);
        k4 = dydt(t + h*a3, ytmp);
        for(int j = 0; j < m; j++)
            ytmp[j] = y[j] + h*(b41*k1[j] + b42*k2[j] + b43*k3[j] + b44*k4[j]);
        k5 = dydt(t + h*a4, ytmp);
        for(int j = 0; j < m; j++)
            ytmp[j] = y[j] + h*(b51*k1[j] + b52*k2[j] + b53*k3[j] + b54*k4[j] + b55*k5[j]);
        k6 = dydt(t + h*a5, ytmp);
        for(int j = 0; j < m; j++)
            ytmp[j] = y[j] + h*(b61*k1[j] + b62*k2[j] + b63*k3[j] + b64*k4[j] + b65*k5[j] + b66*k6[j]);
        k7 = dydt(t + h*a6, ytmp);
        neval += 6;

        r = 0;
        for(int j = 0; j < m; j++) {
            yerr = h*(e1*k1[j] + e2*k2[j] + e3*k3[j] + e4*k4[j] + e5*k5[j] + e6*k6[j] + e7*k7[j]);
            ytol = epsrel * (fabs(y[j]) + h*fabs(k1[j]) + 1e-20);
            r = fmax(r, fabs(yerr/ytol));
        }
        if(r > 1) {
            /* Error too large: try smaller h */
            h *= fmax(0.1, 0.8*pow(r, -0.25));
        }
        else {
            /* Error within tolerance: update result vectors and h */
            t += h;
            for(int j = 0; j < m; j++)
                y[j] = ytmp[j];
            n++;

            k1 = k7;
            h *= fmin(10, 0.8*pow(r, -0.2));
            if(t + h > t1) {    // adjust h if the next step would take us past t1
                h = t1 - t;
                if(h <= hmin) { // if the new h is very small (e.g. because we did this on the last step) then break
                    t = t1;
                    break;
                }
            }
        }
    }

    return n+1;
}
