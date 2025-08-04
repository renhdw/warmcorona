//************************************************************************
//    SIM5 library
//    sim5raytrace.c - step-wise raytracing
//------------------------------------------------------------------------
//    Author:
//    Michal Bursa (bursa@astro.cas.cz)
//    Astronomical Institute of the Czech Academy of Sciences
//************************************************************************


//! \file sim5raytrace.c
//! Raytracing.
//! 
//! Provides routines for step-wise raytracing.



/*
#include "sim5config.h"
#ifndef CUDA
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sim5config.h"
#include "sim5math.h"
#include "sim5kerr.h"
#include "sim5polyroots.h"
#include "sim5photontrace.h"
#endif
*/


//! \cond SKIP
//! (skip doc for these)
#define frac_error(a,b) (fabs((b)-(a))/(fabs(b)+1e-40))
#define vect_copy(v1,v2) {v2[0]=v1[0]; v2[1]=v1[1]; v2[2]=v1[2]; v2[3]=v1[3];}
//! \endcond

//! default maximal relative error of raytracing
#define raytrace_max_error 1e-2


DEVICEFUNC
void raytrace_prepare(double bh_spin, double x[4], double k[4], double f[4], double presision_factor, int options, raytrace_data* rtd)
//! Raytracing with step-wise null-geodesic integration.
//! Makes one step along the geodesic that is tangent to `k` and updates input vectors with new values.
//! The integration method follows Dolence+2009 (http://adsabs.harvard.edu/abs/2009ApJS..184..387D).
//! 
//! The routine automatically controls size of the step based on curvature and required precission. 
//! On each call to raytrace() the routine takes a step size, which is smaller of the two: the internally 
//! chosen step size and step size that is passed on input in `step`.
//!
//! Numerical precision of integration is driven by the precision_factor modifier and raytrace_max_error constant; 
//! roughly, final error (after whole geodesic is integrated) is: 
//! maximal relative error  = (a factor of few) * raytrace_max_error * precision_factor.
//!
//! @param bh_spin black hole spin
//! @param x initial position vector
//! @param k initial direction vector (photon 4-momentum)
//! @param f initial polarization vector (optional, can be NULL)
//! @param precision_factor precision factor
//! @param options additional options
//! @param rtd raytracing data
{
    // read options
    rtd->opt_gr  = !((options & RTOPT_FLAT) == RTOPT_FLAT);
    rtd->opt_pol =  ((options & RTOPT_POLARIZATION) == RTOPT_POLARIZATION);
    rtd->step_epsilon = sqrt(presision_factor)/10.;   // note: precision ~ (step_epsilon)^2; step_epsilon=0.1 gives reasonable precision ~1e-3

    if (rtd->opt_pol) {
        #ifndef CUDA
        fprintf(stderr,"ERR (raytrace_prepare): the polarization vector transport does not work sufficiently well yet. The photon_polarization_vector() routine should be used instead.\n");
        #endif
    }

    // evaluate metric and connection 
    sim5metric m;
    double G[4][4][4];
    rtd->opt_gr ? kerr_metric(bh_spin, x[1], x[2], &m) : flat_metric(x[1], x[2], &m);
    rtd->opt_gr ? kerr_connection(bh_spin, x[1], x[2], G) : flat_connection(x[1], x[2], G);

    // check that k.k=0
    double kk = dotprod(k, k, &m);
    #ifndef CUDA
    if (fabs(kk) > 1e-10) fprintf(stderr,"ERR (kerr_raytrace_prepare): k is not a null vector (k.k=%.3e)\n", kk);
    #endif

    if (rtd->opt_pol) {
        // check that k.f=0
        double kf = dotprod(k, f, &m);
        #ifndef CUDA
        if (fabs(kf) > 1e-10) fprintf(stderr,"ERR (kerr_raytrace_prepare): k and f are not orthogonal (k.f=%.3e)\n", kf);
        #endif
    }

    // set motion constants
    rtd->bh_spin = bh_spin;
    rtd->E = k[0]*m.g00 + k[3]*m.g03;
    rtd->Q  = photon_carter_const(k, &m);
    if (rtd->opt_pol) rtd->WP = photon_wp_const(k, f, &m);

    // set runtime varibles    
    rtd->pass = 0;
    rtd->refines = 0;
    rtd->kt = rtd->E;
    rtd->error = 0.0;

    // evaluate intial derivatives of k and f
    Gamma(G, k, k, rtd->dk);
    if (rtd->opt_pol) Gamma(G, k, f, rtd->df);
}

//! \cond SKIP
#ifdef CUDA
DEVICEFUNC
inline double k_deriv(int j, double _k[4], double G[4][4][4]) {
    int a,b;
    double _dki = 0.0;
    for (a=0;a<4;a++) for (b=a;b<4;b++) _dki -= G[j][a][b]*_k[a]*_k[b];
    return _dki;
}

DEVICEFUNC
inline double f_deriv(int j, double _k[4], double _f[4], double G[4][4][4]) {
    int a,b;
    double _dfi = 0.0;
    for (a=0;a<4;a++) for (b=a;b<4;b++) _dfi -= 0.5*G[j][a][b]*(_k[a]*_f[b] + _k[b]*_f[a]);
    return _dfi;
}
#endif
//! \endcond


DEVICEFUNC
void raytrace(double x[4], double k[4], double f[4], double *step, raytrace_data* rtd)
//! Raytracing with step-wise null-geodesic integration.
//! Makes one step along the geodesic that is tangent to `k` and updates input vectors with new values.
//! The integration method follows Dolence+2009 (http://adsabs.harvard.edu/abs/2009ApJS..184..387D).
//! 
//! The routine automatically controls size of the step based on curvature and required precission. 
//! On each call to raytrace() the routine takes a step size, which is smaller of the two: the internally 
//! chosen step size and step size that is passed on input in `step`.
//!
//! Numerical precision is driven by the precision_factor modifier [see raytrace_prepare()];
//! rtd->error gives the error in the current step and it should be bellow raytrace_max_error*1e-2;
//! so it is adviceable to check rtd->error continuously after each step and stop integration 
//! when the error goes above ~1e-3; at the end of integration one should then check 
//! relative difference in Carter's constant with `raytrace_error()`.
//!
//! @param x position vector
//! @param k direction vector (photon 4-momentum)
//! @param f polarization vector (optional, can be NULL)
//! @param step on input gives maximal step that the routine can take [GM/c2]; on output gives size of step that has actually been taken
//! @param rtd raytracing data
//!
//! @result Position, direction and polarization (if not null) vectors and step size are updated, `rtd` has the relative error.
{
    int i;
    //int iter = 0;
    sim5metric m;
    double G[4][4][4];
    double* dk = rtd->dk;
    double* df = rtd->df;
    double x_orig[4], k_orig[4], f_orig[4]={0,0,0,0};
    double xp[4], kp[4], kp_prev[4], fp[4], fp_prev[4];
    double kk=0.0, kt = rtd->kt;
    float k_frac_error;

    vect_copy(x, x_orig);
    vect_copy(k, k_orig);
    if (rtd->opt_pol) vect_copy(f, f_orig);

    #ifndef CUDA    
    // evaluate derivative of momentum from the geodesic equation; 
    // - this function does the same thing as Gamma() from sim5kerr.c (see there for info), 
    //   except this form takes out the summation over the upper index and uses the symmetry
    //   (also the inline form it runs ~2 times faster - no need for passing G to the function)
    // - for CUDA it cannot be nested, so it is taked out of raytrace()
    inline double k_deriv(int j, double _k[4]) {
        int a,b;
        double _dki = 0.0;
        for (a=0;a<4;a++) for (b=a;b<4;b++) _dki -= G[j][a][b]*_k[a]*_k[b];
        return _dki;
    }

    // evaluate derivative of polarization vector from the geodesic equation
    // - same as for the routine above applies here, except the symmetry assumption
    // - for CUDA it cannot be nested, so it is taked out of raytrace()
    inline double f_deriv(int j, double _k[4], double _f[4]) {
        int a,b;
        double _dfi = 0.0;
        for (a=0;a<4;a++) for (b=a;b<4;b++) _dfi -= 0.5*G[j][a][b]*(_k[a]*_f[b] + _k[b]*_f[a]);
        return _dfi;
    }
    #endif

    // limit step into a reasonable iterval
    // smaller step must be used close to black hole (radial term: 0.05*x[1])
    // smaller step must be used close to polar axis (axial term: pow(1-x[2],0.5))
    //double dl = min(0.05*x[1]*pow(1.-x[2],0.5), *step) * dl_reduction_factor;
    //if (dl < 1e-4) dl = 1e-4;
    double stepsize = rtd->step_epsilon/(fabs(dk[0])/(fabs(k[0])+TINY) + fabs(dk[1])/(fabs(k[1])+TINY) + fabs(dk[2])/(fabs(k[2])+TINY) + fabs(dk[3])/(fabs(k[3])+TINY) + TINY);
    double dl = sim5min(*step, stepsize);
    if (dl < 1e-3) dl = 1e-3;

    rtd->pass++;
    
    // step 1: update of position (Dolence+09, Eq.14a)
    // - remember that x[2] contains cos(theta); the substitution method 
    //   with dm/d\theta=-sin(theta) factor introduces unnecessary inacuracy
    xp[0] = x[0] + k[0]*dl + 0.5*dk[0]*dl*dl;
    xp[1] = x[1] + k[1]*dl + 0.5*dk[1]*dl*dl;
    xp[2] = cos(acos(x[2]) +(k[2]*dl + 0.5*dk[2]*dl*dl));
    xp[3] = x[3] + k[3]*dl + 0.5*dk[3]*dl*dl;

    // this is addition of first two terms of Eq.14d (Dolence+09) 
    if (rtd->opt_pol) {
        for (i=0;i<4;i++) {
            k[i] += 0.5*dk[i]*dl;
            f[i] += 0.5*df[i]*dl;
        }
    } else {
        for (i=0;i<4;i++) k[i] += 0.5*dk[i]*dl;
    }; 

    // update metric and connection
    if (rtd->opt_gr) {
        kerr_metric(rtd->bh_spin, xp[1], xp[2], &m);
        kerr_connection(rtd->bh_spin, xp[1], xp[2], G);
    } else {
        flat_metric(xp[1], xp[2], &m);
        flat_connection(xp[1], xp[2], G);
    }; 

    // step 2: estimate new value for k and f (Dolence+09, Eq.14b)
    if (rtd->opt_pol) {
        for (i=0;i<4;i++) {
            kp[i] = k[i] + 0.5*dk[i]*dl;
            fp[i] = f[i] + 0.5*df[i]*dl;
        }
    } else {
        for (i=0;i<4;i++) kp[i] = k[i] + 0.5*dk[i]*dl;
    }


    // step 3: iteratively improve estimate of momentum based on its derivative
    int k_iter = 0;
    do {
        k_frac_error = 0.0;

        if (rtd->opt_pol) {
            vect_copy(kp, kp_prev);
            vect_copy(fp, fp_prev);
            for (i=0;i<4;i++) {
                // Dolence+09, Eq. (14c,d)
                #ifdef CUDA
                kp[i] = k[i] + 0.5*k_deriv(i,kp_prev,G)*dl;         
                fp[i] = f[i] + 0.5*f_deriv(i,kp_prev,fp_prev,G)*dl;
                #else
                kp[i] = k[i] + 0.5*k_deriv(i,kp_prev)*dl;         
                fp[i] = f[i] + 0.5*f_deriv(i,kp_prev,fp_prev)*dl;
                #endif
                k_frac_error += frac_error(kp[i],kp_prev[i]);
            }
        } else {
            vect_copy(kp, kp_prev);
            for (i=0;i<4;i++) {
                // Dolence+09, Eq. (14c,d)
                #ifdef CUDA
                kp[i] = k[i] + 0.5*k_deriv(i,kp_prev,G)*dl;
                #else
                kp[i] = k[i] + 0.5*k_deriv(i,kp_prev)*dl;
                #endif
                k_frac_error += frac_error(kp[i],kp_prev[i]);
            }
        }

        k_iter++;
	} while (k_frac_error>raytrace_max_error*1e-3 && k_iter<3);


    // precision check
    kt = kp[0]*m.g00 + kp[3]*m.g03;
    kk = fabs(dotprod(kp, kp, &m));
    rtd->error = sim5max(frac_error(kt,rtd->kt), kk);
	if ((k_frac_error>raytrace_max_error*1e-2) || (rtd->error>raytrace_max_error*1e-2)) {
        vect_copy(x_orig, x);
        vect_copy(k_orig, k);
        if (rtd->opt_pol) vect_copy(f_orig, f);
        DEVICEFUNC void raytrace_rk4(double x[4], double k[4], double f[4], double dl, raytrace_data* rtd);
        raytrace_rk4(x, k, f, dl, rtd);
        *step = dl;
        return;
    }

    // assign x, k, f, and dk, df
    if (rtd->opt_pol) {
        for (i=0;i<4;i++) {
            x[i]  = xp[i];
            k[i]  = kp[i];
            f[i]  = fp[i];
            #ifdef CUDA
            dk[i] = k_deriv(i,kp,G);
            df[i] = f_deriv(i,kp,fp,G);
            #else
            dk[i] = k_deriv(i,kp);
            df[i] = f_deriv(i,kp,fp);
            #endif
        }
    } else {
        for (i=0;i<4;i++) {
            x[i]  = xp[i];
            k[i]  = kp[i];
            #ifdef CUDA
            dk[i] = k_deriv(i,kp,G);
            #else
            dk[i] = k_deriv(i,kp);
            #endif
        }
    }

    // assign motion constant for k in the current step
    rtd->kt = kt;

    // assign the actual size of step taken
    *step = dl;
}

DEVICEFUNC
void raytrace1(double x[4], double k[4], double *step, raytrace_data* rtd)
//! Raytracing with step-wise null-geodesic integration.
//! Makes one step along the geodesic that is tangent to `k` and updates input vectors with new values.
//! The integration method follows Dolence+2009 (http://adsabs.harvard.edu/abs/2009ApJS..184..387D).
//! 
//! The routine automatically controls size of the step based on curvature and required precission. 
//! On each call to raytrace() the routine takes a step size, which is smaller of the two: the internally 
//! chosen step size and step size that is passed on input in `step`.
//!
//! Numerical precision is driven by the precision_factor modifier [see raytrace_prepare()];
//! rtd->error gives the error in the current step and it should be bellow raytrace_max_error*1e-2;
//! so it is adviceable to check rtd->error continuously after each step and stop integration 
//! when the error goes above ~1e-3; at the end of integration one should then check 
//! relative difference in Carter's constant with `raytrace_error()`.
//!
//! @param x position vector
//! @param k direction vector (photon 4-momentum)
//! @param f polarization vector (optional, can be NULL)
//! @param step on input gives maximal step that the routine can take [GM/c2]; on output gives size of step that has actually been taken
//! @param rtd raytracing data
//!
//! @result Position, direction and polarization (if not null) vectors and step size are updated, `rtd` has the relative error.
{
    int i;
    //int iter = 0;
    sim5metric m;
    double G[4][4][4];
    double* dk = rtd->dk;
    double x_orig[4], k_orig[4];
    double xp[4], kp[4], kp_prev[4];
    double kk=0.0, kt = rtd->kt;
    float k_frac_error;

    vect_copy(x, x_orig);
    vect_copy(k, k_orig);

    #ifndef CUDA    
    // evaluate derivative of momentum from the geodesic equation; 
    // - this function does the same thing as Gamma() from sim5kerr.c (see there for info), 
    //   except this form takes out the summation over the upper index and uses the symmetry
    //   (also the inline form it runs ~2 times faster - no need for passing G to the function)
    // - for CUDA it cannot be nested, so it is taked out of raytrace()
    inline double k_deriv(int j, double _k[4]) {
        int a,b;
        double _dki = 0.0;
        for (a=0;a<4;a++) for (b=a;b<4;b++) _dki -= G[j][a][b]*_k[a]*_k[b];
        return _dki;
    }

    // evaluate derivative of polarization vector from the geodesic equation
    // - same as for the routine above applies here, except the symmetry assumption
    // - for CUDA it cannot be nested, so it is taked out of raytrace()
    inline double f_deriv(int j, double _k[4], double _f[4]) {
        int a,b;
        double _dfi = 0.0;
        for (a=0;a<4;a++) for (b=a;b<4;b++) _dfi -= 0.5*G[j][a][b]*(_k[a]*_f[b] + _k[b]*_f[a]);
        return _dfi;
    }
    #endif

    // limit step into a reasonable iterval
    // smaller step must be used close to black hole (radial term: 0.05*x[1])
    // smaller step must be used close to polar axis (axial term: pow(1-x[2],0.5))
    //double dl = min(0.05*x[1]*pow(1.-x[2],0.5), *step) * dl_reduction_factor;
    //if (dl < 1e-4) dl = 1e-4;
    double stepsize = rtd->step_epsilon/(fabs(dk[0])/(fabs(k[0])+TINY) + fabs(dk[1])/(fabs(k[1])+TINY) + fabs(dk[2])/(fabs(k[2])+TINY) + fabs(dk[3])/(fabs(k[3])+TINY) + TINY);
    double dl = sim5min(*step, stepsize);
    if (dl < 1e-3) dl = 1e-3;

	double halfdl = 0.5 * dl, halfdl2 = 0.5 * dl * dl;

    rtd->pass++;
    
    // step 1: update of position (Dolence+09, Eq.14a)
    // - remember that x[2] contains cos(theta); the substitution method 
    //   with dm/d\theta=-sin(theta) factor introduces unnecessary inacuracy
    xp[0] = x[0] + k[0]*dl + dk[0] * halfdl2;
    xp[1] = x[1] + k[1]*dl + dk[1] * halfdl2;
    xp[2] = cos(acos(x[2]) +(k[2]*dl + dk[2] * halfdl2));
    xp[3] = x[3] + k[3]*dl + dk[3] * halfdl2;

	for (i=0;i<4;i++) k[i] += dk[i] * halfdl;

    // update metric and connection
    if (rtd->opt_gr) {
        kerr_metric(rtd->bh_spin, xp[1], xp[2], &m);
        kerr_connection1(rtd->bh_spin, xp[1], xp[2], G);
    } else {
        flat_metric(xp[1], xp[2], &m);
        flat_connection(xp[1], xp[2], G);
    }; 

    // step 2: estimate new value for k and f (Dolence+09, Eq.14b)
	for (i=0;i<4;i++) kp[i] = k[i] + dk[i] * halfdl;

    // step 3: iteratively improve estimate of momentum based on its derivative
    int k_iter = 0;
    do {
        k_frac_error = 0.0;

        vect_copy(kp, kp_prev);
        for (i=0;i<4;i++) {
            // Dolence+09, Eq. (14c,d)
            #ifdef CUDA
            kp[i] = k[i] + 0.5*k_deriv(i,kp_prev,G)*dl;
            #else
            kp[i] = k[i] + 0.5*k_deriv(i,kp_prev)*dl;
            #endif
            k_frac_error += frac_error(kp[i],kp_prev[i]);
        }

        k_iter++;
	} while (k_frac_error>raytrace_max_error*1e-3 && k_iter<3);


    // precision check
    kt = kp[0]*m.g00 + kp[3]*m.g03;
    kk = fabs(dotprod(kp, kp, &m));
    rtd->error = sim5max(frac_error(kt,rtd->kt), kk);
	if ((k_frac_error>raytrace_max_error*1e-2) || (rtd->error>raytrace_max_error*1e-2)) {
        vect_copy(x_orig, x);
        vect_copy(k_orig, k);
		double f[4] = {0, 0, 0, 0};
        DEVICEFUNC void raytrace_rk4(double x[4], double k[4], double f[4], double dl, raytrace_data* rtd);
        raytrace_rk4(x, k, f, dl, rtd);
        *step = dl;
        return;
    }

// assign x, k, f, and dk, df
    for (i=0;i<4;i++) {
        x[i]  = xp[i];
        k[i]  = kp[i];
        #ifdef CUDA
        dk[i] = k_deriv(i,kp,G);
        #else
        dk[i] = k_deriv(i,kp);
        #endif
    }

    // assign motion constant for k in the current step
    rtd->kt = kt;

    // assign the actual size of step taken
    *step = dl;
}

//! \cond SKIP
DEVICEFUNC
void raytrace_rk4(double x[4], double k[4], double f[4], double dl, raytrace_data* rtd)
{
    int i;
    sim5metric m;
    double G[4][4][4];
    double xp[4];
    double k1[4], dk1[4], k2[4], dk2[4], k3[4], dk3[4], k4[4], dk4[4];
    double f1[4], df1[4], f2[4], df2[4], f3[4], df3[4], f4[4], df4[4];
	double dl_2 = 0.5 * dl;

    double kt0 = rtd->kt;

    // transform theta-component of coordinate vector to angle
    x[2] = acos(x[2]);

	for (i=0; i<4; i++) xp[i] = x[i];
    rtd->opt_gr ? kerr_connection(rtd->bh_spin, xp[1], cos(xp[2]), G) : flat_connection(xp[1], cos(xp[2]), G);
    if (rtd->opt_pol) {
    	for (i=0; i<4; i++) k1[i] = k[i];
    	for (i=0; i<4; i++) f1[i] = f[i];
        Gamma(G, k1, k1, dk1);
        Gamma(G, k1, f1, df1);
    } else {
    	for (i=0; i<4; i++) k1[i] = k[i];
        Gamma(G, k1, k1, dk1);
    }

	for (i=0; i<4; i++) xp[i] = x[i] + k1[i]*dl_2;
    rtd->opt_gr ? kerr_connection(rtd->bh_spin, xp[1], cos(xp[2]), G) : flat_connection(xp[1], cos(xp[2]), G);
    if (rtd->opt_pol) {
    	for (i=0; i<4; i++) k2[i] = k[i] + dk1[i]*dl_2;
	    for (i=0; i<4; i++) f2[i] = f[i] + df1[i]*dl_2;
        Gamma(G, k2, k2, dk2);
        Gamma(G, k2, f2, df2);
    } else {
    	for (i=0; i<4; i++) k2[i] = k[i] + dk1[i]*dl_2;
        Gamma(G, k2, k2, dk2);
    }

	for (i=0; i<4; i++) xp[i] = x[i] + k2[i]*dl_2;
    rtd->opt_gr ? kerr_connection(rtd->bh_spin, xp[1], cos(xp[2]), G) : flat_connection(xp[1], cos(xp[2]), G);
    if (rtd->opt_pol) {
    	for (i=0; i<4; i++) k3[i] = k[i] + dk2[i]*dl_2;
    	for (i=0; i<4; i++) f3[i] = f[i] + df2[i]*dl_2;
        Gamma(G, k3, k3, dk3);
        Gamma(G, k3, f3, df3);
    } else {
    	for (i=0; i<4; i++) k3[i] = k[i] + dk2[i]*dl_2;
        Gamma(G, k3, k3, dk3);
    }

	for (i=0; i<4; i++) xp[i] = x[i] + k3[i]*dl;
    rtd->opt_gr ? kerr_connection(rtd->bh_spin, xp[1], cos(xp[2]), G) : flat_connection(xp[1], cos(xp[2]), G);
    if (rtd->opt_pol) {
    	for (i=0; i<4; i++) k4[i] = k[i] + dk3[i]*dl;
    	for (i=0; i<4; i++) f4[i] = f[i] + df3[i]*dl;
        Gamma(G, k4, k4, dk4);
        Gamma(G, k4, f4, df4);
    } else {
    	for (i=0; i<4; i++) k4[i] = k[i] + dk3[i]*dl;
        Gamma(G, k4, k4, dk4);
    }

    // update values for momentum and polarization vectors
    if (rtd->opt_pol) {
    	for (i=0; i<4; i++) {
	    	x[i] += dl/6. * ( k1[i] + 2.* k2[i] + 2.* k3[i] +  k4[i]);
	    	k[i] += dl/6. * (dk1[i] + 2.*dk2[i] + 2.*dk3[i] + dk4[i]);
	    	f[i] += dl/6. * (df1[i] + 2.*df2[i] + 2.*df3[i] + df4[i]);
	    }
	} else {
    	for (i=0; i<4; i++) {
	    	x[i] += dl/6. * ( k1[i] + 2.* k2[i] + 2.* k3[i] +  k4[i]);
	    	k[i] += dl/6. * (dk1[i] + 2.*dk2[i] + 2.*dk3[i] + dk4[i]);
	    }
	}

    // transform theta-component of coordinate vector back to angle cosine
	x[2] = cos(x[2]);

    // update values for momentum and polarization vector derivatives
    rtd->opt_gr ? kerr_connection(rtd->bh_spin, x[1], x[2], G) : flat_connection(x[1], x[2], G);
    if (rtd->opt_pol) {
        Gamma(G, k, k, rtd->dk);
        Gamma(G, k, f, rtd->df);
    } else {
        Gamma(G, k, k, rtd->dk);
    }


    kerr_metric(rtd->bh_spin, x[1], x[2], &m);
    double kt1 = k[0]*m.g00 + k[3]*m.g03;

    rtd->error = frac_error(kt1,kt0);
    /*
	if (rtd->error>1e-4) {
	    fprintf(stderr,"WRN: RK4 rtd->error=%.3e (%d/%d) dl=%.2e ke=%.2e kp[1]=%.5e kp[2]=%.5e\n", rtd->error, rtd->pass, rtd->refines, dl, 0.0, k[1], k[2]);
	    if (dl>1e-8) {
            vect_copy(x_orig, x);
            vect_copy(k_orig, k);
            if (rtd->opt_pol) vect_copy(f_orig, f);
            //rtd->refines++;
	        //for (i=0; i<10; i++) raytrace_rk4(x, k, f, dl/10., rtd);
	    } else {
	    fprintf(stderr,"WRN: RK4 step too small rtd->error=%.3e (%d/%d) dl=%.2e ke=%.2e kp[1]=%.5e kp[2]=%.5e\n", rtd->error, rtd->pass, rtd->refines, dl, 0.0, k[1], k[2]);
	    }
	}
    */
}
//! \endcond


DEVICEFUNC
double raytrace_error(double x[4], double k[4], double f[4], raytrace_data* rtd)
//! Raytracing error.
//! Gives relative error in raytracing in terms of relative difference of Carter's constant. 
//! Useful for checking precission of integration.
//! 
//! @param x position vector
//! @param k direction vector
//! @param f polarization vector (optional, can be NULL)
//! @param rtd raytracing data
//!
//! @result Relative error in raytracing.
{
    sim5metric m;
    rtd->opt_gr ? kerr_metric(rtd->bh_spin, x[1], x[2], &m) : flat_metric(x[1], x[2], &m);
    return frac_error(rtd->Q, photon_carter_const(k,&m));
}


#undef frac_error
#undef vect_copy
#undef raytrace_max_error


