//************************************************************************
//    SIM5 library
//    sim5math.h - math macros
//------------------------------------------------------------------------
//    Author:
//    Michal Bursa (bursa@astro.cas.cz)
//    Astronomical Institute
//    Bocni II 1401/1, 141-00 Praha 4, Czech Republic
//************************************************************************


#ifndef _SIM5MATH_H
#define _SIM5MATH_H

#ifdef __cplusplus
	extern "C"{
#endif


//#include "sim5config.h"


//#ifdef SIM5_SINGLE_PREC_MATH
//    #define double   float
//#else
//    #define double   double
//#endif


#define PI      3.14159265359
#define PI2     6.28318530718
#define PI4     12.5663706144
#define PI_half 1.57079632679

#define sqr(a)   ((a) * (a))
#define sqr2(a)  ((a) * (a))
#define sqr3(a)  ((a) * (a) * (a))
#define sqr4(a)  ((a) * (a) * (a) * (a))
#define sqrt3(a) cbrt(a)
#define sqrt4(a) pow(a,0.25)
#define sim5max(a,b) ((a) > (b) ? (a) : (b))
#define sim5min(a,b) ((a) < (b) ? (a) : (b))
#define sim5minmax(val, vmin, vmax) min(vmax,max(val,vmin))
//#define max(a,b) ((a) > (b) ? (a) : (b))
//#define min(a,b) ((a) < (b) ? (a) : (b))
//#define minmax(val, vmin, vmax) min(vmax,max(val,vmin))
#define odd(a) ((a%2==1)?1:0)
#define sign(a) ((a) >= 0.0 ? (+1.0) : (-1.0))
//#define hypot(a,b) (sqrt((a)*(a) + (b)*(b)))
#define deg2rad(a) ((a)/180.0*M_PI)
#define rad2deg(a) ((a)*180.0/M_PI)
#define EE(a) pow(10.0,a)
#define ave(a, b, w) ((1.0-(w))*(a) + (w)*(b))
#define logave(a, b, w) (exp((1.0-(w))*log(a) + (w)*log(b)))
#define inrange(a, min, max) (((a)>=(min))&&((a)<=(max)))

//#define rand  sim5rand()
#define urand sim5urand()




#ifdef CUDA
    typedef double2     sim5complex;
    #define ComplexI    makeComplex(0.0,1.0)
#else
    #include <complex.h>
    #undef I
    #define ComplexI _Complex_I
    //typedef double complex sim5complex;
    typedef double _Complex sim5complex;
#endif


DEVICEFUNC INLINE long sim5round(double num);

DEVICEFUNC INLINE long int factorial(long int n);

DEVICEFUNC INLINE double reduce_angle_pi(double phi);
DEVICEFUNC INLINE double reduce_angle_2pi(double phi);

DEVICEFUNC int ensure_range(double *val, double min, double max, double acc);


DEVICEFUNC void cartesian2spherical1(double x, double y, double z, double Vx, double Vy, double Vz, double* Vr, double* Vh, double* Vf);
DEVICEFUNC void cartesian2spherical2(double cos_h, double sin_f, double cos_f, double Vx, double Vy, double Vz, double* Vr, double* Vh, double* Vf);


DEVICEFUNC INLINE void sim5seed();
DEVICEFUNC INLINE unsigned long long sim5rand();
DEVICEFUNC INLINE double sim5urand();

DEVICEFUNC INLINE
double creal (sim5complex a);


DEVICEFUNC INLINE
double cimag (sim5complex a);






DEVICEFUNC INLINE sim5complex makeComplex(double r, double i);
DEVICEFUNC INLINE sim5complex nullComplex();

#ifdef __cplusplus
	}
#endif

#endif
