//! \file const.h
//! Constants. Most physical constants are in cgs unit.
#ifndef _CONST_H
#define _CONST_H
#include <cmath>
//! \f$\sqrt{2}\f$
const double SQRT2 = std::sqrt(2.);
//! double EPS
const double DEPS = 1e-10;
//! Boltzmann constant in \f$[\rm keV~K]\f$; cf., NIST.
const double K_KEV = 8.6173303e-8;
//! speed of light; in \f$[\rm cm\ s^{-1}]\f$; cf., NIST
const double C = 2.997925e+10;
//! planck constant in \f$[\rm keV\ s]\f$; cf. NIST
const double H = 4.135667662e-18;
//! \f${\rm log}_{10}(1/h^3c^2)\f$
const double logc2h3 = -2. * std::log10(C) - 3. * std::log10(H);
//! \f$1/h^3c^2\f$
const double  one_over_c2h3 = std::pow(10., logc2h3);
//! gravitational constant; in \f$[\rm cm^3\ g^{-1}\ s^{-2}]\f$
const double GRAV = 6.67408e-8; 
//! solar mass in \f$[\rm g]\f$
const double MSUN = 1.98892e+33;
//! Gravitational radius of solar mass black hole in \f$[\rm cm]\f$
const double RG = 1.475e5;
//! classical electron radius, in cm; cf. NIST
const double E_R0 = 2.8179403227e-13;
//! keV to erg
const double KEV2ERG = 1.60218e-9;
//! speed of light; in \f$[\rm cm\ s^{-1}]\f$; cf., NIST
const double LSPEED = 2.9979245e+10;
//! electron mass in \f$[\rm keV]\f$; cf. NIST.
const double ME_KEV = 510.9989461;
//! Apéry's constant 
const double APERY_CONST = 1.202056903159594;
//! Thompson scattering cross section
const double SIGMA_THOMSON = 8. * M_PI * E_R0 * E_R0 / 3.;

//! constant part of normalization of spectra if following superphoton (SP) approach of Schnittman & Krolik 2013.
const double SPECNORM_SP = 4. * APERY_CONST * RG * RG / (C * C * H * H * H);

//! constants while using Gauss-Legendre method
const int NGAUSS = 21;
const unsigned long MAXITER = 1000;
const double CRITERION = 1e-12;

//! constants for calculating Chandrasekhar H functions
const double POL_ALPHA0 = 2.29767;
const double POL_A0 = 1.19736;
const double POL_A1 = 0.61733;
const double POL_Q = 0.68989;
const double POL_ALPHA1 = 1.34864;
const double POL_C = 0.87294;

// minimum weight for superphotons
const double MINWEIGHT = 1e-100;

#endif
