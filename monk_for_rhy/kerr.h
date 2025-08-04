//! \file kerr.h
#ifndef _KERR_H
#define _KERR_H
#include <vector>
#include <functional>
#include <complex>

//! Return the four velocity of a particle on the equatorial plane orbiting the BH with Keplerian velocity.
std::vector<double> keplerumu(const double & a, const double & r);
//!Return the four momentum of a photon reaching the equatorial plane, given constants of motion.
std::vector<double> fourp_disk(const double & a, const double & r, const double & l, const double & q,
    const bool rplus, const bool thetaplus);
//! Return the wave vector of a photon located on the rotation axis;
std::vector<double> fourp_axis(const double & a, const double &h, const double & q,
    const bool rplus, const bool thetaplus);
//! Returns proper area density on the disc plane.
double ds(const double r, void *p);
//! Calcuate the proper area between \f$r_0\f$ and \f$r_1\f$ on the equatorial plane:
double calds(const double a, const double r0, const double r1);
//double calds_gsl(const double a, const double r0, const double r1);
//double calmui(const double & a, const double & l, const double & q, const double & r);
//! Calculates specific intensity of planck spectrum  
void planck_kev(const double tbb, const std::vector<double> & en, std::vector<double> & flux);
//! Calculate Chandrasekhar H function, given characteristic function.
void hfunc(const std::vector<double> & mu, std::vector<double> & res, std::function<double (const double)> phifunc);
//! The characteristic function for radiation in the meridian plane.
double phil(const double mu);
//! The characteristic function for radiation perpendicular to the meridian plane.
double phir(const double mu);
//! Calculates the polarisation degree of radiation from pure scattering disc atmosphere with infinite optical depth.
void calpoldeg(const std::vector<double> & muarr, std::vector<double> & poldeg);
//! Calculates the limb-darkening factor of radiation from pure scattering disc atmosphere with infinite optical depth.
void calchandra_darkening(const std::vector<double> & muarr, std::vector<double> & dark);
//! Evaluates the \f$r-\f$integrand in equation of motion.
double kerr_rintegrand(const double r, const std::array<double, 3> & params);
//! Evaluates the \f$\mu-\f$integrand in equation of motion.
double kerr_muintegrand(const double r, const std::array<double, 3> & params);
//! Evaluates the \f$\mu-\f$integrand in calculating the affine parameter.
double kerr_lambdaintegrand_mu(const double mu, const std::array<double, 3> & params);
//! Evaluates the \f$r-\f$integrand in calculating the affine parameter.
double kerr_lambdaintegrand_r(const double r, const std::array<double, 3> & params);
double kerr_phiintegrand_mu(const double mu, const std::array<double, 3> & params);
double kerr_phiintegrand_r(const double r, const std::array<double, 3> & params);
double kerr_dtintegrand_r_p0(const double r, const std::array<double, 3> & params);
//! Evaluates the function \f$f\f$ (Eq. 15n) in Page & Thorne 1974.
double nt_ffunc(const double r, const double a);
//! Calculates photon wave vector in Boyer-Lindquist frame, given photon emission angle measured by disc fluid.
void kmu_disk_angle(const double a, const double re, const double theta, const double phi, std::array<double, 4> & bl_k);
//void zamokmu_angle(const double a, const double re, const double mu, const double theta, const double phi, std::array<double, 4> & bl_k);

void photon_momentum1(const double a, const double r, const double mu, const double l, const double q, 
    const double krsign, const double kthetasign, std::array<double, 4> & kmu);
//! Sample energy of photons that follow Planckian distribution. 
double sample_bb(const double tbb, std::mt19937_64 & gen);
//! Weighted sample of photon energy that follow Planckian distribution. 
void sample_bb1(const double tbb, std::mt19937_64 & gen, const double ecut, const double ratio, double & en, double & weight);
//! Sample powerlaw distribution.
double sample_pl(const double gamma, const double x0, const double x1, std::mt19937_64 & gen);
//! Calculates two normalised Stokes parameters Q/I and U/I at infinity, given polarisation degree and the Walker-Penrose constant. 
void cal_pol_inf(const double muobs, const double ktheta, const double a, const double l, const double q, 
    const double poldeg, const std::complex<double> wp, double & qfac, double & ufac);
double sample_cutoffpl(const double alpha, const double emin, const double ecut, std::mt19937_64 & gen);
double sample_expo(const double lambda, std::mt19937_64 & gen);
struct ds_params {
    double a;
};
#endif
