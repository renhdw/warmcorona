//! \file diskobj.cpp

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include "diskobj.h"
#include "const.h"
#include "kerr.h"
#include "sim5lib.h"

//! Stefan-Boltzman constant, in [\f$\rm erg~cm^{-2}~s^{-1}~K^{-4}\f$].
const double SB_SIGMA = 5.670400e-05;

//! print information to std::ostream object `out`.
//! @param out output stream. Can be either `stdout` or file stream.
//! @param disk `disk` object.
std::ostream & operator << (std::ostream & out, const diskobj & disk) {
    out << "====================================" << std::endl;
    out << disk.name << std::endl;
    out << "a = " << disk.a << std::endl;
    out << "M = " << disk.m << " Msun" << std::endl;
    out << "Mdot = " << disk.mdot << " Mdot_Edd" << std::endl;
    out << "rms = " << disk.rms << " GM/c^2" << std::endl;
    out << "====================================" << std::endl;
    return out;
}

//! Saves \f$a,\ m, \ mdot\f$, pre-calculate rms, and a few variables to be used in `flux()`, `teff()`, and `teff_kev()`.
//! Following Eq. 14 of Page & Thorne 1974 (hereafter PT74; http://adsabs.harvard.edu/abs/1974ApJ...191..499P).
//! @param aa black hole spin
//! @param mm black hole mass in [Msun]
//! @param mmdot mass accretion rate in [Eddington rate], where the Eddington rate is defined as
//!     \f$2.225475942\times10^{18}\ (M/M_\odot)\ \rm g\ s^{-1}\f$.
nt::nt(const double aa, const double mm, const double mmdot) {
    name = "Novikov-Thorne Disc";
    a = aa;
    m = mm;
    mdot = mmdot;

    // calculate rms;
    double z1,z2,r0, temp;
    double sga = (a >= 0.0) ? +1. : -1.;
    
    if (a < 1.) {
        z1 = 1. + std::pow(1. -a * a, 1. / 3.) * (std::pow(1. + a, 1. / 3.) + std::pow(1. - a, 1. / 3.));
        z2 = std::sqrt(3. * a * a + z1 * z1);
        r0 = 3. + z2 - sga * std::sqrt((3. - z1) * (3. + z1 + 2. * z2));
        rms = r0 + 1e-3;
        x0 = std::sqrt(rms);
        temp = 1. / 3. * std::acos(a);
    } else {
        z1 = 1.;
        z2 = 2.;
        r0 = 1.;
        rms = 1. + 1e-3;
        x0 = 1.;
        temp = 0.;
    }

    x1 = +2. * std::cos(temp - M_PI / 3.);
    x2 = +2. * std::cos(temp + M_PI / 3.);
    x3 = -2. * std::cos(temp);
    // some useful variables for calculating flux/Teff;
    f0base = -1. * x0;

    if (a < 1.) {
        f1base = 3. * (x1 - a) * (x1 - a) / (x1 * (x1 - x2) * (x1 - x3));
        f2base = 3. * (x2 - a) * (x2 - a) / (x2 * (x2 - x1) * (x2 - x3));
        f3base = 3. * (x3 - a) * (x3 - a) / (x3 * (x3 - x1) * (x3 - x2));
    } else {
        f1base = 0.;
        f2base = 0.;
        f3base = -1.5;
    }
}

//! @param r radius in [\f$\rm GM~c^{-2}\f$]
double nt::flux(const double r) const {
    double x, f0, f1, f2, f3, F, flux;
    x = std::sqrt(r);
    f0 = x + f0base - 1.5 * a * std::log(x / x0);
    if (a < 1.) {
        f1 = f1base * std::log((x - x1)/(x0 - x1));
        f2 = f2base * std::log((x - x2)/(x0 - x2));
    } else {
        f1 = f2 = 0.;
    }
    f3 = f3base * std::log((x - x3)/(x0 - x3));
    
    F = 1./(4.*M_PI*r) * 1.5/(x*x*(x*x*x-3.*x+2.*a)) * (f0-f1-f2-f3);
    flux = 9.1721376255e+28 * F * mdot / m;
    return flux;
}


//! @param r radius in \f$[\rm GM~c^{-2}]\f$
double nt::teff(const double r) const {
    double f, tem;
    f = flux(r);
    tem = std::pow(f / SB_SIGMA, 0.25);
    return tem;
}

//! @param r radius in \f$[\rm GM~c^{-2}]\f$
double nt::teff_kev(const double r) const {
    double f, tem;
    f = flux(r);
    tem = std::pow(f / SB_SIGMA, 0.25);
    return tem * K_KEV;
}

//! @param r radius in \f$[\rm GM~c^{-2}]\f$
//! @param lambda \f$\equiv L_z/E_\infty\f$, where \f$L_z\f$ is the angular momentum of 
//!     photon about the symmetry axis and is one constant of motion in Kerr spacetime
double nt::gfactor(const double r, const double lambda) const {
    double x = std::sqrt(r);
    double xB, xC, Omega;
    double redshift;

    xB = 1. + a / (x * x * x);
    xC = 1. - 3./(x*x) + 2.*a/(x*x*x);
    Omega = 1. / x / x / x / xB;

    redshift = sqrt(xC) / xB / (1. - Omega * lambda);
    return redshift;
}

//! Reference: Eq. B37 in Li+09.
//! @param r radius in \f$[\rm GM~c^{-2}]\f$
double nt::vphi(const double r) const {
  double A, omega, Omega, delta, vphi;
  A = r * r * r * r + a * a * r * (r + 2.);
  delta = r * r - 2. * r + a * a;
  omega = 2. * a * r / A;
  Omega = 1. / (a + r * std::sqrt(r));
  vphi  = (Omega - omega) / r / r / std::sqrt(delta) * A;
  return vphi;
}

//! Reference: Equations C20 and C23 of Li+2005 (http://adsabs.harvard.edu/abs/2005ApJS..157..335L)
//! @param r radius in \f$[\rm GM~c^{-2}]\f$.
//! @param l \f$\equiv L_z/E_\infty\f$, where \f$L_z\f$ is the photon angular momentum about the symmetry axis,
//!     and \f$E_\infty\f$ is the photon energy at infinity
//! @param q \f$\equiv \mathcal{Q} / E_\infty\f$, where \f$\mathcal{Q}\f$ is Carter's constant
//! @param signkr sign of \f$k_r\f$
//! @param cosang output; array of cosine and sin of photon emission angles \f$\theta\f$ and \f$\phi\f$
//!     (polar and azimuthal angles, respectively): \f$[{\rm cos}\theta, {\rm sin}\theta, {\rm cos}\phi, {\rm sin}\phi]\f$;
void nt::angles(const double r, const double l, const double q, const double signkr, std::array<double, 4> & cosang) {
    double delta, g, varc, varcsqrtm1, varf, varg, x, xm1, xm2, xm3, xm4, spece, specl;
    x = std::sqrt(r);
    xm1 = 1. / x;
    xm2 = xm1 * xm1;
    xm3 = xm1 * xm2;
    xm4 = xm1 * xm3;
    delta = r * r - 2. * r + a * a;
    g = gfactor(r, l);
    varc = 1. - 3. * xm2 + 2. * a * xm3;
    varcsqrtm1 = 1. / std::sqrt(varc);
    varf = 1. - 2. * a * xm3 + a * a * xm4;
    varg = 1. - 2. * xm2 + a * xm3;
    spece = varg * varcsqrtm1;
    specl = x * varf * varcsqrtm1;
    
    cosang[0] = std::sqrt(q) * g / r;
    cosang[1] = std::sqrt(1. - cosang[0] * cosang[0]);
    cosang[3] = g * (l * spece - specl) / std::sqrt(delta) / cosang[1];
    cosang[2] = signkr * std::sqrt(1. - cosang[3] * cosang[3]);
}

//! This member function assumes that the polarization angle is \f$\pi/2\f$.
//! Reference: 
//! @param r radius in \f$[\rm GM~c^{-2}]\f$
//! @param l see `angles()`.
//! @param q see `angles()`.
std::complex<double> nt::calwp(const double r, const double l, const double q) {
    std::array<double, 4> cosang;
    double cost, sint, cosp, sinp;
    double vphival, Gamma, Delta, A0, X, Y, Eloc, Knorm, K1, K2;
    
    angles(r, l, q, 1., cosang);
    cost = cosang[0];
    sint = cosang[1];
    cosp = cosang[2];
    sinp = cosang[3];
    
    vphival = vphi(r);
    Gamma = 1. / std::sqrt(1. - vphival * vphival);
    Delta = r * r - 2. * r + a * a;
    A0 = r * r * r * r + a * a * r * (r + 2.);
    X = r * r + a * a - a * std::sqrt(Delta) * vphival;
    Y = -1. * a * std::sqrt(Delta) + (r * r + a * a) * vphival;
    Eloc = 1. / gfactor(r, l);
    Knorm = -1. * Eloc * Gamma * r / std::sqrt(A0);
    K1 = Knorm * (sinp * X + sint * Y);
    K2 = Knorm * cost * cosp * X;
    
    std::complex<double> wp(K1, K2);
    return wp;
}

//! We take the polarization angle to be \f$\pi/2\f$.
//! Reference: B46-B47 of Li+2009 (http://adsabs.harvard.edu/abs/2009ApJ...691..847L)
//! @param r radius in \f$[\rm GM~c^{-2}]\f$
//! @param theta polar emission angle
//! @param phi azimuthal emission angle
//! @param l \f$\equiv L_z/E_\infty\f$, where \f$L_z\f$ is the photon angular momentum about the symmetry axis,
//!     and \f$E_\infty\f$ is the photon energy at infinity
std::complex<double> nt::calwp_angle(const double r, const double theta, const double phi, const double l) {
    double vphival, Gamma, Delta, A0, X, Y, Eloc, Knorm, K1, K2;
    double cost, sint, cosp, sinp;
    
    cost = std::cos(theta);
    sint = std::sin(theta);
    cosp = std::cos(phi);
    sinp = std::sin(phi);
    
    vphival = vphi(r);
    Gamma = 1. / std::sqrt(1. - vphival * vphival);
    Delta = r * r - 2. * r + a * a;
    A0 = r * r * r * r + a * a * r * (r + 2.);
    X = r * r + a * a - a * std::sqrt(Delta) * vphival;
    Y = -1. * a * std::sqrt(Delta) + (r * r + a * a) * vphival;
    Eloc = 1. / gfactor(r, l);
    Knorm = -1. * Eloc * Gamma * r / std::sqrt(A0);
    K1 = Knorm * (sinp * X + sint * Y);
    K2 = Knorm * cost * cosp * X;
    
    std::complex<double> wp(K1, K2);
    return wp;
}

/*!
@param r radius in \f$[\rm GM~c^{-2}]\f$.
@param theta polar emission angle
@param phi azimuthal emission angle
@param l \f$\equiv L_z/E_\infty\f$, where \f$L_z\f$ is the photon angular momentum about the symmetry axis, and \f$E_\infty\f$ is the photon energy at infinity
@param polang polarization angle
*/
std::complex<double> nt::calwp_angle(const double r, const double theta, const double phi, const double l, const double polang) {
    double vphival, Gamma, Delta, A0, X, Y, Eloc, Knorm, K1, K2, c1, c2;
    double cost, sint, cosp, sinp, ftheta, fphi;
    
    cost = std::cos(theta);
    sint = std::sin(theta);
    cosp = std::cos(phi);
    sinp = std::sin(phi);
    
    ftheta = std::cos(polang);
    fphi = std::sin(polang);
    
    vphival = vphi(r);
    Gamma = 1. / std::sqrt(1. - vphival * vphival);
    Delta = r * r - 2. * r + a * a;
    A0 = r * r * r * r + a * a * r * (r + 2.);
    X = r * r + a * a - a * std::sqrt(Delta) * vphival;
    Y = -1. * a * std::sqrt(Delta) + (r * r + a * a) * vphival;
    Eloc = 1. / gfactor(r, l);
    Knorm = Eloc * Gamma * r / std::sqrt(A0);
    c1 = cost * cosp * X;
    c2 = (sinp * X + sint * Y);
    K1 = Knorm * (c1 * ftheta - c2 * fphi);
    K2 = -1. * Knorm * (c2 * ftheta + c1 * fphi);
    
    std::complex<double> wp(K1, K2);
    return wp;
}

/*! 
The equation for Walker-Penrose constant evaluated on the disc plane (Equations B46-B47 of Li+2009)
can be re-written into the following form of vector rotation:
 \f[\begin{bmatrix} &K_1^\prime& \\ &-K_2^\prime& \end{bmatrix}
    =
    \begin{bmatrix} c_1 & -c_2 \\ c_2 & c_1 \end{bmatrix} \begin{bmatrix} f^\theta \\
    f^\phi  \end{bmatrix},
 \f] where
 \f{eqnarray}
  K_1^\prime &=& \frac{K_1}{\sqrt{K_1^2 + K_2^2}}, \\
  K_2^\prime &=& \frac{K_2}{\sqrt{K_1^2 + K_2^2}}, \\
  c_1 &=& \frac{E_{\rm loc} \Gamma r {\rm cos}\theta {\rm cos}\phi X}{\sqrt{A_0(K_1^2 + K_2^2)}}, \\
  c_2 &=& \frac{E_{\rm loc} \Gamma r ({\rm sin}\phi X + {\rm sin}\theta Y)}{\sqrt{A_0(K_1^2 + K_2^2)}}.
 \f}
and \f$K_1\f$ and \f$K_2\f$ are the real and imaginary parts of the Walker-Penrose constant, respectively.
This member function calculates the coefficients of the rotation matrix \f$c_1\f$ and \f$c_2\f$.
@param r radius in \f$[\rm GM~c^{-2}]\f$
@param theta polar emission angle
@param phi azimuthal emission angle
@param l \f$\equiv L_z/E_\infty\f$, where \f$L_z\f$ is the photon angular momentum about the symmetry axis,
    and \f$E_\infty\f$ is the photon energy at infinity
@param q \f$\equiv \mathcal{Q} / E_\infty\f$, where \f$\mathcal{Q}\f$ is Carter's constant
@param c1 coefficients of the rotation matrix; see above
@param c2 coefficients of the rotation matrix; see above
*/
void nt::polrotang_disktok(const double r, const double theta, const double phi, const double l, const double q, double & c1, double & c2) {
    double vphival, Gamma, Delta, A0, X, Y, Eloc;
    double cost, sint, cosp, sinp;
    double cnorm;
    
    cost = std::cos(theta);
    sint = std::sin(theta);
    cosp = std::cos(phi);
    sinp = std::sin(phi);
    
    vphival = vphi(r);
    Gamma = 1. / std::sqrt(1. - vphival * vphival);
    Delta = r * r - 2. * r + a * a;
    A0 = r * r * r * r + a * a * r * (r + 2.);
    X = r * r + a * a - a * std::sqrt(Delta) * vphival;
    Y = -1. * a * std::sqrt(Delta) + (r * r + a * a) * vphival;
    Eloc = 1. / gfactor(r, l);
    cnorm = 1. / std::sqrt(A0 * (q + (l - a) * (l - a)));
    
    c1 = Eloc * Gamma * r * cost * cosp * X * cnorm;
    c2 = Eloc * Gamma * r * (sinp * X + sint * Y) * cnorm;
}

/*!
Same as @ref polrotang_disktok(const double, const double, const double, const double, const double, double &, double &),
    but return the rotation angle instead of the coefficients. For detail see the referred function.
@param r radius in \f$[\rm GM~c^{-2}]\f$
@param theta polar emission angle
@param phi azimuthal emission angle
@param l \f$\equiv L_z/E_\infty\f$
@param q \f$\equiv \mathcal{Q} / E_\infty\f$
*/
double nt::polrotang_disktok(const double r, const double theta, const double phi, const double l, const double q) {
    double c1, c2, ang;
    polrotang_disktok(r, theta, phi, l, q, c1, c2);
    
    ang = std::acos(c1);
    if (c2 < 0.)
        ang = 2. * M_PI - ang;
    return ang;
}

/*!
Similarly with @ref polrotang_disktok(const double, const double, const double, const double, const double, double &, double &),
    but for rotation from \f$[K_1^\prime, -K_2^\prime]\f$ to the polarisation vector at infinity.
@param ktheta sign of \f$k^\theta\f$
@param l \f$\equiv L_z/E_\infty\f$
@param q \f$\equiv \mathcal{Q} / E_\infty\f$
@param muobs \f$\mu_{\rm obs}\f$, cosine of polar emission angle of the observer at infinity
@param cospsi cosine of the rotation angle
@param sinpsi sinusoidal of the rotation angle
*/
void nt::polrotang_ktoinf(const double ktheta, const double l, const double q, const double muobs, double & cospsi, double & sinpsi) {
    double signktheta, sintheta, S, T, s2plust2, knorm;
    double sqrT, signs;
    
    signktheta = ktheta / std::abs(ktheta);
    sintheta = std::sqrt(1. - muobs * muobs);
    
    S = l / sintheta - a * sintheta;
    sqrT = q + a * a * muobs * muobs - l * l * muobs * muobs / sintheta / sintheta;
    if (sqrT < 0.) {// && sqrT >= -1e-6)
        sqrT = 0.;
        signs = (S >= 0.) ? 1. : -1.; 
        S = signs * std::sqrt(q + (l - a) * (l - a));
    }
    T = signktheta * std::sqrt(sqrT);
    s2plust2 = q + (l - a) * (l - a);
    knorm = -1. / std::sqrt(s2plust2);
    cospsi = T * knorm;
    sinpsi = S * knorm;
}

/*!
Saves \f$a,\ m, \ mdot\f$, pre-calculate rms, and a few variables to be used in `flux()`, `teff()`, and `teff_kev()`.
Following Eq. 14 of Page & Thorne 1974 (hereafter PT74; http://adsabs.harvard.edu/abs/1974ApJ...191..499P).
@param aa black hole spin
@param mm black hole mass in [Msun]
@param mmdot mass accretion rate in [Eddington rate], where the Eddington rate is defined as
    \f$2.225475942\times10^{18}\ (M/M_\odot)\ \rm g\ s^{-1}\f$.
@param rtrr truncation radius in \f$[\rm GM~c^{-2}]\f$
*/
truncated_nt::truncated_nt(const double aa, const double mm, const double mmdot, const double rtrr) {
    name = "Truncated Novikov-Thorne Disc";
    a = aa;
    m = mm;
    mdot = mmdot;
    rtr = rtrr;
    x0 = std::sqrt(rtr);
    
    double z1, z2, r0, temp;
    double sga = (a >= 0.0) ? +1. : -1.;
    
    if (a < 1.) {
        z1 = 1. + std::pow(1. -a * a, 1. / 3.) * (std::pow(1. + a, 1. / 3.) + std::pow(1. - a, 1. / 3.));
        z2 = std::sqrt(3. * a * a + z1 * z1);
        r0 = 3. + z2 - sga * std::sqrt((3. - z1) * (3. + z1 + 2. * z2));
        rms = r0 + 1e-3;
        temp = 1. / 3. * std::acos(a);
    } else {
        z1 = 1.;
        z2 = 2.;
        r0 = 1.;
        rms = 1. + 1e-3;
        temp = 0.;
    }

    x1 = +2. * std::cos(temp - M_PI / 3.);
    x2 = +2. * std::cos(temp + M_PI / 3.);
    x3 = -2. * std::cos(temp);
    // some useful variables for calculating flux/Teff;
    f0base = -1. * x0;

    if (a < 1.) {
        f1base = 3. * (x1 - a) * (x1 - a) / (x1 * (x1 - x2) * (x1 - x3));
        f2base = 3. * (x2 - a) * (x2 - a) / (x2 * (x2 - x1) * (x2 - x3));
        f3base = 3. * (x3 - a) * (x3 - a) / (x3 * (x3 - x1) * (x3 - x2));
    } else {
        f1base = 0.;
        f2base = 0.;
        f3base = -1.5;
    }
}

void truncated_nt::selftest(const double r) {
    double x, f0, f1, f2, f3, f;
    
    x = std::sqrt(r);
    std::cout << "x0 = " << x0 << std::endl;
    f0 = x + f0base - 1.5 * a * std::log(x / x0);
    if (a < 1.) {
        f1 = f1base * std::log((x - x1)/(x0 - x1));
        f2 = f2base * std::log((x - x2)/(x0 - x2));
    } else {
        f1 = f2 = 0.;
    }
    f3 = f3base * std::log((x - x3)/(x0 - x3));
    
    f = f0 - f1 - f2 - f3;
    
    std::cout << "From Page & Thorne 1974 Equation, f = " << f << std::endl;
    
    double integauss = 0.;
    const long NGAUSS1 = 21;
    double weights[NGAUSS1], abissca[NGAUSS1];
    
    gauleg(rtr, r, abissca, weights, NGAUSS1);
    for (int i = 0; i < NGAUSS1; ++i) {
        integauss += weights[i] * nt_ffunc(abissca[i], a);
    }
    std::cout << "By Gauss-Legendre method, f = " << integauss << std::endl;
}
