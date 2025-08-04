//! \file tridgeo.cpp
#include <cmath>
#include <array>
#include <iostream>
#include "tridgeo.h"
#include "sim5lib.h"

/*!
\brief Constructor.
@param aa BH spin
@param ss corona radius
@param hminn corona minimum height along the rotation axis
@param hmaxx corona maximum height along the rotation axis
 */
zamoslab3d::zamoslab3d(const double aa, const double ss, const double hminn, const double hmaxx) {
    name = "zamoslab3d";
    a = aa;
    s = ss;
    hmin = hminn;
    hmax = hmaxx;
    
    rmax = std::sqrt(s * s + hmax * hmax);
    rmin = hmin;
    
    mumax = 1.;
    mumin = hmin / std::sqrt(s * s + hmin * hmin);
}

/*!
\brief Tell if the photon is inside the corona.
*/
bool zamoslab3d::inside(double rnow, double munow) const {
    double hnow = rnow * munow, snow = std::sqrt(rnow * rnow - hnow * hnow);
    return (hnow >= hmin && hnow <= hmax && snow <= s);
}

/*!
\brief Returns tetrad.
*/
sim5tetrad zamoslab3d::caltet(double r, double mu) const {
    sim5metric met;
    sim5tetrad tet;
    kerr_metric(a, r, mu, &met);
    tetrad_zamo(&met, &tet);
    return tet;
}

/*!
\brief Calculates four velocity.
*/
void zamoslab3d::calumu(const double r, const double mu, std::array<double, 4> & umu, sim5metric & met) const {
    fourvelocity_zamo(&met, umu.data());
}
    
/*!
\brief Constructor.
@param aa BH spin
@param ss corona radius
@param hminn corona minimum height along the rotation axis
@param hmaxx corona maximum height along the rotation axis
*/
kepslab3d::kepslab3d(const double aa, const double ss, const double hminn, const double hmaxx) {
    name = "kepslab3d";
    a = aa;
    s = ss;
    hmin = hminn;
    hmax = hmaxx;
    
    rmax = std::sqrt(s * s + hmax * hmax);
    rmin = hmin;
    
    mumax = 1.;
    mumin = hmin / std::sqrt(s * s + hmin * hmin);
}


bool kepslab3d::inside(double rnow, double munow) const {
    double hnow = rnow * munow, snow = std::sqrt(rnow * rnow - hnow * hnow);
    return (hnow >= hmin && hnow <= hmax && snow <= s);
}

sim5tetrad kepslab3d::caltet(double r, double mu) const {
    double x, Omega;
    sim5metric met;
    sim5tetrad tet;
    x = std::sqrt(r);
    Omega = 1. / (x * x * x + a);
    kerr_metric(a, r, mu, &met);
    tetrad_azimuthal(&met, Omega, &tet);
    return tet;
}

void kepslab3d::calumu(const double r, const double mu, std::array<double, 4> & umu, sim5metric & met) const {
    double x, Omega;
    x = std::sqrt(r);
    Omega = 1. / (x * x * x + a);
    fourvelocity_azimuthal(Omega, &met, umu.data());
}

/*!
\brief Constructor.
@param aa BH spin
@param hh height
@param RR radius
*/
zamosphere3d::zamosphere3d(const double aa, const double hh, const double RR) {
    name = "zamosphere3d";
    a = aa;
    h = hh;
    R = RR;
    rmax = h + R;
    mumax = 1.;
    if (h > 1e-8) {
        rmin = h - R;
        rmax = h + R;
        mumax = 1.;
        mumin = R / h;
    } else {
        rmin = 0.;
        rmax = R;
        mumin = -1.;
    }
}

bool zamosphere3d::inside(double rnow, double munow) const {
    double dissqr = rnow * rnow + h * h - 2. * rnow * h * munow;
    //dissqr = rnow * rnow * (1. - munow * munow) + (rnow * munow - h) * (rnow * munow - h);
    return (dissqr <= R * R);
}

sim5tetrad zamosphere3d::caltet(double r, double mu) const {
    sim5metric met;
    sim5tetrad tet;
    kerr_metric(a, r, mu, &met);
    tetrad_zamo(&met, &tet);
    return tet;
}

void zamosphere3d::calumu(const double r, const double mu, std::array<double, 4> & umu, sim5metric & met) const {
    if (mu < 1.) {
        fourvelocity_zamo(&met, umu.data());
    } else {
        umu[0] = std::sqrt(-1. / met.g00);
        umu[1] = 0.;
        umu[2] = 0.;
        umu[3] = 0.;
    }
}

/*! \brief Constructor.
@param aa BH spin
@param rr outer radius of the shell
@param rinn inner radius of the shell
 */
zamosphere3d_truncated::zamosphere3d_truncated(const double aa, const double rr, const double rinn) {
    name = "zamosphere3d_truncated";
    a = aa;
    R = rr;
    rin = rinn;
    
    rmax = R;
    rmin = rin;
    
    mumax = 1.;
    mumin = -1.;
}

bool zamosphere3d_truncated::inside(double rnow, double munow) const {
    //dissqr = rnow * rnow * (1. - munow * munow) + (rnow * munow - h) * (rnow * munow - h);
    return (rnow <= R && rnow >= rin);
}

sim5tetrad zamosphere3d_truncated::caltet(double r, double mu) const {
    sim5metric met;
    sim5tetrad tet;
    kerr_metric(a, r, mu, &met);
    tetrad_zamo(&met, &tet);
    return tet;
}

void zamosphere3d_truncated::calumu(const double r, const double mu, std::array<double, 4> & umu, sim5metric & met) const {
    if (mu < 1.) {
        fourvelocity_zamo(&met, umu.data());
    } else {
        umu[0] = std::sqrt(-1. / met.g00);
        umu[1] = 0.;
        umu[2] = 0.;
        umu[3] = 0.;
    }
}

bool rotsphere3d::inside(double rnow, double munow) const {
    double dissqr = rnow * rnow + h * h - 2. * rnow * h * munow;
    //dissqr = rnow * rnow * (1. - munow * munow) + (rnow * munow - h) * (rnow * munow - h);
    return (dissqr <= R * R);
}

/*!
\brief Object for corona rotating with given linear velocity
@param aa BH spin
@param hh height
@param RR radius
@param vv linear circular velocity 
*/
rotsphere3d::rotsphere3d(const double aa, const double hh, const double RR, const double vv) {
    name = "rotsphere3d";
    a = aa;
    h = hh;
    R = RR;
    velocity = vv;
    gamma = 1. / std::sqrt(1. - velocity * velocity);
    rmax = h + R;
    mumax = 1.;
    if (h > 1e-8) {
        rmin = h - R;
        rmax = h + R;
        mumax = 1.;
        mumin = R / h;
    } else {
        rmin = 0.;
        rmax = R;
        mumin = -1.;
    }
}

sim5tetrad rotsphere3d::caltet(double r, double mu) const {
    double omega;
    sim5metric met;
    sim5tetrad zamotet, tet;
    kerr_metric(a, r, mu, &met);
    std::array<double, 4> umu_on = {gamma, 0., 0., gamma * velocity}, umu_bl;
    if (mu < 1.) {
        tetrad_zamo(&met, &zamotet);
        on2bl(umu_on.data(), umu_bl.data(), &zamotet);
        omega = umu_bl[3] / umu_bl[0];
        tetrad_azimuthal(&met, omega, &tet);
    } else {
        tetrad_zamo(&met, &tet);
    }
    return tet;
}

void rotsphere3d::calumu(const double r, const double mu, std::array<double, 4> & umu, sim5metric & met) const {
    sim5tetrad zamotet;
    std::array<double, 4> umu_on = {gamma, 0., 0., gamma * velocity};
    if (mu < 1.) {
        tetrad_zamo(&met, &zamotet);
        on2bl(umu_on.data(), umu.data(), &zamotet);
    } else {
        fourvelocity_zamo(&met, umu.data());
    }
}

bool sphere3d_sr::inside(const std::array<double, 3> & xmu) const {
    double rnowsqr;
    rnowsqr = xmu[0] * xmu[0] + xmu[1] * xmu[1] + xmu[2] * xmu[2];
    return (rnowsqr <= radius * radius);
}


bool slab3d_sr::inside(const std::array<double, 3> & xmu) const {
    double snowsqr = xmu[0] * xmu[0] + xmu[1] * xmu[1];
    return (std::abs(xmu[2]) <= h && snowsqr <= s * s);
}
