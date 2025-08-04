//! \file gridgeod.cpp
//! This file contains several functions that sample disc photons, then send the sampled photons to
//! routines that handle single (super)photons (such as \ref tridsca_sp in \ref scatter.cpp).

#include <algorithm>
#include <vector>
#include <cmath>
#include <iostream>
#include <random>
#include <memory>
#include "diskobj.h"
#include "sim5lib.h"
#include "utils.h"
#include "quadroots.h"
#include "tridgeo.h"
#include "geoinf.h"
#include "const.h"
#include "detector.h"
#include "kerr.h"
#include "scatter.h"

/*
#ifdef USEHEAL
//! Healpix pixelization order
const int HEALORDER = 2;
#endif
*/

// calculates redshift at lamp-post
// @param a black hole spin
// @param h lamp-post height
// @param re photon emitting radius from disc
/*
double lampg(const double & a, const double & h, const double & re) {
    nt agn(a, 1., 1.);
    double diskg = agn.gfactor(re, 0.);
    double deltal =  h * h - 2. * h + a * a;
    double ut = std::sqrt((h * h + a * a) / deltal);
    return ut * diskg;
}
*/


/*
// Workhorse for calculating spectrum for on-axis lamp-post. Basiclly, it pixelate the \f$4\pi\f$ sky in the local frame of
// the lamp-post. Given the symmetry with respective to the azimuthal angle \f$\phi\f$,
// we just need to sample the poloidal angle \f$\theta\f$ and assign the appropiate solid angle
// factor.
//
// For each value of \f$\theta\f$, we solve for \f$r_e\f$, the radius where the photon hits the
// disc on the equatorial plane, under the help of functions in `lpray.cpp`. This goes as follows:
// first we calculate Carter's constant \f$\mathcal{Q}\f$, given spin, height, and emission angle,
// following Eq. A84 of Dovciak+2004. For any geodesic crossing the lamp-post on z-axis, the photon's
// angular momentum with respect to z is apparently 0, which greatly simplifies the calculation.
// One consequence is that the lamp-post is located at the turning point along \f$\theta-\f$ integration.
// So we just integrate \f$\mu={\rm cos}(\theta)\f$ from 0 to 1 to get \f$\tau_\mu\f$.
//
// Finally with \f$r_e\f$ we can calculate the spectrum in the rest frame of lamp-post. Note that the
// specific intensity \f$I_\nu\f$ is conserved in the following way (MTW):
// \f{equation}{g^3 I_\nu \equiv const.\f}
// Hence the observed flux \f$F_\nu\f$ is
// \f{equation}{F_E = 2E^3 \int \frac{{\rm cos}\theta d\Omega}{e^{E/gkT} - 1},\f}
// where \f$E\f$ is the photon energy, \f$\theta\f$ is the incident angle of the photon,
// \f$d\Omega = {\rm sin}\theta d\theta d\phi\f$, all of the above
// are measured in rest frame of the lamp-post. Then \f$g\f$ is the redshift factor, and \f$T_{\rm eff}\f$
// is the effective temperature of the NT disc (measured in rest frame of the disc), and just a
// function of \f$r_e\f$ given \f$m, mdot, a\f$. Given symmetry about \f$\phi\f$, we can rewrite the
// equation above as
// \f{equation}{F_E = 4\pi E^3 \int \frac{{\rm sin}\theta{\rm cos}\theta d\theta}{e^{E/gkT} - 1}.\f}
// @param en energy array
// @param flux flux array for output
// @param ne size of en
// @param a dimensionless spin
// @param m black hole mass in Msun
// @param mdot mass accretion rate in Eddington rate
// @param h height of lamp-post in \f$R_g\equiv GM/c^3\f$
void onaxisgrid(const long & ntheta, diskobj & disk, const double & a, const double & height, const double & vr,
    std::vector<double> & rearr, std::vector<double> & garr, std::vector<double> & domegaarr, bool todisk) {

    rearr.clear();
    garr.clear();
    domegaarr.clear();

    double vrrest;
    // disk related variables
    double gfac, eo;
    double q, re; // integral; emitting radius
    // if kr > 0 or not;
    bool plus;
    // integrals
    double inte, sint;
    int ncase;

    // vector for un-stationary case
    std::vector<double> umu(4), er(4), kmu(4), k_mu(4), etheta(4);  // four-velocity of fluid, r-tetrad, and contravariant wave vector
    // test if I was doing write
    sim5metric met;
    kerr_metric(a, height, 1., &met);

    // calculate vr in rest frame
    vrrest = std::sqrt(-1. * met.g11 / met.g00) * vr;
    //std::cout << "vr = " << vr << std::endl;
    //std::cout << "vrrest = " << vrrest << std::endl;
    if ((vrrest) > 1) {
      std::cerr << "Rest frame radial velocity greater than 1!" << std::endl;
    }

    // pixelization
    double dtheta = M_PI / double(ntheta);
    std::vector<double> theta(ntheta);
    unsigned long n = {0};
    std::generate(theta.begin(), theta.end(), [&n, dtheta]{ return double(n++) * dtheta + dtheta / 2.;});
    wdoublevec(theta, "theta.dat");

    double domega_norm = 2. * M_PI * dtheta;

    for (auto it = theta.cbegin(); it != theta.cend(); ++it) {
        sint = std::sin(*it);
        domegaarr.push_back(sint * domega_norm);

        if (std::abs(vr) < 1e-6) {
            q = calq(a, height, sint * sint);
            plus = ((*it <= PI / 2.) ? true : false);
            eo = std::sqrt(-1. / met.g00);
        } else {
            q = calq(a, height, vr, *it, umu, er, kmu, etheta);
            plus = (kmu[1] > 0.);
            if (!todisk) {
                kmu[1] *= -1.;
                kmu[2] *= -1.;
            }
            eo = -1. * dotprod(umu.data(), kmu.data(), &met);
        }

        inte = intet(q, a);
        quadroots roots(a, 0., q);
        roots.sort();

        re = solver(inte, a, height, 0., q, plus, roots, ncase);

        if (re >= disk.rms) {
            gfac = eo * disk.gfactor(re, 0.);
            garr.push_back(gfac);
            rearr.push_back(re);
        } else {
            garr.push_back(0.);
            rearr.push_back(0.);
        }
    }
}
*/

/*
//! Return metric, tetrad, and four-velocity of an stationary observer at a
//! off-axis lamp-post.
//!
//! @param a BH spin
//! @param height lamp-post height
//! @param s x-coordinate of the lamp-post
//! @param met output; Kerr metric (type: sim5metric)
//! @param tet output; ZAMO tetrad (type: sim5tetrad)
//! @param four_v output; four velocity of the observer at rest frame
void lpoffaxis(const double & a, const double & height, const double & s, sim5metric & met,
	sim5tetrad & tet, std::vector<double> & four_v) {
    // BL coordinate r and m=cos(theta)
    double r, mu;

    r = sqrt(s * s + height * height);
    mu = height / r;

    kerr_metric(a, r, mu, &met);
    tetrad_zamo(&met, &tet);

    // zamo four-velocity
    fourvelocity_zamo(&met, four_v.data());
}
*/

/*
//! Return metric, tetrad, and four-velocity of an observer at a
//! off-axis lamp-post with azimuthal motion about the z-axis, with angular
//! velocity \f$\Omega\f$.
//!
//! @param a BH spin
//! @param height lamp-post height
//! @param s x-coordinate of the lamp-post
//! @param omega angular velocity for azimuthal motion
//! @param met output; Kerr metric (type: sim5metric)
//! @param tet output; ZAMO tetrad (type: sim5tetrad)
//! @param four_v output; four velocity of the observer at rest frame
void lpoffaxis(const double & a, const double & height, const double & s, const double & omega, sim5metric & met,
	sim5tetrad & tet, std::vector<double> & four_v) {
    // take care to see if four_v has four elements;
    if (four_v.size() != 4)
        std::cerr << "four_v must have four elements!" << std::endl;
    // BL coordinate r and m=cos(theta)
    double r, mu;

    r = sqrt(s * s + height * height);
    mu = height / r;

    kerr_metric(a, r, mu, &met);
    tetrad_azimuthal(&met, omega, &tet);

    // zamo four-velocity
    fourvelocity_azimuthal(omega, &met, four_v.data());
}
*/

/*
//! Integrate all pixels across the 4pi sky to obtain the spectrum. This is the workhorse for
//! all kinds of off-axis corona. One just needs to feed the function with
//! disk info (type diskobj), local tetrad (type sim5tetrad), and metric (type sim5metric).
//!
//! As for the details, see document of lpaxis(). There is
//! yet one major difference: for off-axis lamp-post the turning point on theta is not the lamp-post
//! itself, such that one has to handle the theta-integral more carefully.
//! @param en energy array
//! @param flux output flux array
//! @param disk `diskobj` object. As it it a virtual base class, this is a general interface for
//!         any kinds of disk.
//! @param ntheta theta array size (the polar angle for photon wave vector as observed in rest frame
//! @param nphi phi array size (the azimuthal angle for photon wave vector as observed in rest frame
//! @param tet local tetrad
//! @param met metric
//! @param four_v four velocity of the observer in Boyer-Lindquist frame
//! @param caldphi
void offaxisgrid(const diskobj & disk, sim5tetrad & tet, long & ntheta, const long & nphi,
    std::vector<double> & rearr, std::vector<double> & garr, std::vector<double> & domegaarr) {

    rearr.clear();
    garr.clear();
    domegaarr.clear();

    double gfac, eo, diskg;

    // wave vector in rest frame and bl frame
    double on_k[4], bl_k[4];
    std::vector<double> four_v(4);
    double knorm;
    // constants of motions
    double l, q, re; //, dphimu, dphir, diskphi;
    // signs of ktheta, kr
    int signtheta, signr, ncase;
    // integrals;
    double intmu;

    // pixelization
    double dphi = 2. * PI / (double)(nphi);
    double sintnow;
    std::vector<double> phi(nphi);
    unsigned long n = {0};
    std::generate(phi.begin(), phi.end(), [&n, dphi]{ return double(n++) * dphi + dphi/ 2.;});

    double dtheta = M_PI / (double)(ntheta);
    std::vector<double> theta(ntheta);
    n = 0;
    std::generate(theta.begin(), theta.end(), [&n, dtheta]{ return double(n++) * dtheta + dtheta / 2.;});

	// calculate sin(theta)
    std::vector<double> sintheta(ntheta);
    std::transform(theta.cbegin(), theta.cend(), sintheta.begin(), [](double x){return std::sin(x);});

    for (long itheta = 0; itheta < ntheta; ++itheta) {
        sintnow = sintheta[itheta];
        for (long iphi = 0; iphi < nphi; ++iphi) {
            domegaarr.push_back(sintnow * dphi * dtheta);
            on_k[0] = 1.;
            on_k[1] = cos(theta[itheta]);
            on_k[2] = sin(theta[itheta]) * cos(phi[iphi]);
            on_k[3] = sin(theta[itheta]) * sin(phi[iphi]);

            on2bl(on_k, bl_k, &tet);
            // sign of r, theta
            knorm = tet.metric.g00 * bl_k[0] + tet.metric.g03 * bl_k[3];
            signr = (bl_k[1] > 0 ? 1 : -1);
            signtheta = (bl_k[2] > 0 ? 1 : -1);

            // photon energy measured by zamo observer
            eo = dotprod(bl_k, tet.e[0], &tet.metric) / knorm;

            // calculate constants of motion
            photon_motion_constants(tet.metric.a, tet.metric.r, tet.metric.m, bl_k, &l, &q);
            quadroots roots(tet.metric.a, l, q);
            roots.sort();

            if (q < 0.0)
                re = 0.;
            else {
                intmu = calintmu(tet.metric.a, l, q, tet.metric.m, signtheta);
                re = solver(intmu, tet.metric.a, tet.metric.r, l, q, (signr > 0 ? true : false),
                    roots, ncase);
            }

#ifdef DEBUG
            std::cout << "Integral over mu = " << intmu << std::endl;
            std::cout << "re = " << re << std::endl;
            std::cout << " ============================== " << std::endl << std::endl;
#endif
            if (re >= disk.rms) {
                diskg = disk.gfactor(re, l);
                gfac = eo * diskg;
                garr.push_back(gfac);
                rearr.push_back(re);
            } else {
                garr.push_back(0.);
                rearr.push_back(0.);
            }
        }
    }
    wdoublevec(theta, "theta.dat");
    wdoublevec(phi, "phi.dat");
}
*/

/*
//! Returns array of impact factors alpha, beta, mu-coordinate of the intersection between
//! the geodesic and the slab, and redshift. For I/O simplicity, all arrays are 1D.
void calmuarr_inftoslab(const double a, const double h, const double smax, const double muobs,
    const twodgeo & twodg, const long nalpha, const long nbeta, const double alphamax,
    double & domega, std::vector<double> & muarr, std::vector<double> & garr) {

    std::vector<double> alpha(nalpha), beta(nbeta);
    std::vector<double> fourv(4), kmu(4);
    geoinf geo;
    sim5tetrad tet;
    double r, munow, signr, signtheta, gslab;
    bool success;
    
    muarr.clear();
    garr.clear();

    long n = 0;
    double dalpha = 2. * alphamax / (double)(nalpha);
    std::generate(alpha.begin(), alpha.end(), [&n, dalpha, alphamax]{ return double(n++) * dalpha - alphamax + dalpha / 2.;});

    double dbeta = 2. * alphamax / (double)(nbeta);
    n = 0;
    std::generate(beta.begin(), beta.end(), [&n, dbeta, alphamax]{ return double(n++) * dbeta - alphamax + dbeta / 2.;});
    domega = dalpha * dbeta;

    for (auto ia = alpha.cbegin(); ia != alpha.cend(); ++ia) {
        for (auto ib = beta.cbegin(); ib != beta.cend(); ++ib) {
            std::cout << "alpha = " << *ia << ", beta = " << *ib << std::endl;
            geo = geoinf(a, muobs, *ia, *ib);
            success = false;
            
            try {
                r = calslabrad_sim5(geo, h, smax, signr, signtheta);
                std::cout << "r = " << r << std::endl;
                muarr.push_back(r);
                munow = h / r;
                if (munow > geo.muplus) 
                    munow = geo.muplus;
                tet = twodg.caltet(r, munow);
                photon_momentum(geo.a, r, munow, geo.l, geo.q, signr, signtheta, kmu.data());
                gslab = -1. * dotprod(tet.e[0], kmu.data(), &tet.metric);
                garr.push_back(1. / gslab);
                success = true;
            }
            catch (char const * error) {
                success = false;
                std::cerr << error << std::endl;
            }
            
            if (!success) {
                muarr.push_back(0.);
                garr.push_back(0.);
            }
        }
    }
}
*/

/*
void calmuarr_disktoslab(const double a, const double re, const double h, const double smax,
    const twodgeo & twodg, const long ntheta, const long nphi, std::vector<double> & muarr,
    std::vector<double> & garr, std::vector<double> & theta, std::vector<double> & domegaarr) {
    
    kepslab disktwodg(a);
    geodisk geo;
    sim5tetrad disktet, slabtet;
    
    std::vector<double> phi(nphi), kmudisk(4), kmuslab(4);
    double dtheta, dphi, mu, rsign, tsign, rnow, edisk, eslab;
    long n = 0;
    
    disktet = disktwodg.caltet(re, 0.);
    
    theta.clear();
    theta.resize(ntheta);
    dtheta = M_PI / 2. / (double)(ntheta);
    dphi = 2. * M_PI / (double)(nphi);
    std::generate(theta.begin(), theta.end(), [&n, dtheta]{ return (double)(n++) * dtheta + dtheta / 2.;});
    
    n = 0;
    std::generate(phi.begin(), phi.end(), [&n, dphi]{ return (double)(n++) * dphi + dphi / 2.;});
    
    muarr.clear();
    garr.clear();
    domegaarr.clear();
    
    for (auto it = theta.cbegin(); it != theta.cend(); ++it) {
        for (auto ip = phi.cbegin(); ip != phi.cend(); ++ip) {
            std::cout << "theta = " << *it << ", phi = " << *ip << std::endl;
            geo = geodisk(a, re, *it, *ip);
            std::cout << "geo.kmu[2] = " << geo.kmu[2] << std::endl;
            try {
                rnow = geo.calslabrad_sim5(h, smax, rsign, tsign);
                //std::cout << "mu = " << mu << std::endl;
                mu = h / rnow;
                if (mu > geo.muplus) 
                    mu = geo.muplus;
                kmudisk = geo.kmu;
                kmudisk[1] *= -1.;
                kmudisk[2] *= -1.;
                edisk = -1. * dotprod(kmudisk.data(), disktet.e[0], &disktet.metric);
                
                std::cout << "rnow = " << rnow << ", mu = " << mu << ", ktheta = " << kmudisk[2] << std::endl;
                photon_momentum(a, rnow, mu, geo.l, geo.q, rsign, tsign, kmuslab.data());
                slabtet = twodg.caltet(rnow, mu);
                eslab = -1. * dotprod(kmuslab.data(), slabtet.e[0], &slabtet.metric);
                
                muarr.push_back(mu);
                garr.push_back(edisk / eslab);
                domegaarr.push_back(std::sin(*it) * dtheta * dphi);
            }
            
            catch(char const* error) {
                muarr.push_back(0.);
                garr.push_back(0.);
                domegaarr.push_back(0.);
                std::cerr << error << std::endl;
            }
        }
    }
}
*/

/*!
\brief Pre-calculates several values related with photon geodesic, for photons emitted by disc on the equatorial plane.

The disc is divided into \f$N_r\f$ radial bins between \f$r_{\rm in}\f$ and \f$r_{\rm out}\f$. At the centroid of
each radial bin, we pixelate the sky in the disc fluid rest frame by dividing the local polar emission angle between
0 and \f$\theta_{\rm max}\f$ into \f$N_\theta\f$ grids and dividing the local azimuthal emission angle between 0 and \f$2\pi\f$
into \f$N_\phi\f$ grids. At last we calculate the geodesic of photon emitted in the center of each pixel. In total we
have \f$N_{\rm total} = N_r * N_\theta * N_\phi\f$ geodesics.
@param a BH spin
@param rin \f$r_{\rm in}\f$ in \f$[\rm GM~^{-2}]\f$
@param rout \f$r_{\rm out}\f$ in \f$[\rm GM~^{-2}]\f$
@param thetamax \f$\theta_{\rm max}\f$
@param nr \f$N_r\f$: number of radial bins
@param ntheta \f$N_\theta\f$: number of bins in local polar emission angle
@param nphi \f$N_\phi\f$: number of bins in local polar emission angle
@param muobsarr output; vector with \f$N_{\rm total}\f$ elements. If the photon escape to infinity, then the corresponding 
    element contains the inclination at infinity. It the photon arrives at the disc, the element is the radius on the disc. If the
    photon enters the BH event horizon, the element is -2. The index of muobsarr is related with indices of \f$r,\theta,\phi\f$ by
    \f$i= i_r * (N_\theta + i_\theta * N_\phi) + i_\phi\f$;
@param radius output; vector with \f$N_r\f$ elements. Centroid of radial bins.
@param theta output; vector with \f$N_\theta\f$ elements. Centroid of \f$\theta\f$ bins.
@param phi output; vector with \f$N_\phi\f$ elements. Centroid of \f$\phi\f$ bins.
@param rr1arr output; vector with \f$N_{\rm total}\f$ elements. This vector contains \f$r_1\f$, the largest root of \f$R(r)=0\f$.
    For definition of \f$R(r)\f$, see \ref kerr_rintegrand.
@param muplusarr output; vector with \f$N_\phi\f$ elements. This vector contains \f$\mu_+\f$, that largest turning point in \f$\mu\f$.
@param larr vector of \f$l\equiv L_z/E_\infty\f$ with \f$N_\phi\f$ elements.
@param qarr vector of \f$l\equiv \mathcal{Q}/E_\infty\f$ with \f$N_\phi\f$ elements.
@param kmuobssignarr vector with \f$N_\phi\f$ elements containing sign of \f$k^\theta\f$ at infinity or disc.
 */
void calmuobsarr(const double a, const double rin, const double rout, const double thetamax, const long nr, const long ntheta, 
    const long nphi, std::vector<double> & muobsarr, std::vector<double> & radius, std::vector<double> & theta,
    std::vector<double> & phi, std::vector<double> & rr1arr, std::vector<double> & muplusarr, 
    std::vector<double> & larr, std::vector<double> & qarr, std::vector<double> & kmuobssignarr) {
    
    double dlogr = std::log(rout / rin) / (double)(nr), logrmin = std::log(rin);
    double x, omega, sint, muinf;
    sim5metric met;
    sim5tetrad tet;
    std::array<double, 4> on_k, bl_k;
    on_k[0] = 1.;
    geoinf geo;
    long ngrid = nr * ntheta * nphi;
    std::vector<double> logradius(nr);
    
    theta.clear();
    theta.resize(ntheta);
    phi.clear();
    phi.resize(nphi);
    radius.clear();
    radius.resize(nr);
    muobsarr.resize(ngrid);
    rr1arr.resize(ngrid);
    muplusarr.resize(ngrid);
    larr.resize(ngrid);
    qarr.resize(ngrid);
    kmuobssignarr.resize(ngrid);
    
    long n = 0;
    std::generate(logradius.begin(), logradius.end(), [&n, dlogr, logrmin]{return double(n++) * dlogr + logrmin + dlogr / 2.;});
    std::transform(logradius.cbegin(), logradius.cend(), radius.begin(), [](double x){return std::exp(x);});
    
    double dphi = 2. * M_PI / (double)(nphi);
    double dtheta = thetamax * M_PI / (double)(ntheta);
    n = 0;
    std::generate(theta.begin(), theta.end(), [&n, dtheta]{return (double)(n++) * dtheta + dtheta / 2.;});
    n = 0;
    std::generate(phi.begin(), phi.end(), [&n, dphi]{return (double)(n++) * dphi + dphi / 2.;});
    
    n = 0;
    for (auto ir = radius.cbegin(); ir != radius.cend(); ++ir) {
        x = std::sqrt(*ir);
        omega = 1. / (x * x * x + a);
        kerr_metric(a, *ir, 0., &met);
        tetrad_azimuthal(&met, omega, &tet);
        
        for (auto it = theta.cbegin(); it != theta.cend(); ++it) {
            sint = std::sin(*it);
            on_k[2] = std::cos(*it);
            
            for (auto ip = phi.cbegin(); ip != phi.cend(); ++ip) {
                on_k[1] = sint * std::cos(*ip);
                on_k[3] = sint * std::sin(*ip);
                on2bl(on_k.data(), bl_k.data(), &tet);
                geo = geoinf(a, *ir, 0., bl_k);
                rr1arr[n] = geo.rr1;
                muplusarr[n] = geo.muplus;
                larr[n] = geo.l;
                qarr[n] = geo.q;
                
                muinf = geo.calfate(*ir, 0., bl_k[1], bl_k[2], rin);
                if (geo.status == 1) {
                    muobsarr[n] = muinf;
                    kmuobssignarr[n] = geo.kmuobssign;
                } else {
                    muobsarr[n] = -2.;
                    kmuobssignarr[n] = -1.;
                }
                ++n;
            }
        }
    }
    
    wdoublevec(muobsarr, "muobs.dat");
    wdoublevec(radius, "radius.dat");
    wdoublevec(theta, "theta.dat");
    wdoublevec(phi, "phi.dat");
}
                        
/*
void griddisk_3dgeo(const double a, const double m, const double mdot, const long nr, const long ntheta, const long nphi, const tridgeo & trid, 
    const detector & det) {
    std::vector<double> muobsarr, radius, theta, phi, rr1arr, muplusarr;
    std::vector<int> nmuturnarr;
    long indexnow;
    calmuobsarr(a, nr, ntheta, nphi, muobsarr, radius, theta, phi, rr1arr, muplusarr, nmuturnarr);
    for (long ir = 0; ir < radius.size(); ++ir) {
        for (long it = 0; it < theta.size(); ++it) {
            for (long ip = 0; ip < phi.size(); ++ip) {
                long indexnow = ir * (theta.size() + it * phi.size()) + ip;
                if ((kr >= 0. && radius[ir] > trid.rmax) (kr < 0. && rr1arr[indexnow] > trid.rmax)){
                    // photon package goes to infinity
                    det.receive(muobsarr[indexnow], pack);
*/

/*
void griddisk_nt(const double a, const double rout, const double m, const double mdot, 
    const long nr, const long ntheta, const long nphi, detector & det, const bool pol) {
    
    double cost, dphi, dtheta, dr, one_over_dr, domega, gfac;
    std::vector<double> muobsarr, radius, theta, phi, rr1arr, muplusarr, larr, qarr, poldegarr, muarr, darkarr, kmuobssignarr;
    long indexnow;
    
    nt disk(a, m, mdot);
    gammapack pack;
    
    std::cout << "Calculating geodesics..." << std::endl;
    calmuobsarr(a, disk.rms, rout, 0.5, nr, ntheta, nphi, muobsarr, radius, theta, phi, rr1arr, muplusarr, larr, qarr, kmuobssignarr);
    std::cout << "Done!" << std::endl;
    dphi = phi[1] - phi[0];
    dtheta = theta[1] - theta[0];
    dr = std::sqrt(radius[1] / radius[0]);
    one_over_dr = 1. / dr;
    
    if (pol) {
        muarr.resize(ntheta);
        darkarr.resize(ntheta);
        std::transform(theta.cbegin(), theta.cend(), muarr.begin(), [](double x){return std::cos(x);});
        poldegarr.resize(ntheta);
        calpoldeg(muarr, poldegarr);
        calchandra_darkening(muarr, darkarr);
        wdoublevec(darkarr, "dark.dat");
    }
        
    std::cout << "Calculating spectrum..." << std::endl;
    for (size_t ir = 0; ir < radius.size(); ++ir) {
        //std::cout << ir + 1 << "-th out of " << radius.size() << " radial bins\r";
        for (size_t it = 0; it < theta.size(); ++it) {
            cost = std::cos(theta[it]);
            domega = cost * dphi * (std::cos(theta[it] - dtheta / 2.) - std::cos(theta[it] + dtheta / 2.));
            
            for (size_t ip = 0; ip < phi.size(); ++ip) {
                indexnow = ir * theta.size() * phi.size() + it * phi.size() + ip;
                if (det.detected(muobsarr[indexnow])) {
                    gfac = disk.gfactor(radius[ir], larr[indexnow]);
                    if (pol) {
                        pack = gammapack(radius[ir], radius[ir] * one_over_dr, radius[ir] * dr, 
                            det.ne, det.emin / gfac, det.emax / gfac, disk, theta[it], phi[ip], larr[indexnow], poldegarr[it]);
                        pack.multiply(darkarr[it]);
                    } else {
                        pack = gammapack(radius[ir], radius[ir] * one_over_dr, radius[ir] * dr, 
                            det.ne, det.emin / gfac, det.emax / gfac, disk);
                    }
                        
                    pack.multiply(domega);
                    if (pol) {
                        //std::cout << "rnow = " << radius[ir] << ", tnow = " << theta[it] << ", pnow = " << phi[ip] << std::endl;
                        pack.calpolflux(muobsarr[indexnow], kmuobssignarr[indexnow], a, larr[indexnow], qarr[indexnow]);
                    }
                    det.receive_align(muobsarr[indexnow], pack);
                }
            }
        }
        showprogress((double)(ir + 1) / (double)(nr));
    }
    std::cout << std::endl;
    std::cout << "Done!" << std::endl;
}
*/

/*
//! Note that we include polarization anyway, as the scattering cross section depends on polarization degree, 
//! hence there is no switch for polarization calculation
//! @param a bh spin
//! @param m bh mass in solar mass
//! @param mdot mass accretion rate in Eddington rate
//! @param te temperature of thermal electron in [keV]
//! @param mfp scattering mean free path, in [rg]
//! @param trid `tridgeo` object; 
//! @param nr number of radial bin on the disc plane
//! @param ntheta number of theta bin
//! @param nphi number of phi bin
//! @param npack number of photon packages 
//! @param det `detector` object
void grid3dcorona(const double a, const double rout, const double m, const double mdot, const double te, 
    const double mfp, const tridgeo & trid, const long nr, const long ntheta, 
    const long nphi, const long npack, detector & det) {
    
    bool hit;
    int status;
    long n0scatter, indexnow;
    long ngrid = nr * ntheta * nphi;
    double knorm;
    double dlogr, dr, one_over_dr, rlo, rhi, dphi, dtheta, gfacdisk;
    double x, omega, sint, domega, gfac, domeganow;
    //double rhit, muhit;
    long nhit = 0;
    
    // grid vectors
    std::vector<double> muobsarr, rr1arr, muplusarr, larr, qarr, kmuobssignarr;
    // vectors record the information of scattering
    std::vector<double> muinfarr(npack), poldegnowarr(npack), lfinalarr(npack), qfinalarr(npack);
    std::vector<int> statusarr(npack), nscatterarr(npack), nmuturnarr1(npack);
    std::vector<std::complex<double>> wpnowarr(npack);
    std::vector<std::vector<double>> framegfacarr(npack), loboostarr(npack), inloboostarr(npack), scattermuarr(npack);
    
    std::complex<double> wpnow;
    std::array<double, 4> on_k, bl_k;
    std::vector<double> radius(nr), logradius(nr), theta(ntheta), phi(nphi), deten_emass(det.en),
        muarr(ntheta), darkarr(ntheta), poldegarr(ntheta);
        
    theta.resize(ntheta);
    phi.resize(nphi);
    radius.resize(nr);
    muobsarr.resize(ngrid);
    rr1arr.resize(ngrid);
    muplusarr.resize(ngrid);
    larr.resize(ngrid);
    qarr.resize(ngrid);
    kmuobssignarr.resize(ngrid);
    
    nt disk(a, m, mdot);
    
    sim5metric met;
    sim5tetrad tet;
    geoinf geo;
    gammapack pack;
    
    // random number generator
    std::random_device rd;
    std::mt19937_64 gen(rd());
    
    dlogr = std::log(rout / disk.rms) / (double)(nr);
    dr = std::exp(dlogr);
    one_over_dr = 1. / dr;
    
    dphi = 2. * M_PI / (double)(nphi);
    dtheta = M_PI / 2. / (double)(ntheta);
    
//    std::transform(det.en.cbegin(), det.en.cend(), deten_emass.begin(), 
//        std::bind2nd(std::multiplies<double>(), 1. / ME_KEV));
    
    std::cout << "Calculating geodesics..." << std::endl;
    calmuobsarr(a, disk.rms, rout, 0.5, nr, ntheta, nphi, muobsarr, radius, theta, phi, rr1arr, muplusarr, larr, qarr, kmuobssignarr);
    std::cout << "Done!" << std::endl;
    std::transform(theta.cbegin(), theta.cend(), muarr.begin(), [](double x){return std::cos(x);});
    calpoldeg(muarr, poldegarr);
    calchandra_darkening(muarr, darkarr);
    
    on_k[0] = 1.;
    
    std::cout << "Calculating spectrum..." << std::endl;
    for (size_t ir = 0; ir < radius.size(); ++ir) {
        rlo = radius[ir] * one_over_dr;
        rhi = radius[ir] * dr;
        x = std::sqrt(radius[ir]);
        omega = 1. / (x * x * x + a);
        kerr_metric(a, radius[ir], 0., &met);
        tetrad_azimuthal(&met, omega, &tet);
        
        for (size_t it = 0; it < theta.size(); ++it) {
            //std::cout << "theta = " << theta[it] << std::endl;
            sint = std::sin(theta[it]);
            on_k[2] = muarr[it];
            domega = muarr[it] * dphi * (std::cos(theta[it] - dtheta / 2.) - std::cos(theta[it] + dtheta / 2.));
            
            for (size_t ip = 0; ip < phi.size(); ++ip) {
                //std::cout << "ir = " << ir << ", it = " << it << ", ip = " << ip << std::endl;
                indexnow = ir * theta.size() * phi.size() + it * phi.size() + ip;
                n0scatter = 0;
                on_k[1] = sint * std::cos(phi[ip]);
                on_k[3] = sint * std::sin(phi[ip]);
                on2bl(on_k.data(), bl_k.data(), &tet);
                // normalize wave vector such that k_t = -1;
                knorm = -1. / (met.g00 * bl_k[0] + met.g03 * bl_k[3]);
                std::transform(bl_k.cbegin(), bl_k.cend(), bl_k.begin(), std::bind2nd(std::multiplies<double>(), knorm));
                
                geo = geoinf(a, radius[ir], 0., bl_k);
                
                hit = false;
                // If it is not possible for the photon to reach the corona
                if (geo.ispossible(trid, radius[ir], 0., bl_k[1], bl_k[2])) {
                    status = 0;
                    
                    //std::cout << "finding crossing point" << std::endl;
                    //if (ir == 6 && it == 48 && ip == 65) {
                    // status = geo.findxtrid(trid, radius[ir], 0., bl_k[1], bl_k[2], rhit, muhit, bl_k);
                    //std::cout << "finding crossing point ended" << std::endl;
                    //}
                    if (status != 1) {
                        hit = false;
                    } else {
                        // first we correct kmu; as we called findxtrid which utilizes step-by-step raytrace,
                        // kmu may no longer be a null vector (this is not guaranteed either) due to accumulation of errors.
                        // We `correct` this by leaving four-coordinate unchanged (note that numerical error accumulates in 
                        // both wave vector and coordinate) and call photon_momentum() to re-calculate the wave-vector.
                        // A more precise method would be to find the intersection by solving equations of motions.
                        // we start from the point where the photon first enters the corona
                        // the correction to photon wave vector is done within geoinf.findxtrid()
                        
                        // now we have intersection point: pos, and photon wavevector kmunow
                        // starting here we'll have <npack> photon packs
                        // for scattered photons we put the detector inside the routine as compton recoil is
                        // energy dependent actaully should see the difference afterwards by comparing with elastic thompson
                        // scattering 
                        
                        wpnow = disk.calwp_angle(radius[ir], theta[it], phi[ip], larr[indexnow]);
                        gfacdisk = 1. / disk.gfactor(radius[ir], larr[indexnow]);
                        std::cout << "re = " << radius[ir] << "theta = " << theta[it] << ", phi = " << phi[ip] << std::endl;
                        return;
                        for (long i = 0; i < npack; ++i) {
                            framegfacarr[i].push_back(gfacdisk);
                            try {
                                std::cout << "starting scattering process" << std::endl;
                                //tridsca(a, 1. / mfp, te, poldegarr[it], rhit, muhit, wpnow, bl_k, trid, gen, muinfarr[i], poldegnowarr[i], lfinalarr[i],
                                //    qfinalarr[i], statusarr[i], nscatterarr[i], nmuturnarr1[i], wpnowarr[i], framegfacarr[i], loboostarr[i], 
                                //    inloboostarr[i], scattermuarr[i]);
                                std::cout << "scattering process finished" << std::endl;
                            }
                            catch (const char * error) {
                                std::cout << "error in scattering process" << std::endl;
                                std::cout << error << std::endl;
                            }
                                
                        }
                    }
                }
                        
                domeganow = domega;
                if (hit) {
                    // calculate number of 0
                    for (long i = 0; i < npack; ++i) {
                        if (nscatterarr[i] == 0) {
                            ++n0scatter;
                        } else {
                            if (statusarr[i] == 1 && det.detected(muinfarr[i])) {
                            }
                        }
                    }
                    domeganow = domega * (double)(n0scatter) / (double)(npack);
                } 
                
                if ((!hit || n0scatter > 0) && muobsarr[indexnow] > 0. && det.detected(muobsarr[indexnow])) {
                    gfac = disk.gfactor(radius[ir], larr[indexnow]);
                    pack = gammapack(radius[ir], rlo, rhi, det.ne, det.emin / gfac, det.emax / gfac, disk, theta[it], 
                        phi[ip], larr[indexnow], poldegarr[it]);
                    pack.multiply(domeganow);
                    pack.calpolflux(muobsarr[indexnow], kmuobssignarr[indexnow], a, larr[indexnow], qarr[indexnow]);
                    det.receive_align(muobsarr[indexnow], pack);
                }
                showprogress((double)(indexnow) / (double)(ngrid));
            }
        }
    }
    std::cout << "Nhit = " << nhit << std::endl;
    std::cout << "Done!" << std::endl;
}
*/

/*
// Integration for self-irradiation
void selfirr(const diskobj & disk, const double fcol, 
    const double rtr, const long nr, const long ntheta, const long nphi, const long ne, 
    const double emin, const double emax) {
    
    double a = disk.a;
    double dphi, dtheta, domega, sint;
    long n;
    double re = 2.5, x = std::sqrt(re), rf, xf, omega, omega_rf, knorm;
    double de;
    double en_re, en_rf, teff, gfac, specnorm, krsign;
    std::array<double, 4> on_k, bl_k;
    std::array<double, 4> umu, umu_rf;
    
    geoinf geo;
    sim5metric met, met_rf;
    sim5tetrad tet;
    
    dphi = 2. * M_PI / (double)(nphi);
    dtheta = M_PI / 2. / (double)(ntheta);
    
    std::vector<double> theta(ntheta), muarr(ntheta), 
        phi(nphi), en(ne), flux(ne);
    std::vector<double> rfarr, rayrfarr;
    std::vector<double> rarrray, muarrray, rarrint, muarrint;
    
    n = 0;
    std::generate(theta.begin(), theta.end(), [&n, dtheta]{return (double)(n++) * dtheta + dtheta / 2.;});
    std::transform(theta.cbegin(), theta.cend(), muarr.begin(), [](double x){return std::cos(x);});
    n = 0;
    std::generate(phi.begin(), phi.end(), [&n, dphi]{return (double)(n++) * dphi + dphi / 2.;});
    
    de = (emax - emin) / (double)(ne);
    n = 0;
    std::generate(en.begin(), en.end(), [&n, de, emin]{ return (double)(n++) * de + de / 2. + emin;});
    
    omega = 1. / (x * x * x + a);
    kerr_metric(a, re, 0., &met);
    tetrad_azimuthal(&met, omega, &tet);
    
    double c2h3 = std::pow(10., logc2h3);
    // spectrum normalization
    specnorm = KEV2ERG * c2h3 / pow(fcol, 4);

    on_k[0] = 1.;
    
    for (long itheta = 0; itheta < ntheta; ++itheta) {
        on_k[2] = muarr[itheta];
        domega = muarr[itheta] * dphi * (std::cos(theta[itheta] - dtheta / 2.) - std::cos(theta[itheta] + dtheta / 2.));
        sint = std::sin(theta[itheta]);
        for (long iphi = 0; iphi < nphi; ++iphi) {
            on_k[1] = sint * std::cos(phi[iphi]);
            on_k[3] = sint * std::sin(phi[iphi]);
            on2bl(on_k.data(), bl_k.data(), &tet);
            // normalize wave vector such that k_t = -1;
            knorm = -1. / (met.g00 * bl_k[0] + met.g03 * bl_k[3]);
            std::transform(bl_k.cbegin(), bl_k.cend(), bl_k.begin(), std::bind2nd(std::multiplies<double>(), knorm));
            
            geo = geoinf(a, re, 0., bl_k);
            try {
                rf = geo.solver_selfirr(re, bl_k[1], krsign);
            }
            catch(const char * error) {
                rf = -1.;
            }
            
            if (rf > 0.) {
                rfarr.push_back(rf);
                xf = std::sqrt(rf);
                omega_rf = 1. / (xf * xf * xf + a);
                kerr_metric(a, rf, 0., &met_rf);
                // calculate redshift
                fourvelocity_azimuthal(omega, &met, umu.data());
                fourvelocity_azimuthal(omega_rf, &met_rf, umu_rf.data());
                
                en_re = (-1. * umu[0] + geo.l * umu[3]);
                en_rf = (-1. * umu_rf[0] + geo.l * umu_rf[3]);
                
                gfac = en_re / en_rf;
                teff = disk.teff_kev(rf) * fcol * gfac;
                
                for (long ie = 0; ie < ne; ++ie) {
                    flux[ie] += domega * en[ie] * en[ie] * en[ie] / (std::exp(en[ie] / teff) - 1.);
                }
            }
        }
    }
    std::transform(flux.begin(), flux.end(), flux.begin(),
        std::bind2nd(std::multiplies<double>(), specnorm));
    
    wdoublevec(flux, "flux.dat");
    wdoublevec(en, "en.dat");
}
*/
    
/*
void griddisk_nt_sp_pol(const double a, const double rout, const double m, const double mdot, const double fcol,
    const long nr, const long ntheta, const long nphi, const long nphoton, std::ofstream & enfile,
    std::ofstream & weightfile, std::ofstream & muinffile, std::ofstream & qweightfile, std::ofstream & uweightfile) {
    
    double cost, dphi, dtheta, dr, one_over_dr, domega, gfac, vphi, gamma;
    double area, teff, temnow, ennow, one_over_ut, Omega, kerrx;
    double qfac, ufac, qweight, uweight, muinf, weightnorm, weightnow;
    std::vector<double> muobsarr, radius, theta, phi, rr1arr, muplusarr, larr, qarr, poldegarr, muarr, darkarr, kmuobssignarr;
    long indexnow;
    std::mt19937_64 gen;
    std::complex<double> wp;
    
    nt disk(a, m, mdot);
    
    std::cout << "Calculating geodesics..." << std::endl;
    calmuobsarr(a, disk.rms, rout, 0.5, nr, ntheta, nphi, muobsarr, radius, theta, phi, rr1arr, muplusarr, larr, qarr, kmuobssignarr);
    std::cout << "Done!" << std::endl;
    dphi = phi[1] - phi[0];
    dtheta = theta[1] - theta[0];
    dr = std::sqrt(radius[1] / radius[0]);
    one_over_dr = 1. / dr;
    
    muarr.resize(ntheta);
    darkarr.resize(ntheta);
    std::transform(theta.cbegin(), theta.cend(), muarr.begin(), [](double x){return std::cos(x);});
    poldegarr.resize(ntheta);
    calpoldeg(muarr, poldegarr);
    calchandra_darkening(muarr, darkarr);
    
    std::cout << "Calculating spectrum..." << std::endl;
    for (size_t ir = 0; ir < radius.size(); ++ir) {
        teff = disk.teff_kev(radius[ir]);
        temnow = teff * fcol / ME_KEV;
        area = calds(a, radius[ir] * one_over_dr, radius[ir] * dr);
        vphi = disk.vphi(radius[ir]);
        gamma = 1. / std::sqrt(1. - vphi * vphi);
        kerrx = std::sqrt(radius[ir]);
        Omega = 1. / (kerrx * kerrx * kerrx + a);
        one_over_ut = std::sqrt(1. - 2. / radius[ir] + 4. * Omega * a / radius[ir] - Omega * Omega * (radius[ir] * radius[ir] + 
            a * a + 2. * a * a / radius[ir]));
        weightnorm = m * m * area * gamma * one_over_ut * teff * teff * teff / fcol / (double)(nphoton);

        //std::cout << ir + 1 << "-th out of " << radius.size() << " radial bins\r";
        for (size_t it = 0; it < theta.size(); ++it) {
            cost = std::cos(theta[it]);
            domega = cost * dphi * (std::cos(theta[it] - dtheta / 2.) - std::cos(theta[it] + dtheta / 2.));
            weightnow = weightnorm * darkarr[it] * domega;
            
            for (size_t ip = 0; ip < phi.size(); ++ip) {
                indexnow = ir * theta.size() * phi.size() + it * phi.size() + ip;
                gfac = disk.gfactor(radius[ir], larr[indexnow]);

                if (muobsarr[indexnow] > 0.) {
                    // construct local superphotons
                    // first randomize energy
                    muinf = muobsarr[indexnow];
                    wp = disk.calwp_angle(radius[ir], theta[it], phi[ip], larr[indexnow]);
                    cal_pol_inf(muobsarr[indexnow], kmuobssignarr[indexnow], a, larr[indexnow], qarr[indexnow], 
                        poldegarr[it], wp, qfac, ufac);
                    qweight = qfac * weightnow;
                    uweight = ufac * weightnow;
                    for (long iphoton = 0; iphoton < nphoton; ++iphoton) {
                        ennow = sample_bb(temnow, gen) * gfac;
                        enfile.write(reinterpret_cast<char*>(&ennow), sizeof(ennow));
                        weightfile.write(reinterpret_cast<char*>(&weightnow), sizeof(weightnow));
                        muinffile.write(reinterpret_cast<char*>(&muinf), sizeof(muinf));
                        
                        qweightfile.write(reinterpret_cast<char*>(&qweight), sizeof(qweight));
                        uweightfile.write(reinterpret_cast<char*>(&uweight), sizeof(uweight));
                    }
                }
            }
        }
        showprogress((double)(ir + 1) / (double)(nr));
    }
    std::cout << std::endl;
    std::cout << "Done!" << std::endl;
}
*/


/*
//! Note that we include polarization anyway, as the scattering cross section depends on polarization degree, 
//! hence there is no switch for polarization calculation
//! This routine uses the ``superphoton'' approach. And to increase SNR at high energy, we 
//! include bias for Compton scattering.
//! @param a bh spin
//! @param m bh mass in solar mass
//! @param mdot mass accretion rate in Eddington rate
//! @param te temperature of thermal electron in [keV]
//! @param mfp scattering mean free path, in [rg]
//! @param trid `tridgeo` object; 
//! @param nr number of radial bin on the disc plane
//! @param ntheta number of theta bin
//! @param nphi number of phi bin
//! @param npack number of photon packages 
//! @param det `detector` object
void grid3dcorona_sp(const bool pol, const double a, const double rin, const double rout, 
    const double m, const double mdot, const double fcol, const double te, 
    const double mfp, const tridgeo & trid, const long nr, const long ntheta, const long nphi, 
    const long nphoton, std::vector<std::ofstream> & ofiles, std::ofstream & nscafile) {
    
    bool hit;
    int status;
    long indexnow, ngrid = nr * ntheta * nphi;
    double knorm;
    double dlogr, dr, one_over_dr, rnow, dphi, dtheta, gdisk, ghit, vphi, omega, area, muobsnow, weightnorm, weightnorm0;
    double x, teff, temnow, temnow_emass, sint, domega, gfac, gamma, one_over_ut, ennow, enhit, qfac, ufac, 
        qweight, uweight;
    double rhit, muhit;
    int num0 = 0;
    
    // grid vectors
    std::vector<double> muobsarr, rr1arr, muplusarr, larr, qarr, kmuobssignarr;
    
    std::complex<double> wpnow;
    std::array<double, 4> on_k, bl_k, umu, kmu_hit;
    std::vector<double> radius(nr), logradius(nr), theta(ntheta), phi(nphi), muarr(ntheta);
    
    std::vector<double> darkarr(ntheta), poldegarr(ntheta);
        
    theta.resize(ntheta);
    phi.resize(nphi);
    radius.resize(nr);
    muobsarr.resize(ngrid);
    rr1arr.resize(ngrid);
    muplusarr.resize(ngrid);
    larr.resize(ngrid);
    qarr.resize(ngrid);
    kmuobssignarr.resize(ngrid);
    
    nt disk(a, m, mdot);
    std::cout << "SPECNORM_SP = " << SPECNORM_SP << std::endl;
    
    sim5metric met, metnow;
    sim5tetrad tet;
    geoinf geo;
    
    // random number generator
    std::random_device rd;
    std::mt19937_64 gen(rd());
    
    dlogr = std::log(rout / disk.rms) / (double)(nr);
    dr = std::exp(dlogr * 0.5);
    one_over_dr = 1. / dr;
    
    dphi = 2. * M_PI / (double)(nphi);
    dtheta = M_PI / 2. / (double)(ntheta);
    
    std::cout << "Calculating geodesics..." << std::endl;
    calmuobsarr(a, disk.rms, rout, 0.5, nr, ntheta, nphi, muobsarr, radius, theta, phi, rr1arr, muplusarr, larr, qarr, kmuobssignarr);
    std::cout << "Done!" << std::endl;
    std::transform(theta.cbegin(), theta.cend(), muarr.begin(), [](double x){return std::cos(x);});
    
    if (pol) {
        calpoldeg(muarr, poldegarr);
        calchandra_darkening(muarr, darkarr);
    }
    
    on_k[0] = 1.;
    
#pragma omp parallel
{
#pragma omp single
{
    std::cout << "Scattering..." << std::endl;
    for (size_t ir = 0; ir < radius.size(); ++ir) {
        std::cout << "ir = " << ir << std::endl;
        rnow = radius[ir];
        teff = disk.teff_kev(rnow);
        temnow = teff * fcol;
        temnow_emass = temnow / ME_KEV;
        area = calds(a, rnow * one_over_dr, rnow * dr);
        vphi = disk.vphi(rnow);
        gamma = 1. / std::sqrt(1. - vphi * vphi);
        
        x = std::sqrt(rnow);
        omega = 1. / (x * x * x + a);
        
        one_over_ut = std::sqrt(1. - 2. / rnow + 4. * omega * a / rnow - omega * omega * (rnow * rnow + 
            a * a + 2. * a * a / rnow));
        
        kerr_metric(a, rnow, 0., &met);
        tetrad_azimuthal(&met, omega, &tet);
        
        weightnorm0 = m * m * area * gamma * one_over_ut * teff * teff * teff / fcol / (double)(nphoton);

        for (size_t it = 0; it < theta.size(); ++it) {
            //std::cout << "theta = " << theta[it] << std::endl;
            sint = std::sin(theta[it]);
            on_k[2] = muarr[it];
            domega = muarr[it] * dphi * (std::cos(theta[it] - dtheta / 2.) - std::cos(theta[it] + dtheta / 2.));
            weightnorm = weightnorm0 * domega * darkarr[it];
            
            for (size_t ip = 0; ip < phi.size(); ++ip) {
                //std::cout << "ir = " << ir << ", it = " << it << ", ip = " << ip << std::endl;
                indexnow = ir * theta.size() * phi.size() + it * phi.size() + ip;
                on_k[1] = sint * std::cos(phi[ip]);
                on_k[3] = sint * std::sin(phi[ip]);
                on2bl(on_k.data(), bl_k.data(), &tet);
                // normalize wave vector such that k_t = -1;
                knorm = -1. / (met.g00 * bl_k[0] + met.g03 * bl_k[3]);
                std::transform(bl_k.cbegin(), bl_k.cend(), bl_k.begin(), std::bind2nd(std::multiplies<double>(), knorm));
                gdisk = disk.gfactor(radius[ir], larr[indexnow]);
                
                geo = geoinf(a, radius[ir], 0., bl_k);
                if (geo.l != geo.l || geo.q != geo.q) {
                    std::cout << "Inside grid3dcorona_sp; Checking geo.l:" << std::endl;
                    std::cout << "r = " << radius[ir] << ", mue = 0" << std::endl;
                    std::cout << "bl_k = " << bl_k << std::endl;
                    break;
                }
                if (pol) {
                    if (muobsarr[indexnow] > 0.) {
                        wpnow = disk.calwp_angle(radius[ir], theta[it], phi[ip], larr[indexnow]);
                        cal_pol_inf(muobsarr[indexnow], kmuobssignarr[indexnow], a, larr[indexnow], qarr[indexnow], poldegarr[it],
                            wpnow, qfac, ufac);
                        
                        if (qfac != qfac || ufac != ufac) {
                            std::cout << "Inside gridgeod()" << std::endl;
                            std::cout << "wpnow = " << wpnow << std::endl;
                            std::cout << "muobsarr[indexnow] = " << muobsarr[indexnow] << std::endl;
                            std::cout << "larr[indexnow] = " << larr[indexnow] << std::endl;
                            std::cout << "qarr[indexnow] = " << qarr[indexnow] << std::endl;
                            std::cout << "poldegarr[it] = " << poldegarr[it] << std::endl;
                            std::cout << "zero scattering superphoton:" << std::endl;
                            std::cout << "qfac or ufac is NaN" << std::endl;
                        }
                    } else {
                        qfac = 0.;
                        ufac = 0.;
                    }
                }
                
                hit = false;
                // If it is not possible for the photon to reach the corona
                if (!geo.ispossible(trid, radius[ir], 0., bl_k[1], bl_k[2])) {
                    hit = false;
                } else {
                    //std::cout << "finding crossing point" << std::endl;
                    status = geo.findxtrid(trid, radius[ir], 0., bl_k[1], bl_k[2], disk.rms, rhit, muhit, kmu_hit);
                    
                    if (status != 1) {
                        hit = false;
                    } else {
                        // first we correct kmu; as we called findxtrid which utilizes step-by-step raytrace,
                        // kmu may no longer be a null vector (this is not guaranteed either) due to accumulation of errors.
                        // We `correct` this by leaving four-coordinate unchanged (note that numerical error accumulates in 
                        // both wave vector and coordinate) and call photon_momentum() to re-calculate the wave-vector.
                        // A more precise method would be to find the intersection by solving equations of motions.
                        // we start from the point where the photon first enters the corona
                        // the correction to photon wave vector is done within geoinf.findxtrid()
                        
                        // now we have intersection point: pos, and photon wavevector kmunow
                        // starting here we'll have <npack> photon packs
                        // for scattered photons we put the detector inside the routine as compton recoil is
                        // energy dependent actaully should see the difference afterwards by comparing with elastic thompson
                        // scattering 
                        
                        hit = true;
                        kerr_metric(a, rhit, muhit, &metnow);
                        trid.calumu(rhit, muhit, umu, metnow);
                        ghit = dotprod(umu.data(), kmu_hit.data(), &metnow);
                        gfac = -1. * ghit * gdisk;
                        //std::cout << "gfac = " << gfac << std::endl;
                        for (long i = 0; i < nphoton; ++i) {
                            ennow = sample_bb(temnow_emass, gen);
                            enhit = ennow * gfac;
                            superphoton sp(enhit);
                            sp.weight = 1.;
                            sp.rnow = rhit;
                            sp.munow = muhit;
                            sp.kmu = kmu_hit;
                            if (pol) {
                                sp.poldeg = poldegarr[it];
                                sp.wp = wpnow;
                            }
#pragma omp task
                            tridsca_sp(pol, a, rin, mfp, te, weightnorm, qfac, ufac, muobsarr[indexnow], dr,
                                sp, trid, ofiles, nscafile);
                        }
                    }
                }
                
                if ((!hit) && muobsarr[indexnow] > 0) {
                    muobsnow = muobsarr[indexnow];
                    for (long i = 0; i < nphoton; ++i) {
                        ennow = sample_bb(temnow_emass, gen) * gdisk;
#pragma omp critical
{
                        ofiles[0].write(reinterpret_cast<char*>(&ennow), sizeof(ennow));
                        ofiles[1].write(reinterpret_cast<char*>(&weightnorm), sizeof(weightnorm));
                        ofiles[2].write(reinterpret_cast<char*>(&muobsnow), sizeof(muobsnow));
                        nscafile.write(reinterpret_cast<char*>(&num0), sizeof(num0));
                        if (pol) {
                            qweight = weightnorm * qfac;
                            uweight = weightnorm * ufac;
                            ofiles[3].write(reinterpret_cast<char*>(&qweight), sizeof(qweight));
                            ofiles[4].write(reinterpret_cast<char*>(&uweight), sizeof(uweight));
                        }
                    }
}
                }
                //showprogress((double)(indexnow) / (double)(ngrid));
            }
        }
    }
}
}
    std::cout << "Done!" << std::endl;
}
*/

/*!
\brief Write to a file the parameters required to reconstruct energy and polarisation spectrum of thin disc on the equatorial plane.

Although this function has `nt` in its name, it is intended for any kind of emissivity profile.
We write the parameters into a binary parameter file. All data are in double-precision.
The file contains a *header* with 6 numbers, followed by N segments where N is the number of pixels that arrives at infinity.
The structure of the file is thus as follows:
\f{equation}
{a, r_{\rm out}, 0.,N_r, N_\theta, N_\phi, [r, g, \mu_\infty, w, {\rm cos}\Delta\psi, {\rm sin}\Delta\psi]\times N},
\f}
where in the *header*: \f$r_{\rm out}\f$ is the disc outer edge radius, \f$N_r, N_\theta, N_\phi\f$ are number of bins (see \ref calmuobsarr).
And in each segment: \f$r\f$ is the emission radius, \f$g\equiv E_\infty/E_{\rm disc}\f$ is the redshift factor, \f$\mu_\infty\f$ is the photon
inclination at infinity, \f$w\f$ is the photon statistical weight, and \f$\Delta\psi\f$ is the rotation angle of polarisation vector from disc to
infinity. Here \f$w = {\rm cos}\theta_{\rm em} d\mathcal{S} d\Omega/
u^t\f$, where , \f$\theta_{\rm em}\f$ is the local emission angle,
\f$\mathcal{S}\f$ is the proper area, and \f$d\Omega\f$ is the solid angle,
and \f$u^t\f$ is the factor to convert local time to distant observer's time. The photon generation rate per unit distant observer's time
\f$\dot{N} = \frac{4\zeta(3) k_{\rm B}^3}{c^2h^3}\frac{f_{\rm limb} T_{\rm eff}^3}{f_{\rm col}u^t} w\f$,
where \f$\zeta()\f$ is the Riemann zeta function, \f$f_{\rm limb}\f$ is the limb-darkening factor, \f$T_{\rm eff}\f$ is the effective temperature,
and \f$f_{\rm col}\f$ is the color correction factor.

@param pol the option to turn on/off polarisation calculation. Note that the directionarity is also affected by this option. If off,
    we assume isotropic emission from the disc in the disc fluid rest frame. Otherwise the angular dependence follows Chandrasekhar 1960
    (http://adsabs.harvard.edu/abs/1960ratr.book.....C) assuming a pure scattering atmosphere with infinite optical depth, which can be 
    calculated using \ref calchandra_darkening.
@param a BH spin
@param rout disk outer disc radius, in \f$[\rm GM~c^{-2}]\f$
@param thetamax the maximum polar emission angle \f$\theta_{\rm max}\f$; see \ref calmuobsarr.
@param nr number of radial bins \f$N_r\f$
@param ntheta number of \f$\theta\f$ bins
@param nphi number of \f$\phi\f$ bins
@param paramfile the name of the output parameter file.
 */
void griddisk_nt_sp_pol_writeparam(const bool pol, const double a, const double rout, 
    const double thetamax, const long nr, const long ntheta, const long nphi, std::string paramfile) {
    
    double cost, dphi, dtheta, dr, one_over_dr, domega, vphi, gamma;
    double area, one_over_ut, Omega, kerrx, muinf, weightnorm, weightnow, gfac;
    double c1, c2, cospsiinf, sinpsiinf, cospsi, sinpsi;
    std::vector<double> muobsarr, radius, theta, phi, rr1arr, muplusarr, larr, qarr, kmuobssignarr;
    long indexnow;
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::complex<double> wp;
    auto dsize = sizeof(double);
    nt disk(a, 1., 1.);
    
    std::cout << "Calculating geodesics..." << std::endl;
    calmuobsarr(a, disk.rms, rout, thetamax, nr, ntheta, nphi, muobsarr, radius, theta, phi, rr1arr, muplusarr, larr, qarr, kmuobssignarr);
    std::cout << "Done!" << std::endl;
    dphi = phi[1] - phi[0];
    dtheta = theta[1] - theta[0];
    dr = std::sqrt(radius[1] / radius[0]);
    one_over_dr = 1. / dr;
    
    std::ofstream ofile(paramfile, std::ios::out | std::ofstream::binary);
    
    std::vector<double> input_params = {a, rout, 0.,
        (double)(nr), (double)(ntheta), (double)(nphi)};
    ofile.write((char*)&input_params[0], input_params.size() * dsize);
    
    std::cout << "Calculating spectrum..." << std::endl;
    for (size_t ir = 0; ir < radius.size(); ++ir) {
        area = calds(a, radius[ir] * one_over_dr, radius[ir] * dr);
        vphi = disk.vphi(radius[ir]);
        gamma = 1. / std::sqrt(1. - vphi * vphi);
        kerrx = std::sqrt(radius[ir]);
        Omega = 1. / (kerrx * kerrx * kerrx + a);
        one_over_ut = std::sqrt(1. - 2. / radius[ir] + 4. * Omega * a / radius[ir] - Omega * Omega * (radius[ir] * radius[ir] + 
            a * a + 2. * a * a / radius[ir]));
        weightnorm = area * gamma * one_over_ut;

        //std::cout << ir + 1 << "-th out of " << radius.size() << " radial bins\r";
        for (size_t it = 0; it < theta.size(); ++it) {
            cost = std::cos(theta[it]);
            domega = std::abs(cost) * dphi * (std::cos(theta[it] - dtheta / 2.) - std::cos(theta[it] + dtheta / 2.));
            weightnow = weightnorm * domega;
            
            for (size_t ip = 0; ip < phi.size(); ++ip) {
                indexnow = ir * theta.size() * phi.size() + it * phi.size() + ip;

                if (muobsarr[indexnow] > -1.5) {
                    // construct local superphotons
                    // first randomize energy
                    muinf = muobsarr[indexnow];
                    gfac = disk.gfactor(radius[ir], larr[indexnow]);
                    
                    if (pol) {
                        disk.polrotang_disktok(radius[ir], theta[it], phi[ip], larr[indexnow], qarr[indexnow], c1, c2);
                        disk.polrotang_ktoinf(kmuobssignarr[indexnow], larr[indexnow], qarr[indexnow], muinf, cospsiinf, sinpsiinf);
                        cospsi = c1 * cospsiinf - c2 * sinpsiinf;
                        sinpsi = c1 * sinpsiinf + c2 * cospsiinf;
                        //std::cout << cospsi * cospsi + sinpsi * sinpsi << std::endl;
                    } else {
                        cospsi = 0.;
                        sinpsi = 0.;
                    }
                    ofile.write(reinterpret_cast<char*>(&radius[ir]), dsize);
                    ofile.write(reinterpret_cast<char*>(&theta[it]), dsize);
                    ofile.write(reinterpret_cast<char*>(&gfac), dsize);
                    ofile.write(reinterpret_cast<char*>(&muinf), dsize);
                    ofile.write(reinterpret_cast<char*>(&weightnow), dsize);
                    ofile.write(reinterpret_cast<char*>(&cospsi), dsize);
                    ofile.write(reinterpret_cast<char*>(&sinpsi), dsize);
                }
            }
        }
        showprogress((double)(ir + 1) / (double)(nr));
    }
    std::cout << std::endl;
    ofile.close();
}

/*!
\brief Register blackbody photons, given blackbody temperatures, weight, inclination, and normalised Stokes parameters.

The weight written into the file is \f$w^\prime = \frac{1}{N_{\rm photon}} \frac{f_{\rm limb} {\rm cos}\theta_{\rm em} d\mathcal{S} d\Omega
 T_{\rm eff}^3}{f_{\rm col}u^t}\f$. To obtain the photon generation rate: \f$\dot{N}=\frac{4\zeta(3) k_{\rm B}^3}{c^2h^3} w^\prime\f$.

@param pol switch for polarization calculations
@param nphoton nubmer of photons 
@param tem color temperature as observed by the observer
@param weightnorm normalization of weight
@param qfac normalised Stokes parameter Q
@param ufac normalised Stokes parameter U
@param muinf \f$\mu_{\rm obs}\f$
@param gen random number generator
@param ofiles vector of output files; see \ref register_sp_nt.
*/
void register_sp_bbody(const bool pol, const long nphoton, const double tem, const double weightnorm, const double qfac, const double ufac, const
    double muinf, std::mt19937_64 & gen, std::vector<std::ofstream> & ofiles) {
    
    size_t double_size = sizeof(double);
    double eninf, muinf1 = muinf;
    double weight = weightnorm / (double)(nphoton);
    double qweight = qfac * weight, uweight = ufac * weight;
    
    if (muinf > -1.5) {
        for (long i = 0; i < nphoton; ++i) {
            eninf = sample_bb(tem, gen);
            ofiles[0].write(reinterpret_cast<char*>(&eninf), double_size);
            ofiles[1].write(reinterpret_cast<char*>(&weight), double_size);
            ofiles[2].write(reinterpret_cast<char*>(&muinf1), double_size);
            if (pol) {
                ofiles[3].write(reinterpret_cast<char*>(&qweight), double_size);
                ofiles[4].write(reinterpret_cast<char*>(&uweight), double_size);
            }
        }
    }
}

/*!
\brief Overloaded \ref register_sp_bbody(const bool,const long,const double,const double, const double,const double, const double, 
    std::mt19937_64 &, std::vector<std::ofstream> &).
    
In this one the blackbody photon is weighted sampled to ensure relatively more photons at high energy tail.
@param pol switch for polarization calculations
@param nphoton nubmer of photons 
@param tem color temperature as observed by the observer
@param weightnorm normalization of weight
@param qfac weight * qfac = Stokes parameter Q
@param ufac weight * ufac = Stokes parameter U
@param muinf \f$\mu_{\rm obs}\f$
@param ecut energy above which more photons will be sampled
@param ratio over-sampling ratio
@param gen random number generator
@param ofiles vector of output files; see \ref register_sp_nt.
*/
void register_sp_bbody(const bool pol, const long nphoton, const double tem, const double weightnorm, const double qfac, const double ufac, const
    double muinf, const double ecut, const double ratio, std::mt19937_64 & gen, std::vector<std::ofstream> & ofiles) {
    
    size_t double_size = sizeof(double);
    double eninf, sample_weight, muinf1 = muinf, weightnow, qweightnow, uweightnow;
    double weight = weightnorm / (double)(nphoton);
    double qweight = qfac * weight, uweight = ufac * weight;
    
    if (muinf > 0.) {
        for (long i = 0; i < nphoton; ++i) {
            sample_bb1(tem, gen, ecut, ratio, eninf, sample_weight);
            weightnow = weight * sample_weight;
            ofiles[0].write(reinterpret_cast<char*>(&eninf), double_size);
            ofiles[1].write(reinterpret_cast<char*>(&weightnow), double_size);
            ofiles[2].write(reinterpret_cast<char*>(&muinf1), double_size);
            if (pol) {
                qweightnow = sample_weight * qweight;
                uweightnow = sample_weight * uweight;
                ofiles[3].write(reinterpret_cast<char*>(&qweightnow), double_size);
                ofiles[4].write(reinterpret_cast<char*>(&uweightnow), double_size);
            }
        }
    }
}

/*!
\brief With the geodesic information saved in the parameter file produced by \ref griddisk_nt_sp_pol_writeparam, this function
calculates the values for Novikov-Thorne disc profile.

For each pixel that arrives at infinity, this function first calcualtes the effective temperature and the statistical weight 
given \f$a, M, \dot{M}, f_{\rm col}\f$. We assume Novikov-Thorne temperature profile (Page & Thorne 1974). 
If the polarisation option is switched off, the normalised Stokes parameters are also calculated and the limb-darkening factor is for
semi-infinite pure scattering disc (Chandrasekhar 1960; see also \ref calchandra_darkening). If the option is on we assume isotropic emission.
We multiply the weight \f$w\f$ calculated by \ref griddisk_nt_sp_pol_writeparam with a factor to obtain \f$w^\prime\f$:
\f$ w^\prime = \frac{1}{N_{\rm photon}} \frac{f_{\rm limb} T_{\rm eff}^3}{f_{\rm col}} w \f$.
The temperature, weight \f$w^\prime\f$, and Stokes parameters are passed to \ref register_sp_bbody.

@param pol whether to turn on polarisation
@param wsample whether to weighted sample blackbody spectrum
@param parafile the input parameter file
@param m BH mass in solar mass
@param mdot mass accretion rate in Eddington rate
@param fcol color correction factor
@param nphoton number of photons
@param gen random number generator
@param ofiles vector of std::ofstream, pointing to the output files:
    - ofiles[0]: dimensionless energies of the superphotons
    - ofiles[1]: statistical weights of the superphotons
    - ofiles[2]: photon inclination angle at infinity
    - If `pol` is turned on, `ofiles` has two additional elements:
    - ofiles[3]: Stokes parameter Q
    - ofiles[4]: Stokes parameter U
    
 */
void register_sp_nt(const bool pol, const bool wsample, const std::string parafile, const double m, const double mdot, const double fcol,
    const long nphoton, std::mt19937_64 & gen, std::vector<std::ofstream> & ofiles) {
    
    std::vector<double> params, theta, muarr, darkarr, poldegarr;
    double a, teff, tcol_emass, qfac, ufac, flimb, dtheta, weightnorm, ftheta, fphi, cos2psi, sin2psi, dit, poldegnow, theta1;
    rdoublevec(params, parafile);
    long ntheta, n = 0, nparam = 7, nseg, offset = 6, pos, it0, it1;
    nseg = (params.size() - offset) / nparam;
    
    // ntheta = (long)(params[4]);
    ntheta = 501;
    dtheta = 0.5 * M_PI / (double)(ntheta - 1);
    
    a = params[0];
    nt disk(a, m, mdot);
    std::cout << "a = " << a << std::endl;
    std::cout << "m = " << m << std::endl;
    
    if (pol) {
        theta.resize(ntheta);
        std::generate(theta.begin(), theta.end(), [&n, dtheta]{return (double)(n++) * dtheta;});
        muarr.resize(ntheta);
        darkarr.resize(ntheta);
        std::transform(theta.cbegin(), theta.cend(), muarr.begin(), [](double x){return std::abs(std::cos(x));});
        poldegarr.resize(ntheta);
        calpoldeg(muarr, poldegarr);
        calchandra_darkening(muarr, darkarr);
        //wdoublevec(theta, "poltheta.dat");
        //wdoublevec(muarr, "polmu.dat");
        //wdoublevec(poldegarr, "polpol.dat");
        //std::cout << "theta = " << theta << std::endl;
        //std::cout << "poldeg = " << poldegarr << std::endl;
    } else {
        qfac = 0.;
        ufac = 0.;
    }
    std::cout << "nseg = " << nseg << std::endl;
    
    for (long i = 0; i < nseg; ++i) {
        pos = i * nparam + offset;
        teff = disk.teff_kev(params[pos]);
        tcol_emass = fcol * teff / ME_KEV * params[pos + 2];;
        if (pol) {
            theta1 = (params[pos + 1] <= 0.5 * M_PI) ? params[pos + 1] : M_PI - params[pos + 1];
            dit = (theta1 - 0.5 * dtheta) / dtheta;
            it0 = (long)(std::floor(dit));
            it1 = it0 + 1;
            flimb = (darkarr[it0] * (theta1 - theta[it0]) + darkarr[it1] * (theta[it1] - theta1)) / dtheta;
            poldegnow = (poldegarr[it0] * (theta1 - theta[it0]) + poldegarr[it1] * (theta[it1] - theta1)) / dtheta;
            //std::cout << "thetanow = " << theta1 << ", poldegnow = " << poldegnow << std::endl;
            
            ftheta = -1. * params[pos + 6];
            fphi = params[pos + 5];
            sin2psi = 2. * ftheta * fphi;
            cos2psi = 2. * ftheta * ftheta - 1.;
            qfac = poldegnow * cos2psi;
            ufac = poldegnow * sin2psi;
            
        } else {
            flimb = 1.;
            qfac = 0.;
            ufac = 0.;
        }
        //std::cout << "muinf = " << params[pos + 3] << std::endl;
        weightnorm = params[pos + 4] * disk.m * disk.m  * flimb * teff * teff * teff / fcol;
        if (wsample) {
            register_sp_bbody(pol, nphoton, tcol_emass, weightnorm, qfac, ufac, params[pos + 3], 10., 10., gen, ofiles);
        } else {
            register_sp_bbody(pol, nphoton, tcol_emass, weightnorm, qfac, ufac, params[pos + 3], gen, ofiles);
        }
        //showprogress((double)(i + 1) / (double)(nseg));
    }
    std::cout << std::endl;
}

/*!
\brief Samples disc photon with \ref calmuobsarr, and according to photon destination writes photon info into the following two files:
- `disc_params.dat`: for photons arriving at infinity without entering the corona; similar with \ref griddisk_nt_sp_pol_writeparam.
- `sca_params.dat`: for photons entering the corona: the parameter file has a *header*, followed by \f$N\f$ segments where \f$N\f$ is the number of geodesic that enters the
    corona. The structure of the file: 
    \f{equation}a, r_{\rm out}, ctype, N_r, N_\theta, N_\phi,
    [r, \theta, \mu_\infty, w, r_{\rm hit},\mu_{\rm hit}, k_{\rm hit}^t, k_{\rm hit}^r, k_{\rm hit}^\theta,
    k_{\rm hit}^\phi, g, c_1, c_2, |WP|, {\rm cos}\Delta\psi, {\rm sin}\Delta\psi]\times N,
    \f}
    here *ctype* is the corona type, \f$r\f$ is the
    emission radius, \f$\theta\f$ is the local emission angle, \f$w\f$ is the statistical weight (same as in \ref griddisk_nt_sp_pol_writeparam).
    \f$r_{\rm hit}\f$ and \f$\mu_{\rm hit}\f$ are the coordinates where the photon enters the corona, \f$\boldsymbol{k}_{\rm hit}\f$ is the
    photon wave vector in Boyer-Lindquist frame while entering the corona, \f$g\equiv E_\infty/E_{\rm disc}\f$ is the redshift factor,
    \f$c_1\f$ and \f$c_2\f$ are coefficients of the rotation matrix from polarsation vector at the disc to 
    \f$k^\prime\f$ (see \ref nt::polrotang_disktok), \f$|WP|\f$ is the modulus of the Walker-Penrose constant, and \f$\Delta\psi\f$ is the angle between
    \f$k^\prime\f$ to the polarisation vector at infinity (see \ref nt::polrotang_ktoinf). The file is in binary format with double precision.
    
On *ctype*:
- 0: \ref zamosphere3d
- 1: \ref zamoslab3d
- 2: \ref zamosphere3d_truncated
- 3: \ref rotsphere3d
    
@param rin inner edge radius of the thin disc. If rin < rms then rin = rms.
@param rout outer boundary of thin disc
@param thetamax maximum local emission angle (see \ref calmuobsarr).
@param trid \ref tridgeo object; for different types of corona
@param nr number of radial bin on the disc plane
@param ntheta number of theta bin
@param nphi number of phi bin
*/
void grid3dcorona_sp_writeparam(const double rin, const double rout, const double thetamax, const tridgeo & trid, 
    const long nr, const long ntheta, const long nphi) {
    
    bool hit;
    int status, gtype;
    long indexnow;
    double knorm;
    double rinnow, dlogr, dr, one_over_dr, rnow, dphi, dtheta, gdisk, vphi, omega, area, weightnorm, weightnorm0;
    double x, sint, domega, gamma, one_over_ut;
    double rhit, muhit, wpnorm;
    double cospsi, sinpsi, c1, c2, cospsiinf, sinpsiinf;
    double d_coronatype;
    
    if (trid.name == "zamosphere3d") {
        gtype = 0;
    } else if (trid.name == "zamoslab3d") {
        gtype = 1;
    } else if (trid.name == "zamosphere3d_truncated") {
        gtype = 2;
    } else if (trid.name == "rotsphere3d") {
        gtype = 3;
    } else if (trid.name == "kepslab3d") {
        gtype = 5;
    } else {
        std::cerr << "Invalid corona type!" << std::endl;
        return;
    }
    d_coronatype = (double)(gtype);
    
    auto dsize = sizeof(double);
    
    // grid vectors
    std::vector<double> muobsarr, rr1arr, muplusarr, larr, qarr, kmuobssignarr;
    
    std::complex<double> wpnow;
    std::array<double, 4> on_k, bl_k, kmu_hit;
    std::vector<double> radius(nr), logradius(nr), theta(ntheta), phi(nphi), muarr(ntheta);
    
    std::ofstream disk_pfile("disk_params.dat", std::ios::out | std::ios::binary);
    std::ofstream sca_pfile("sca_params.dat", std::ios::out | std::ios::binary);
    
    // write params
    std::vector<double> input_params = {trid.a, rout, d_coronatype, 
        (double)(nr), (double)(ntheta), (double)(nphi)};
    disk_pfile.write((char*)&input_params[0], input_params.size() * dsize);
    sca_pfile.write((char*)&input_params[0], input_params.size() * dsize);
    
    /*
    std::cout << "trid.a = " << trid.a << std::endl;
    std::ofstream logfile("input.log", std::ios::out);
    logfile << "Corona type = " << trid.name << std::endl << 
        "rmin = " << trid.rmin << std::endl << "rmax = " << trid.rmax << std::endl;
    logfile << "a = " << trid.a << std::endl << "rout = " << rout << std::endl;
    logfile << "nr = " << nr << std::endl << "ntheta = " << ntheta << std::endl << "nphi = " << nphi << std::endl;
    logfile.close();
    */
    
    nt disk(trid.a, 1., 1.);
    sim5metric met;
    sim5tetrad tet;
    geoinf geo;
    
    rinnow = (rin < disk.rms) ? disk.rms : rin;
    
    dlogr = std::log(rout / rinnow) / (double)(nr);
    dr = std::exp(dlogr * 0.5);
    one_over_dr = 1. / dr;
    
    std::cout << "Calculating geodesics to infinity..." << std::endl;
    calmuobsarr(trid.a, rin, rout, thetamax, nr, ntheta, nphi, muobsarr, radius, theta, phi, rr1arr, muplusarr, larr, qarr, kmuobssignarr);
    std::cout << "Done!" << std::endl;
    std::transform(theta.cbegin(), theta.cend(), muarr.begin(), [](double x){return std::cos(x);});
    std::cout << "muarr = " << muarr << std::endl;
    
    dtheta = theta[1] - theta[0];
    dphi = phi[1] - phi[0];
    
    on_k[0] = 1.;
    
    std::cout << "Calculating geodesics to corona..." << std::endl;
    for (size_t ir = 0; ir < radius.size(); ++ir) {
        //std::cout << "ir = " << ir << std::endl;
        rnow = radius[ir];
        area = calds(trid.a, rnow * one_over_dr, rnow * dr);
        vphi = disk.vphi(rnow);
        gamma = 1. / std::sqrt(1. - vphi * vphi);
        
        x = std::sqrt(rnow);
        omega = 1. / (x * x * x + trid.a);
        
        one_over_ut = std::sqrt(1. - 2. / rnow + 4. * omega * trid.a / rnow - omega * omega * (rnow * rnow + 
            trid.a * trid.a + 2. * trid.a * trid.a / rnow));
        
        kerr_metric(trid.a, rnow, 0., &met);
        tetrad_azimuthal(&met, omega, &tet);
        
        weightnorm0 = area * gamma * one_over_ut;

        for (size_t it = 0; it < theta.size(); ++it) {
            sint = std::sin(theta[it]);
            on_k[2] = muarr[it];
            domega = std::abs(muarr[it]) * dphi * (std::cos(theta[it] - dtheta / 2.) - std::cos(theta[it] + dtheta / 2.));
            weightnorm = weightnorm0 * domega;
            
            for (size_t ip = 0; ip < phi.size(); ++ip) {
                //std::cout << "ir = " << ir << ", it = " << it << ", ip = " << ip << std::endl;
                indexnow = ir * theta.size() * phi.size() + it * phi.size() + ip;
                on_k[1] = sint * std::cos(phi[ip]);
                on_k[3] = sint * std::sin(phi[ip]);
                on2bl(on_k.data(), bl_k.data(), &tet);
                // normalize wave vector such that k_t = -1;
                knorm = -1. / (met.g00 * bl_k[0] + met.g03 * bl_k[3]);
                std::transform(bl_k.cbegin(), bl_k.cend(), bl_k.begin(), std::bind2nd(std::multiplies<double>(), knorm));
                gdisk = disk.gfactor(radius[ir], larr[indexnow]);
                
                geo = geoinf(trid.a, radius[ir], 0., bl_k);
                if (geo.l != geo.l || geo.q != geo.q) {
                    std::cout << "Inside grid3dcorona_sp; Checking geo.l:" << std::endl;
                    std::cout << "r = " << radius[ir] << ", mue = 0" << std::endl;
                    std::cout << "bl_k = " << bl_k << std::endl;
                    break;
                }
                
                disk.polrotang_disktok(radius[ir], theta[it], phi[ip], larr[indexnow], qarr[indexnow], c1, c2);
                disk.polrotang_ktoinf(kmuobssignarr[indexnow], larr[indexnow], qarr[indexnow], muobsarr[indexnow], cospsiinf, sinpsiinf);
                cospsi = c1 * cospsiinf - c2 * sinpsiinf;
                sinpsi = c1 * sinpsiinf + c2 * cospsiinf;
                wpnorm = std::sqrt(qarr[indexnow] + (larr[indexnow] - trid.a) * (larr[indexnow] - trid.a));
                
                hit = false;
                // If it is not possible for the photon to reach the corona
                if (!geo.ispossible(trid, radius[ir], 0., bl_k[1], bl_k[2])) {
                    hit = false;
                } else {
                    status = geo.findxtrid(trid, radius[ir], 0., bl_k[1], bl_k[2], rin, rhit, muhit, kmu_hit);
                    //catch(const char * error) {
                    //    std::cout << error << std::endl;
                    //}
                    
                    if (status != 1) {
                        hit = false;
                    } else {
                        //std::cout << "muobsarr[indexnow] = " <<  muobsarr[indexnow] << std::endl;
                        hit = true;
                        //kerr_metric(trid.a, rhit, muhit, &metnow);
                        //trid.calumu(rhit, muhit, umu, metnow);
                        //ghit = dotprod(umu.data(), kmu_hit.data(), &metnow);
                        //gfac = -1. * ghit * gdisk;
                        std::vector<double> sca_params = {radius[ir], theta[it], muobsarr[indexnow],
                            weightnorm, rhit, muhit, kmu_hit[0], kmu_hit[1], kmu_hit[2], kmu_hit[3], gdisk,
                            c1, c2, wpnorm, cospsi, sinpsi};
                        sca_pfile.write((char*)&sca_params[0], sca_params.size() * dsize);
                    }
                }
                
                if ((!hit) && muobsarr[indexnow] > -1.5) {
                    std::vector<double> disk_params = {radius[ir], theta[it], gdisk, muobsarr[indexnow], weightnorm,
                        cospsi, sinpsi};
                    disk_pfile.write((char*)&disk_params[0], disk_params.size() * dsize);
                }
                showprogress((double)(indexnow) / (double)(nr * ntheta * nphi));
            }
        }
    }
    std::cout << "Done!" << std::endl;
    disk_pfile.close();
    sca_pfile.close();
}

/*!
\brief With the `sca_params.dat` produced by \ref register_sp_tridsca, propagate the photons entering the corona until they arrive at infinity
/strike at the disc/enter the BH event horizon.

@param pol whether to include polarisation. If turned off, all photons are assumed to be unpolarised when being scattered by electron
@param wsample whether to weighted sample blackbody radiation
@param progress whether to print progress
@param parafile the path to `sca_params.dat`
@param m BH mass in solar mass
@param mdot mass accretion rate in Eddington rate (see \ref nt::nt(const double,const double,const double)).
@param fcol color correction factor
@param te electron temperature in \f$[\rm keV]\f$
@param mfp Thompson mean free path
@param dr step size in \f$f\f$ while propagating inside the corona
@param nphoton number of photon per geodesic
@param ofiles vector of \ref std::ofstream that point to output files; all files are in binary format and contain data in double precision.
    ofile[8] and ofiles[9] are written only if `pol` is turned on.
    - ofiles[0]: photon energy at infinity \f$E_\infty\f$
    - ofiles[1]: photon statistical weight \f$w^\prime\f$ (see \ref register_sp_nt)
    - ofiles[2]: 
        - inclination at infinity, if the photon finally arrives at infinity
        - disc radius, if the photon hits the disc;
    - ofiles[3]: \f$l\equiv L_z/E_\infty\f$, where \f$L_z\f$ is photon angular momentum about black hole spin axis
    - ofiles[4]: \f$\mathcal{Q}/E_\infty^2\f$, where \f$\mathcal{Q}\f$ is one constant of motion found by Carter 1968
    - ofiles[5]: sign of \f$k^r\f$
    - ofiles[6]: sign of \f$k^\theta\f$
    - ofiles[7]: only if `pol` is turned on; Stokes parameter Q \f$(q w^\prime)\f$, where \f$q\f$ is one normalised Stokes parametere
    - ofiles[8]: only if `pol` is turned on; Stokes parameter U \f$(u w^\prime)\f$, where \f$u\f$ is another normalised Stokes parametere
@param nscafile one \ref std::ofstream object that point to a binary file that contains photon scattering number in `int` format (note that
the actual size of `int` type depends on compiler and/or architect)
 */
void register_sp_tridsca(const bool pol, const bool wsample, const bool progress, 
    const std::string parafile, const double m, const double mdot, const double fcol,
    const double te, const double mfp, const double dr, const long nphoton, std::vector<std::ofstream> & ofiles_inf,
    std::ofstream & nscafile_inf, std::vector<std::ofstream> & ofiles_disc, std::ofstream & nscafile_disc) {
    
    double ecut = 10., ratio = 10., ennow, weightnow, theta1, dit;
    
    std::vector<double> params, theta, muarr, darkarr, poldegarr;
    double a, teff, tcol_emass, flimb, dtheta, weightnorm, poldegnow, K1, K2, ftheta, fphi, qfac, ufac, cos2psi, sin2psi;
    rdoublevec(params, parafile);
    long ntheta, n = 0, nparam = 16, it0, it1, nseg, offset = 6;
    long gtype;
    nseg = (params.size() - offset) / nparam;
    std::array<double, 4> kmunow;
    std::complex<double> wpnow;
    
    std::unique_ptr<tridgeo> tridptr(nullptr);
    std::vector<double> gparams;
    rdoublevec(gparams, "../gparams.dat");
    a = params[0];
    gtype = (long)(std::round(params[2]));
    std::cout << "gtype = " << gtype << std::endl;
    
    switch(gtype) {
        case 0 : {
            tridptr.reset(new zamosphere3d(a, gparams[0], gparams[1]));
            break;
        }
        case 1 : {
            tridptr.reset(new zamoslab3d(a, gparams[0], gparams[1], gparams[2]));
            break;
        }
        case 2 : {
            tridptr.reset(new zamosphere3d_truncated(a, gparams[0], gparams[1]));
            break;
        }
        case 3 : {
            tridptr.reset(new rotsphere3d(a, gparams[0], gparams[1], gparams[2]));
            break;
        }
        case 5 : {
            tridptr.reset(new kepslab3d(a, gparams[0], gparams[1], gparams[2]));
            break;
        }
        default : {
            std::cerr << "Invalid gtype value of " << gtype << std::endl;
        }
    }
    //zamosphere3d trid(a, gparams[0], gparams[1]);
    nt disk(a, m, mdot);
    //std::cout << "total photons = " << nphoton1 * nseg << std::endl;
    
    std::ofstream logfile("input.log", std::ios::out);
    logfile << "Corona type = " << tridptr->name << std::endl << 
        "rmin = " << tridptr->rmin << std::endl << "rmax = " << tridptr->rmax << std::endl;
    logfile << "te = " << te << std::endl << "mfp = " << mfp << std::endl;
    logfile << "pol = " << pol << std::endl << "a = " << tridptr->a << std::endl << "m = " << disk.m << std::endl
        << "mdot = " << disk.mdot << std::endl << "fcol = " << fcol << std::endl;
    logfile << "ntheta = " << ntheta << std::endl << "nphoton = " << nphoton << std::endl;
    logfile.close();
    
    // ntheta = (long)(params[4]);
    ntheta = 501;
    dtheta = 0.5 * M_PI / (double)(ntheta - 1);
    
    if (pol) {
        theta.resize(ntheta);
        std::generate(theta.begin(), theta.end(), [&n, dtheta]{return (double)(n++) * dtheta;});
        muarr.resize(ntheta);
        darkarr.resize(ntheta);
        std::transform(theta.cbegin(), theta.cend(), muarr.begin(), [](double x){return std::abs(std::cos(x));});
        poldegarr.resize(ntheta);
        calpoldeg(muarr, poldegarr);
        calchandra_darkening(muarr, darkarr);
    }
#pragma omp parallel
{
#pragma omp single
{
    std::mt19937_64 gen(std::random_device {}());
    for (long i = 0; i < nseg; ++i) {
        auto pos = i * nparam + offset;
        std::copy(params.cbegin() + pos + 6, params.cbegin() + pos + 10, kmunow.begin());
        
        teff = disk.teff_kev(params[pos]);
        //std::cout << "(T / Tin)^3 = " << teff / tin << std::endl;
        tcol_emass = fcol * teff / ME_KEV * params[pos + 10];;
        if (pol) {
            theta1 = (params[pos + 1] <= 0.5 * M_PI) ? params[pos + 1] : M_PI - params[pos + 1];
            dit = (theta1 - 0.5 * dtheta) / dtheta;
            it0 = (long)(std::floor(dit));
            it1 = it0 + 1;
            flimb = (darkarr[it0] * (theta1 - theta[it0]) + darkarr[it1] * (theta[it1] - theta1)) / dtheta;
            poldegnow = (poldegarr[it0] * (theta1 - theta[it0]) + poldegarr[it1] * (theta[it1] - theta1)) / dtheta;
            //std::cout << "param[i * nparam + 1] = " << params[i * nparam + 1] << std::endl;
            //std::cout << "theta[it] = "  << theta[it] << std::endl;
            
            K1 = -1. * params[pos + 13] * params[pos + 12];
            K2 = -1. * params[pos + 13] * params[pos + 11];
            wpnow.real(K1);
            wpnow.imag(K2);
            ftheta = -1. * params[pos + 15];
            fphi = params[pos + 14];
            sin2psi = 2. * ftheta * fphi;
            cos2psi = 2. * ftheta * ftheta - 1.;
            qfac = poldegnow * cos2psi;
            ufac = poldegnow * sin2psi;
            
        } else {
            flimb = 1.;
            qfac = 0.;
            ufac = 0.;
        }
        weightnorm = params[pos + 3] * disk.m * disk.m  * flimb * teff * teff * teff / fcol / (double)(nphoton);
        
        for (long ip = 0; ip < nphoton; ++ip) {
            superphoton sp;
            if (wsample) {
                sample_bb1(tcol_emass, gen, ecut, ratio, ennow, weightnow);
                sp.en = ennow;
            } else {
                sp.en = sample_bb(tcol_emass, gen);
                weightnow = 1.;
            }
                
            sp.kmu = kmunow;
            sp.rnow = params[pos + 4];
            sp.munow = params[pos + 5];
            sp.poldeg = poldegnow;
            sp.wp = wpnow;
            sp.nscatter = 0;
            sp.weight = weightnow;
        }
        if (progress) {
            showprogress((double)(i) / (double)(nseg));
        }
   }
    std::cout << std::endl;
}
}
}

