//! \file detector.cpp
//! Implement ABC `detector` and other classes derived from it.
#include <algorithm>
#include <cmath>
#include <iomanip>
#include "detector.h"
#include "utils.h"

//! @param nee size of energy bin. Linear space if positive; otherwise logarithm space.
//! @param eminn lower boundary of energy bin 
//! @param emaxx upper boundary of energy bin 
detector::detector(const long nee, const double eminn, const double emaxx) {
    double de1;
    long n = 0;
    
    emin = eminn;
    emax = emaxx;
    
    linear = (nee > 0);
    std::cout << "linear = " << linear << std::endl;
    ne = std::abs(nee);
    
    en.resize(ne);
    
    if (linear) {
        de = (emax - emin) / (double)(ne);
        de1 = de;
        n = 0;
        std::generate(en.begin(), en.end(), [&n, eminn, de1]{return (double)(n++) * de1 + de1 / 2. + eminn;});
    } else {
        double logemin, logemax;
        std::vector<double> logen(ne);
        logemin = std::log(emin);
        logemax = std::log(emax);
        de = (logemax - logemin) / (double)(ne);
        de1 = de;
        n = 0;
        std::generate(logen.begin(), logen.end(), [&n, logemin, de1]{return (double)(n++) * de1 + de1 / 2. + logemin;});
        std::transform(logen.cbegin(), logen.cend(), en.begin(), [](double x){return std::exp(x);});
    }
}

/*
//! For detail of the `photon package` scheme, see Schnittman & Krolik 2013.
//! @param muobs \f$\mu_{\rm obs} = {\rm cos}(\theta_{\rm obs})\f$, where \f$\theta_{\rm obs}\f$ is the photon inclination angle at infinity
//! @param pack `gammapack` object
void singledet::receive(const double muobs, const gammapack & pack) {
    double gfac = -1. / pack.k_t, c0, c1, index, floor, fluxtemp;
    size_t i0, i1;
    // interpolate energy grid from photon pack to receive grid
    if (muobs > mumin && muobs <= mumax) {
        for (size_t j = 0; j < ne; ++j) {
            index = (en[j] - pack.emin * gfac) / pack.de / gfac;
            floor = std::floor(index);
            i0 = (size_t)(floor);
            i1 = i0 + 1;
            if (i1 >= pack.en.size()) {
                std::cerr << "energy grid out of range!" << std::endl;
            }
            c1 = index - floor;
            c0 = (1. - c1);
            fluxtemp = c0 * pack.flux[i0] + c1 * pack.flux[i1];
            flux[j] += fluxtemp;
        }
    }
}

//! For detail of the `photon package` scheme, see Schnittman & Krolik 2013.
//! @param muobs \f$\mu_{\rm obs} = {\rm cos}(\theta_{\rm obs})\f$, where \f$\theta_{\rm obs}\f$ is the photon inclination angle at infinity
//! @param pack `gammapack` object
void singledet::receive_align(const double muobs, const gammapack & pack) {
    if (ne != pack.en.size()) {
        std::cerr << "size not consistent!" << std::endl;
        return;
    }
    if (muobs >= mumin && muobs <= mumax) {
        std::transform(flux.cbegin(), flux.cend(), pack.flux.cbegin(), flux.begin(), std::plus<double>());
    }
}
*/

//! The factor is \f$4\pi/d\Omega\f$.
void singledet::norm() {
    double fourpi_over_domega = 2. / (mumax - mumin);
    std::transform(flux.begin(), flux.end(), flux.begin(), std::bind2nd(std::multiplies<double>(), fourpi_over_domega));
}

//! After the normalisation, the flux unit is converted from \f$[\rm photons\ s^{-1}\ bin^{-1}]\f$ to \f$[\rm photons\ s^{-1}\ keV^{-1}]\f$.
void singledet::norm_flux() {
    if (linear) {
        std::transform(flux.begin(), flux.end(), flux.begin(), std::bind2nd(std::multiplies<double>(), 1. / de));
    } else {
        double dloge = std::log10(en[1] / en[0]);
        double sqrt_Eratio = std::pow(10., 0.5 * dloge);
        sqrt_Eratio = sqrt_Eratio - 1. / sqrt_Eratio;
        
        std::vector<double> dearr(en.size());
        std::transform(en.cbegin(), en.cend(), dearr.begin(), std::bind2nd(std::multiplies<double>(), sqrt_Eratio));
        std::transform(flux.begin(), flux.end(), dearr.begin(), flux.begin(), std::divides<double>());
    }
}

//! After the normalisation, the flux unit is converted from \f$[\rm photons\ s^{-1}\ bin^{-1}]\f$ to \f$[\rm photons\ s^{-1}\ keV^{-1}]\f$.
void singledet_sp::norm_flux() {
    if (linear) {
        std::transform(flux.begin(), flux.end(), flux.begin(), std::bind2nd(std::multiplies<double>(), 1. / de));
    } else {
        double dloge = std::log10(en[1] / en[0]);
        double sqrt_Eratio = std::pow(10., 0.5 * dloge);
        sqrt_Eratio = sqrt_Eratio - 1. / sqrt_Eratio;
        
        std::vector<double> dearr(en.size());
        std::transform(en.cbegin(), en.cend(), dearr.begin(), std::bind2nd(std::multiplies<double>(), sqrt_Eratio));
        std::transform(flux.begin(), flux.end(), dearr.begin(), flux.begin(), std::divides<double>());
    }
}

void singledet::write() {
    wdoublevec(en, "en.dat");
    wdoublevec(flux, "flux.dat");
}

//! @param filename filename
void singledet::writename(std::string filename) {
    wdoublevec(en, "en.dat");
    wdoublevec(flux, filename);
}

void singledet::print() const {
    std::cout << std::setprecision(4) << "Covering mu in (" << mumin << ", " << mumax << "]" << std::endl;
}

void multidet::init() {
    binsize = 1. / (double)(nbin);
    bincen.resize(nbin);
    long n = 0;
    std::generate(bincen.begin(), bincen.end(), [&n, binsize=binsize] {return (double)(n++) * binsize + binsize / 2.;});
    dets.resize(nbin);
    
    for (size_t ibin = 0; ibin < nbin; ++ibin)
        dets[ibin] = singledet(bincen[ibin] - binsize / 2., bincen[ibin] + binsize / 2., ne, emin, emax);
}


/*
//! @param muobs cosine of photon inclination angle at infinity
//! @param pack the photon package object
void multidet::receive(const double muobs, const gammapack & pack) {
    for (size_t ibin = 0; ibin < nbin; ++ibin) {
        if (dets[ibin].detected(muobs))
            dets[ibin].receive(muobs, pack);
    }
}

//! @param muobs \f$\mu_{\rm obs} = {\rm cos}(\theta_{\rm obs})\f$, where \f$\theta_{\rm obs}\f$ is the photon inclination angle at infinity
//! @param pack `gammapack` object
void multidet::receive_align(const double muobs, const gammapack & pack) {
    for (size_t ibin = 0; ibin < nbin; ++ibin) {
        if (dets[ibin].detected(muobs))
            dets[ibin].receive_align(muobs, pack);
    }
}
*/


void multidet::multiply(const double fac) {
    for (size_t ibin = 0; ibin < nbin; ++ibin) {
        dets[ibin].multiply(fac);
    }
}
    
void multidet::norm() {
    for (size_t ibin = 0; ibin < nbin; ++ibin) {
        dets[ibin].norm();
    }
}

void multidet::write() {
    wdoublevec(bincen, "bincenmu.dat");
    for (size_t ibin = 0; ibin < nbin; ++ibin) {
        dets[ibin].writename("flux_" + std::to_string(ibin) + ".dat");
    }
}

void multidet::print() const {
    std::cout << nbin << " detectors in total:" << std::endl;
    for (size_t i = 0; i < nbin; ++i) {
        std::cout << i << "-th detector: ";
        dets[i].print();
    }
}

//! This class has additional arrays for Stokes parameters and polarisation degree/angle compared with `singledet` class
void poldet::init() {
    qflux.resize(ne);
    uflux.resize(ne);
    poldeg.resize(ne);
    polang.resize(ne);
}

void poldet::norm() {
    double fourpi_over_domega = 2. / (mumax - mumin);
    double polflux, xiq, ang;
    for (size_t ie = 0; ie < ne; ++ie) {
        polflux = std::sqrt(qflux[ie] * qflux[ie] + uflux[ie] * uflux[ie]);
        poldeg[ie] = polflux / flux[ie];
        xiq = qflux[ie] / polflux;
        ang = std::acos(xiq);
        polang[ie] = (uflux[ie] >= 0. ? ang / 2. : M_PI - ang / 2.);
        flux[ie] *= fourpi_over_domega;
    }
}

void poldet::print() const {
    std::cout << "Polarimeter; " << std::setprecision(4) << "covering mu in (" << mumin << ", " << mumax << "]" << std::endl;
}

void poldet::write() {
    wdoublevec(en, "en.dat");
    wdoublevec(flux, "flux.dat");
    wdoublevec(qflux, "qflux.dat");
    wdoublevec(uflux, "uflux.dat");
    
    wdoublevec(poldeg, "poldeg.dat");
    wdoublevec(polang, "polang.dat");
}

//! @param head if `head == abc`, then the flux file will be `flux_abc.dat`
void poldet::writename(const std::string head) {
    wdoublevec(en, "en.dat");
    wdoublevec(flux, "flux_" + head + ".dat");
    wdoublevec(poldeg, "poldeg_" + head + ".dat");
    wdoublevec(polang, "polang_" + head + ".dat");
}

/*
//! Receive a `photon package` at \f$\mu_{\rm obs}\f$. For detail of the `photon package` scheme, see Schnittman & Krolik 2013.
//! This member function works if the energy bins of `pack` and the detector agree (but it does not check this!).
//! Note that we no longer use `photon package`.
//! @param muobs \f$\mu_{\rm obs} = {\rm cos}(\theta_{\rm obs})\f$, where \f$\theta_{\rm obs}\f$ is the photon inclination angle at infinity
//! @param pack `gammapack` object
void poldet::receive_align(const double muobs, const gammapack & pack) {
    if (pack.pol == false) {
        std::cerr << "Photon pack without polarization info!" << std::endl;
        return;
    }
    
    if (detected(muobs)) {
        std::transform(flux.cbegin(), flux.cend(), pack.flux.cbegin(), flux.begin(), std::plus<double>());
        std::transform(qflux.cbegin(), qflux.cend(), pack.qflux.cbegin(), qflux.begin(), std::plus<double>());
        std::transform(uflux.cbegin(), uflux.cend(), pack.uflux.cbegin(), uflux.begin(), std::plus<double>());
    }
}
*/

//! For details of the `superphoton` scheme, see Dolence+2009. The `superphoton` class is defined 
//! and implemented in `superphoton.h` and `superphoton.cpp`
//! @param sparr array of superphotons
void singledet_sp::receive_sp(const std::vector<superphoton> & sparr) {
    std::vector<double> enarr, warr;
    for (size_t i = 0; i < sparr.size(); ++i) {
        if (detected(sparr[i].munow)) {
            enarr.push_back(sparr[i].en);
            warr.push_back(sparr[i].weight);
        }
    }
    calflux(enarr, warr);
    //hist_flux(linear, enarr, warr, en, flux);
}

//! The factor is \f$4\pi/d\Omega\f$.
void singledet_sp::norm() {
    double fourpi_over_domega = 2. / (mumax - mumin);
    std::transform(flux.begin(), flux.end(), flux.begin(), std::bind2nd(std::multiplies<double>(), fourpi_over_domega));
}

void singledet_sp::write() {
    wdoublevec(en, "en.dat");
    wdoublevec(flux, "flux.dat");
}

//! @param filename flux filename
void singledet_sp::writename(std::string filename) {
    wdoublevec(en, "en.dat");
    wdoublevec(flux, filename);
}
    
//! Basically a simple histogram for array that are evenly spaced on either linear or log scale.
//! Currently not optimised; should find a faster algorithm/library.
//! The flux array has the unit of \f$[\rm counts\ s^{-1}\ bin^{-1}]\f$.
//! This function has a boolen argument as the condition for the photon to be detected. One can overload function with
//! specified condition.
//! @param enarr array of photon energy
//! @param warr array of photon weight
//! @param cond bool array; whether the corresponding (super)photon is detected
void singledet_sp::calflux(const std::vector<double> & enarr, const std::vector<double> & warr, const std::vector<bool> & cond) {
    // find out the index
    long index;
    if (linear) {
        for (size_t it = 0; it < enarr.size(); ++it) {
            index = (long) (std::floor((enarr[it] - emin) / de));
            if (cond[it] && index >= 0 && (size_t)(index) < flux.size()) {
                flux[index] += warr[it];
            }
        }
        
    } else {
        double dloge = std::log10(en[1] / en[0]);
        double logemin = std::log10(en[0]) - 0.5 * dloge;
        
        for (size_t it = 0; it < enarr.size(); ++it) {
            index = (size_t) (std::floor((std::log10(enarr[it]) - logemin) / dloge));
            if (cond[it] && index >= 0 && (size_t)(index) < flux.size()) {
                flux[index] += warr[it];
            }
        }
    }
}

//! The condition for detection is \f$\mu\in(\mu_{\rm min}, \mu_{\rm max}]\f$.
//! @param enarr array of photon energy
//! @param warr array of photon weight
//! @param muobsarr array of \f$\mu_{\rm obs}\f$
void singledet_sp::calflux(const std::vector<double> & enarr, const std::vector<double> & warr, const std::vector<double> & muobsarr) {
    std::vector<bool> cond(muobsarr.size());
    for (size_t it = 0; it < muobsarr.size(); ++it)
        cond[it] = detected(muobsarr[it]);
    calflux(enarr, warr, cond);
}

//! The condition is that the number of scattering equals `nsca` if `nsca >= 0`, and the number of scattering is no less 
//!     than `nsca` if `nsca < 0`.
//! @param enarr array of photon energy
//! @param warr array of photon weight
//! @param nscaarr array of scattering number
//! @param nsca the desired number of scattering
void singledet_sp::calflux(const std::vector<double> & enarr, const std::vector<double> & warr, const std::vector<int> & nscaarr, const int nsca) {
    std::vector<bool> cond(nscaarr.size());
    for (size_t it = 0; it < nscaarr.size(); ++it) {
        if (nsca >= 0) {
            cond[it] = (nscaarr[it] == nsca);
        } else {
            cond[it] = (nscaarr[it] >= -1 * nsca);
        }
    }
    calflux(enarr, warr, cond);
}

//! The condition for detection is \f$\mu\in(\mu_{\rm min}, \mu_{\rm max}]\f$ and number of scattering == nsca at the same time.
//! @param enarr array of photon energy
//! @param warr array of photon weight
//! @param muobsarr array of \f$\mu_{\rm obs}\f$
//! @param nscaarr array of scattering number
//! @param nsca the desired number of scattering
void singledet_sp::calflux(const std::vector<double> & enarr, const std::vector<double> & warr, const std::vector<double> & muobsarr, 
    const std::vector<int> & nscaarr, const int nsca) {
    
    std::vector<bool> cond(muobsarr.size());
    for (size_t it = 0; it < muobsarr.size(); ++it) {
        if (nsca >= 0) {
            cond[it] = detected(muobsarr[it]) && (nscaarr[it] == nsca);
        } else {
            cond[it] = detected(muobsarr[it]) && (nscaarr[it] >= -1 * nsca);
        }
    }
    calflux(enarr, warr, cond);
}

//! @param enarr array of photon energy
//! @param warr array of photon weight
void singledet_sp::calflux(const std::vector<double> & enarr, const std::vector<double> & warr) {
    std::vector<bool> cond(enarr.size(), true);
    calflux(enarr, warr, cond);
}
