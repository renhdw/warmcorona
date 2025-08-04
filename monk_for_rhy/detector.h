//! \file detector.h
#ifndef _DETECTOR_H
#define _DETECTOR_H
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include "superphoton.h"

//! Abstract base class for all detector classes.
class detector {
public:
    //! size of energy array
    size_t ne; 
    //! bisize of energy bin
    double de; 
    //! lower boundary of energy bin
    double emin; 
    //! lower boundary of energy bin
    double emax; 
    //! whether the energy bin is in linear or logarithm scale
    bool linear; 
    //! energy array
    std::vector<double> en; 
    
    //! (virtual) destructor
    virtual ~detector(){}; 
    //! Constructor
    detector(){}; 
    //! Constructor
    detector(const long ne, const double, const double);
    
    /*
    //! virtual method for detect a photon package with \b different energy grid
    virtual void receive(const double muobs, const gammapack & pack) = 0;
    //! virtual method for detect a photon package with \b identical energy grid
    virtual void receive_align(const double muobs, const gammapack & pack) = 0;
    //! virtual method for multiplying the flux array with a given factor
    */
    virtual void multiply(const double) = 0;
    //! virtual method for writing arrays to disc
    virtual void write() = 0;
    //! virtual method, to see if the photon is detected
    virtual bool detected(const double) = 0;
    //! virtual method to normalise the flux with respect to solid angle
    virtual void norm() = 0;
    //! virtual method to normalise the flux with respect to energy bin width
    virtual void norm_flux() = 0;
    //! virtual method to print information to stdout
    virtual void print() const = 0;
};

//! Single detector for `photon package`. This detector detects photons in a predefined inclination bin.
class singledet : public detector {
public:
    //! lower boundary of \f$\mu\equiv{\rm cos}\theta_{\rm obs}\f$
    double mumin;
    //! upper boundary of \f$\mu\equiv{\rm cos}\theta_{\rm obs}\f$
    double mumax;
    //! flux array
    std::vector<double> flux; 
    
    //! Destructor
    ~singledet(){}; 
    //! Constructor
    singledet(){}; 
    //! Constructor
    singledet(const double tmin, const double tmax, const long nee, const double eminn, const double emaxx) :
        detector(nee, eminn, emaxx), mumin(tmin), mumax(tmax) {flux.resize(ne);}; 
        
    /*
    //! Receive a `photon package` with inclination at infinity \f$\mu_{\rm obs}\f$, when the `photon package` has \b different energy grid with this object.
    void receive(const double muobs, const gammapack & pack); 
    //! Receive a `photon package` with inclination at infinity \f$\mu_{\rm obs}\f$, when the `photon package` has \b identical energy grid with this object.
    void receive_align(const double muobs, const gammapack & pack);
    //! Multiply flux with a given factor.
    */
    void multiply(const double fac) {
        std::transform(flux.begin(), flux.end(), flux.begin(), std::bind2nd(std::multiplies<double>(), fac));} 
    //! Write files to disc with default filename
    void write(); 
    //! Write files to disc with specified filename
    void writename(std::string filename); 
    //! Check if the photon is detected
    bool detected(const double muobs) { return (muobs > mumin && muobs <= mumax);}; 
    //! Normalise the spectrum with respect to solid angle.
    void norm();
    //! Normalise the spectrum with respect to energy bin width.
    void norm_flux();  
    //! print inclination bin range to stdout.
    void print() const;
};

//! single detector for superphoton
class singledet_sp : public detector {
public:
    //! lower boundary of \f$\mu\equiv{\rm cos}\theta_{\rm obs}\f$
    double mumin; 
    //! upper boundary of \f$\mu\equiv{\rm cos}\theta_{\rm obs}\f$
    double mumax; 
    //! flux array
    std::vector<double> flux; 
    
    //! Destructor
    ~singledet_sp(){};
    //! Constructor 
    singledet_sp(){};
    //! Constructor 
    singledet_sp(const double mumin, const double mumax, const unsigned long nee, const double eminn, const double emaxx) :
        detector(nee, eminn, emaxx), mumin(mumin), mumax(mumax) {flux.resize(ne);}; 
        
    //! Receive a `superphoton`.
    void receive_sp(const std::vector<superphoton> & sparr); 
    //! Calculate flux of photons with given condition. 
    void calflux(const std::vector<double> & enarr, const std::vector<double> & weightarr, const std::vector<bool> & cond);
    //! Calculate flux of photons in the inclination bin of this object.
    void calflux(const std::vector<double> & enarr, const std::vector<double> & weightarr, const std::vector<double> & muobsarr);
    //! Calculate flux of photons with given numbers of scattering.
    void calflux(const std::vector<double> & enarr, const std::vector<double> & warr, const std::vector<int> & nscaarr, const int nsca);
    //! Calculate flux of photons in the inclination of this object and with given nubmer of scattering at the same time.
    void calflux(const std::vector<double> & enarr, const std::vector<double> & warr, const std::vector<double> & muobsarr, 
        const std::vector<int> & nscaarr, const int nsca);
    //! Calculate flux of all photons.
    void calflux(const std::vector<double> & enarr, const std::vector<double> & weightarr);
    //! Multiply the flux with a given factor
    void multiply(const double fac) {
        std::transform(flux.begin(), flux.end(), flux.begin(), std::bind2nd(std::multiplies<double>(), fac));} 
    //! Write files to disc with default filename
    void write(); 
    //! Write files to disc with specified filenames
    void writename(std::string filename); 
    //! Check if the photon is detected
    bool detected(const double muobs) { return (muobs > mumin && muobs <= mumax);}; 
    //! Normalise the spectrum with respect to solid angle
    void norm();
    //! Normalise the spectrum with respect to energy bin width.
    void norm_flux();
    //! Normalise the spectrum with respect to solid angle
    void print() const {};
};

//! An set of several `singledet` with seamless coverage of half sky (i.e., with inclination from \f$0\f$ to \f$\pi/2\f$).
class multidet : public detector {
public:
    //! number of inclination bin
    unsigned long nbin; 
    //! size of inclination bin
    double binsize;
    //! array of inclination bin center
    std::vector<double> bincen;
    //! array `singledet`
    std::vector<singledet> dets;
    
    //! Destructor
    ~multidet(){};
    //! Constructor
    multidet(){};
    //! Constructor
    multidet(unsigned const long nbinn, const unsigned long nee, const double eminn, const double emaxx) :
        detector(nee, eminn, emaxx), nbin(nbinn) {init();};
        
    //! Init the object. Find out the boundaries of the single detectors, construct them, and save them to array.
    void init();
    //! Multiply flux for each single detector.
    void multiply(const double fac);
    //! Detect the photon for each single detector
    bool detected(const double muobs) { return (muobs > 0. && muobs <= 1.);};
    //! For each single detector write files to disc.
    void write();
    //! For each single detector normalise the flux with respect to the corresponding solid angle
    void norm();
    //! For each single detector print inclination range
    void print() const;
};

//! Polarisation detector.
class poldet : public singledet {
public:
    //! Destructor
    ~poldet(){};
    //! Constructor
    poldet(){};
    //! Constructor
    poldet(const double tmin, const double tmax, const unsigned long nee, const double eminn, const double emaxx) :
        singledet(tmin, tmax, nee, eminn, emaxx){init();};
    
    //! Flux of Stokes parameters \f$Q\f$
    std::vector<double> qflux;
    //! Flux of Stokes parameters \f$U\f$
    std::vector<double> uflux;
    //! Polarisation degree array
    std::vector<double> poldeg;
    //! Polarisation angle array
    std::vector<double> polang;
        
    //! Init the object.
    void init();
    //! Normalise spectrum with respect to solid angle
    void norm();
    //! Print info to stdout
    void print() const;
    //! Write files to disc with default filenames
    void write();
    //! Write files to disc with specified filenames
    void writename(const std::string head);
};

#endif
