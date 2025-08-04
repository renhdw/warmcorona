//! \file diskobj.h
//! This fine includes definitions of virtual `diskobj` class that can be further derived to store information
//! of different types of accretion disc. One derived class, `nt` (for Novikov-Thorne disc),
//! is also included.

#ifndef _DISKOBJ_H
#define _DISKOBJ_H
#include <string>
#include <array>
#include <complex>
#include <random>

//! \brief The abstract base class for all disc objects.
class diskobj {
public:
    //! disc name
    std::string name; 
    //! BH spin
    double a; 
    //! BH mass in solar mass
    double m;
    //! Mass accretion rate in Eddington accretion rate; \f$ Mdot = 2.225475942\times 10^{18} m \  mdot\ \rm [g~s^{-1}]\f$.
    double mdot;    
    //! Radius of innermost stable orbit, in \f$\rm [GM\ c^{-2}]\f$.
    double rms;

    //! virtual destructor
    virtual ~diskobj(){};
    //diskobj(const double aa, const double mm, const double mmdot);
    // member functions
    //! Virtual method for calculating redshift factor
    virtual double gfactor(const double r, const double lambda) const = 0;
    //! Virtual method for calculating local flux density
    virtual double flux(const double r) const = 0;
    //! Virtual method for calculating effective temperature in \f$[\rm K]\f$
    virtual double teff(const double r) const = 0;
    //! Virtual method for calculating effective temperature in \f$[\rm keV]\f$
    virtual double teff_kev(const double r) const = 0;
    //! Virtual method for calculating linear circular velocity of the disc fluid
    virtual double vphi(const double r) const = 0;
    //! Virtual method for calculating photon emission angle in disc fluid rest frame given constants of motion.
    virtual void angles(const double r, const double l, const double q, const double signr, std::array<double, 4> & cosang) = 0;
    //! Virtual method for calculating Walker-Penrose constant given constants of motion.
    virtual std::complex<double> calwp(const double r, const double l, const double q) = 0;
    //! Virtual method for calculating Walker-Penrose constant given photon emission angle in disc fluid rest frame
    virtual std::complex<double> calwp_angle(const double r, const double theta, const double phi, const double l) = 0;
};

std::ostream & operator << (std::ostream & out, const diskobj & disk);

//! Novikov-Thorne disc. We assume that the disc extends down to ISCO and their is zero torque at the inner edge of the disc.
class nt : public diskobj {
public:
    //! Constructor.
    nt(){};
    //! Destructor.
    ~nt(){};
    //! Constructor.
    nt(const double aa, const double mm, const double mmdot);
    //! Calculates redshift factor \f$g=E_{\infty}/E_{\rm disc}\f$.
    double gfactor(const double r, const double lambda) const;
    //! Calculates flux density in the local rest frame of the disc fluid; in [\f$\rm erg~s^{-1}~cm^{-2}\f$]; cf., PT74 Eq. 11b.
    double flux(const double r) const;
    //! Calculates effective temperature, in \f$\rm [K]\f$.
    double teff(const double r) const;
    //! Calculates effective temperature, in \f$[\rm keV]\f$.
    double teff_kev(const double r) const;
    //! Calculates velocity along \f$\phi\f$- direction; measured by a ZAMO observer.
    double vphi(const double r) const;
    //! Calculates photon emission angle in disc fluid rest frame given constants of motion.
    void angles(const double r, const double l, const double q, const double signr, std::array<double, 4> & cosang);
    //! Calculates Walker-Penrose constant given constants of motion.
    std::complex<double> calwp(const double r, const double l, const double q);
    //! Calculates Walker-Penrose constant given photon emission angle in disc fluid rest frame.
    std::complex<double> calwp_angle(const double r, const double theta, const double phi, const double l);
    //! Walker-Penrose constant given photon emission angle in disc fluid rest frame.
    std::complex<double> calwp_angle(const double r, const double theta, const double phi, const double l, const double polang);
    //! The rotation angle from polarisation vector on the disc to an unit vector related with the Walker-Penrose constant, given constants of motion.
    double polrotang_disktok(const double r, const double theta, const double phi, const double l, const double q);
    //! The coefficients from polarisation vector on the disc to an unit vector related with the Walker-Penrose constant, given constants of motion.
    void polrotang_disktok(const double r, const double theta, const double phi, const double l, const double q, double & c1, double & c2);
    //! Cosine and sinusoidal of the rotation angle from an unit vector to the polarisation vector at infinity.
    void polrotang_ktoinf(const double ktheta, const double l, const double q, const double muobs, double & cospsi, double & sinpsi);
public:
    // some useful variables;
    double x0, x1, x2, x3, f0base, f1base, f2base, f3base; /*!< private variables for calculating flux density.*/
};

//! The thin disc with inner edge radius larger than rms and zero torque at inner disc
class truncated_nt : public nt {
public:
    //! disc inner edge radius, in \f$[\rm GM~c^{-2}]\f$
    double rtr;
    //! constructor
    truncated_nt(const double aa, const double mm, const double mmdot, const double rtrr);
    //! constructor
    truncated_nt(){};
    //! destructor
    ~truncated_nt(){};
    //! testunit
    void selftest(const double r);
};
    
#endif
