//! \file superphoton.h
//! Defines three classes: \ref superphoton, \ref superphoton_sr, and \ref superphoton_sr_pol.
#ifndef _SUPERPHOTON_H
#define _SUPERPHOTON_H
#include <vector>
#include <array>
#include <cmath>
#include <random>
#include <complex>
#include <iostream>

/*
//! pack of superphoton
class superphotonpack {
public:
    double weight = 1., re, thetae, phie;
    std::vector<double> scattermuarr, loboostarr, inloboostarr, framegfacarr;
    std::array<double, 4> kmu;
    std::array<double, 3> xmu;
    
    ~superphotonpack(){};
    superphotonpack(){};
    superphotonpack(const double rre, const double ttheta, const double pphi) : re(rre),
        thetae(ttheta) {};
};
*/

/*! \brief Defines the `superphoton` class.
For details of the `superphoton` scheme see e.g., Dolence+2009 (http://adsabs.harvard.edu/abs/2009ApJS..184..387D).
 Compared with ``ordinary`` photon, superphoton also has its statistical weight and can ``split`` into multiple superphotons.
 */
class superphoton {
public:
    //! Superphoton energy at infinity; constant of motion
    double en;
    //! Statistical weight; has the physical meaning of superphoton generation rate per unit time in distant observer's frame.
    double weight;
    //! Polarisation degree; constant of motion
    double poldeg;
    //! Boyer-Lindquist coordinate \f$r\f$
	double rnow;
    //! Cosine of Boyer-Lindquist coordinate \f$\theta\f$
    double munow;
    //! Photon wave vector in BL frame
    std::array<double, 4> kmu;
    //! Number of scattering
    int nscatter;
    //! Walker-Penrose constant
    std::complex<double> wp;
    
    //! Destructor
    ~superphoton(){};
    //! Constructor without any argument
    superphoton(){};
    //! Constructor with photon energy
    superphoton(const double enn) : en(enn), weight(1.), nscatter(0) {};
    //! Copy constructor. If this function is not explicitly defined, the there will be problem passing argument in \ref tridsca_sp in parallel computation.
    superphoton(const superphoton & sp) : en(sp.en), weight(sp.weight), poldeg(sp.poldeg), 
        rnow(sp.rnow), munow(sp.munow), kmu(sp.kmu), nscatter(sp.nscatter), wp(sp.wp) {};
    // {std::cout << "sp0.poldeg = " << sp.poldeg << ", poldeg = " << poldeg << std::endl;};
    //! Overloaded assignment operator
    superphoton & operator = (const superphoton &);
};

/*!
\brief The superphoton class in flat spacetime. This class is also without polarisation transportation.
 */
class superphoton_sr {
public:
    //! Superphoton dimensionless energy
    double en;
    //! Statistical weight; has the physical meaning of superphoton generation rate per unit time
    double weight;
    //! Photon position, in Cartesian coordinate
    std::array<double, 3> xmu;
    //! Photon wave vector
    std::array<double, 4> kmu;
    //! Number of scattering
    int nscatter;
    
    //! Destructor
    ~superphoton_sr(){};
    //! Constructor without any argument
    superphoton_sr(){};
    //! Constructor with photon energy
    superphoton_sr(const double enn) : en(enn), weight(1.), nscatter(0) {};
};

//! Superphoton class in flat spacetime, but with polarisation info as compared with \ref superphoton_sr.
class superphoton_sr_pol {
public:
    //! superphoton dimensionless energy
    double en;
    //! pol. degree
    double poldeg;
    //! Statistical weight; has the physical meaning of superphoton generation rate per unit time
    double weight;
    //! Photon position, in Cartesian coordinate
    std::array<double, 3> xmu;
    //! Photon wave vector
    std::array<double, 4> kmu;
    //! Polarisation vector
    std::array<double, 4> fmu;
    //! Number of scattering
    int nscatter;
    
    //! Destructor
    ~superphoton_sr_pol(){};
    //! Constructor without any arguments
    superphoton_sr_pol(){};
    //! Constructor with energy and pol. degree
    superphoton_sr_pol(const double enn, const double poldegg) : en(enn), poldeg(poldegg), weight(1.), nscatter(0) {};
};
    
#endif
