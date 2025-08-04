//! \file tridgeo.h
//! Here we define an ABC `tridgeo` as a base class for classes of all kinds of 3D geometries.
#ifndef _TRIDGEO_H
#define _TRIDGEO_H
#include <string>
#include <array>
#include "sim5lib.h"

/*!
An abstract base class for all kinds of 3D corona. We expect that a general corona, regardless of its shape,
would be able to tell us the following: 1) if a photon is inside the corona given position; 2) the four-velocity of
corona fluid anywhere in the corona; 3) the local tetrad attached to corona fluid; 4) the size. Therefore we define the corresponding
virtual functions that will be implemented by the derived classes. 
The corona should have finite size, such that it has finite lower and upper bounds in \f$r\f$ and \f$\mu\f$.
These can be used to tell if the photon can reach the corona given current position and wave vector.
*/
class tridgeo {
public:
    //! geometry name
    std::string name;
    //! BH spin
    double a;
    //! Minimum \f$r\f$ of corona
    double rmin;
    //! Maximum \f$r\f$ of corona
    double rmax;
    //! Minimum \f$\mu\f$
    double mumin;
    //! Maximum \f$\mu\f$
    double mumax;
    //! Tell if the photon is inside
    virtual bool inside(const double r, const double mu) const = 0;
    //! Calculate the tetrad attached to the corona fluid
    virtual sim5tetrad caltet(const double, const double) const = 0;
    //! Calculate the four velocity of the corona fluid
    virtual void calumu(const double, const double, std::array<double, 4> &, sim5metric &) const = 0;
    //! Virtual destructor
    virtual ~tridgeo(){};
    //! Return the length scale of the corona
    virtual double length() const = 0;
};

//! ZAMO slab
class zamoslab3d : public tridgeo {
public:
    //! Corona radius
    double s;
    //! Minimum height
    double hmin;
    //! Maximum height
    double hmax;
    zamoslab3d(const double, const double, const double, const double);
    bool inside(const double r, const double mu) const;
    sim5tetrad caltet(const double, const double) const;
    void calumu(const double, const double, std::array<double, 4> &, sim5metric &) const;
    ~zamoslab3d(){};
    double length() const { return std::sqrt((hmax - hmin) * s);};
};

//! Slab that is co-rotating with the underlying Keplerian disc.
class kepslab3d : public tridgeo {
public:
    //! Corona radius
    double s;
    //! Minimum height
    double hmin;
    //! Maximum height
    double hmax;
    kepslab3d(const double, const double, const double, const double);
    bool inside(const double r, const double mu) const;
    sim5tetrad caltet(const double, const double) const;
    void calumu(const double, const double, std::array<double, 4> &, sim5metric &) const;
    double length() const { return hmax - hmin;};
    ~kepslab3d() {};
};
    
//! ZAMO sphere
class zamosphere3d : public tridgeo {
public:
    //! Height
    double h;
    //! Radius
    double R;
    zamosphere3d(const double, const double, const double);
    bool inside(const double r, const double mu) const;
    sim5tetrad caltet(const double, const double) const;
    void calumu(const double, const double, std::array<double, 4> &, sim5metric &) const;
    ~zamosphere3d() {};
    double length() const { return R;};
};

//! ZAMO spherical shell.
class zamosphere3d_truncated : public tridgeo {
public:
    //! Radius of outer boundary
    double R;
    //! Radius of inner boundary
    double rin;
    zamosphere3d_truncated(const double, const double, const double);
    bool inside(const double r, const double mu) const;
    sim5tetrad caltet(const double, const double) const;
    void calumu(const double, const double, std::array<double, 4> &, sim5metric &) const;
    ~zamosphere3d_truncated() {};
    double length() const { return R - rin;};
};

//! Spherical corona that is rotating with uniform linear circular velocity
class rotsphere3d: public tridgeo {
public:
    //! Height
    double h;
    //! Radius
    double R;
    //! Linear circular velocity
    double velocity;
    //! Lorentz factor, \f$\gamma=1/\sqrt{1 - v^2}\f$.
    double gamma;
    rotsphere3d(const double, const double, const double, const double);
    ~rotsphere3d() {};
    bool inside(const double r, const double mu) const;
    sim5tetrad caltet(const double, const double) const;
    void calumu(const double, const double, std::array<double, 4> &, sim5metric &) const;
    double length() const { return R;};
};

//! Abstract base class for corona in flat spacetime. No tetrad/\f$U^\mu\f$ calculation compared with \ref tridgeo
class tridgeo_sr {
public:
    //! Tell if the photon is inside. Accept std:array<double,3> argument, in accordance with \ref superphoton_sr.
    virtual bool inside(const std::array<double, 3> & xmu) const = 0;
    //! Return the length scale of the corona
    virtual double length() const = 0;
    //! Virtual constructor
    virtual ~tridgeo_sr(){};
};

//! Spherical corona in flat spacetime
class sphere3d_sr : public tridgeo_sr {
public:
    //! Radius
    double radius;
    //! Destructor
    ~sphere3d_sr(){};
    //! Constructor
    sphere3d_sr(const double rr) : radius(rr) {};
    //! Constructor without argument
    sphere3d_sr(){};
    bool inside(const std::array<double, 3> & xmu) const;
    double length() const { return radius;};
};

//! Slab corona in flat spacetime
class slab3d_sr : public tridgeo_sr {
public:
    //! Thickness
    double h;
    //! Radius
    double s;
    //! Destructor
    ~slab3d_sr(){};
    //! Constructor without argument
    slab3d_sr(const double hh, const double ss) : h(hh), s(ss) {};
    //! Constructor
    slab3d_sr(){};
    bool inside(const std::array<double, 3> & xmu) const;
    double length() const { return std::sqrt(s * h);};
};

//! co-rotating slab
class coroslab3d: public tridgeo {
public:
    //! Corona radius
    double s;
    //! Minimum height
    double hmin;
    //! Maximum height
    double hmax;
    coroslab3d(const double, const double, const double, const double);
    ~coroslab3d(){};
    
    bool inside(const double r, const double mu) const;
    sim5tetrad caltet(const double, const double) const;
    void calumu(const double, const double, std::array<double, 4> &, sim5metric &) const;
    double length() const { return std::sqrt((hmax - hmin) * s);};
};
    
#endif
