//! \file geoinf.h
#ifndef _GEOINF_H
#define _GEOINF_H
#include <iostream>
#include <array>
#include "tridgeo.h"
//! Class for photon geodesic in Kerr spacetime. Implemented various member functions for calculations related with photon motion in Kerr spacetime.
class geoinf{
public:
    //! BH spin
    double a;
    //! \f$a^2\f$
    double a2;
    //! radius of BH event horizon on the equatorial plane
    double rh;
    //! \f$\mu_\infty\f$, photon inclination at infinity
    double muobs;
    //! sign of \f$k^r\f$ at infinity
    double krobssign;
    //! sign of \f$k^\theta\f$ at infinity
    double kmuobssign;
    //! one impact parameter
    double alpha;
    //! another impact parameter
    double beta;
    //! \f$l\equiv L_z/E_\infty\f$, where \f$L_z\f$ is photon angular momentum about BH rotation axis and is one constant of motion, and \f$E_\infty\f$ is photon energy at infinity
    double l;
    //! \f$\mathcal{Q}/E^2_\infty\f$, where \f$\mathcal{Q}\f$ is another constant of motion (see Carter 1968)
    double q;
    //! A constant factor in evaluating \f$\mu-\f$integral; \f$q=1/\sqrt{l^2+q}\f$.
    double taumunorm;
    //! If \f$q\geq 0\f$: \f$=\int_0^{\mu_+} \mu^2d\mu/\sqrt{M(\mu)}\f$; otherwise 0.
    double dlambda_mu0;
    //! \f$r-\f$part of \f$\lambda\f$ integral between \f$r_2\f$ and \f$r_3\f$, the second and third largest roots of \f$R(r)=0\f$.
    double dlambda_r2r3;
    //! If \f$a\neq 0\f$: \f$\mu_+^2=(\sqrt{(l^2 + q - a^2)^2 + 4a^2q} -(l^2 + q - a^2)) / 2a^2\f$; otherwise \f$\mu_+^2=q/(l^2+q)\f$.
    double muplus2;
    //! \f$\mu_+\f$, the upper bound of \f$\mu\f$ allowed region
    double muplus;
    //! \f$\mu_-\f$. If \f$q<0\f$, \f$\mu_-\f$ is the lower bound of the allowed region above the equatorial plane (for the allowed region below the equatorial plane: \f$\mu\in[-\mu_+,-\mu_-]\f$).
    double muminus = 0.;
    //! If \f$a\neq 0\f$: \f$\mu_-^2=(\sqrt{(l^2 + q - a^2)^2 + 4a^2q} +(l^2 + q - a^2)) / 2a^2\f$; otherwise \f$\mu_-^2=0\f$.
    double muminus2;
    //! \f$m_\mu = \mu_+^2/(\mu_+^2 + \mu_-^2)\f$
    double mmu;
    //! If \f$q\geq 0\f$: \f$=\int_0^{\mu_+} d\mu/\sqrt{M(\mu)}\f$; otherwise \f$=\int_{\mu_-}^{\mu_+} d\mu/\sqrt{M(\mu)}\f$
    double taumu0;
    //! \f$\int_{\mu_\infty}^{\mu_+} d\mu/\sqrt{M(\mu)}\f$
    double taumuinf;
    //! obsolete; the \f$\mu-\f$ integral along the whole geodesic
    double taumu;
    //! \f$r-\f$integral between \f$r_2\f$ and \f$r_3\f$; \f$=\int_{r_2}^{r_3} dr/\sqrt{R(r)}\f$.
    double taur2r3;
    //! number of real roots of \f$R(r)=0\f$; see \ref quadroots
    int nrr;
    /*! \brief 
    If \f$R(r)=0\f$ has real roots, then *rr1* is \f$r_1\f$, the largest real root;
    otherwise *rr1* is \f$u_1\f$, the real part of the first pair of complex roots
    */
    double rr1;
    /*! \brief
    If \f$R(r)=0\f$ has real roots, then *rr2* is \f$r_2\f$, the second largest real root;
    otherwise *rr2* is \f$v_1\f$, the imaginary part of the first pair of complex roots
    */
    double rr2;
    /*! \brief If \f$R(r)=0\f$ has four real roots, then *rr3* is \f$r_3\f$, the third largest real root. 
     If \f$R(r)=0\f$ has two real roots, then *rr3* is
     \f$A=\sqrt{(r_1-u)^2 + v^2}\f$. If \f$R(r)=0\f$ has no real roots, then *rr3* is \f$u_2\f$, the real part of the second pair of complex roots.
     */
    double rr3;
    /*! \brief If \f$R(r)=0\f$ has four real roots, then *rr4* is \f$r_4\f$, the third largest real root. 
     If \f$R(r)=0\f$ has two real roots, then *rr4* is
     \f$B=\sqrt{(r_2-u)^2 + v^2}\f$. If \f$R(r)=0\f$ has no real roots, then *rr4* is \f$v_2\f$, the imaginary part of the second pair of complex roots.
     */
    double rr4; 
    //! If \f$R(r)=0\f$ has two real roots, then the complex roots are \f$u\pm iv\f$.
    double u;
    //! If \f$R(r)=0\f$ has two real roots, then the complex roots are \f$u\pm iv\f$.
    double v;
    //! The radius where the photon reaches the equatorial plane.
    double re;
    //! If two real roots: \f$m_4=((A+B)^2 - (r_3-r_4)^2)/(4AB)\f$. If four real roots: \f$m_4 = (r_2-r_3)(r_1-r_4) / (r_2-r_4) / (r_1 - r_3)\f$.
    double m4;
    //! A constant factor in \f$r-\f$ integral. If two real roots: \f$=1/\sqrt{AB}\f$. If four real roots: \f$=2/\sqrt{(r_1-r_3)(r_2-r_4)}\f$.
    double taurnorm;
    //! The \f$r-\f$ integral between \f$r_1\f$ and \f$\infty\f$.
    double taurinf;
    //! \f$=\sqrt{(r_1-r_3)(r_2-r_4)}\f$. Assigned only when four real roots.
    double temp;
    //! obsolete; number of turning points in \f$\mu\f$
    int nmuturn;
    //! obsolete; number of turning points in \f$r\f$
    int nrturn;
    //! an integer indicating the final destination of the photon as found out by \ref calfate
    int status;
    //double muturn = 0., rturn = 0.;
    //! if the photon enters the BH event horizon
    bool fallinbh = false;
    //! if \f$\mu_\infty\f$ is solved
    bool ismuobsset = false;
    //! if re is solved
    bool isresolved = false;
    //bool lowen = false;
    //! true if a=0
    bool spin0;

    geoinf(){};
    geoinf(const double aa, const double ll, const double qq);
    geoinf(const double aa, const double mu, const double alphaa, const double betaa);
    geoinf(const double aa, const double r, const double mu, std::array<double, 4> & kmu);
    //inline bool areturnsfound() const {return turnsfound;};
    bool ispossible(const tridgeo & trid, const double r, const double mu, const double kr, const double kmu);
    bool isallowed(const double r, const double mu) const;
    /*
    void correctkmu(const double re, const double mue, const double kre, const double kmue,
        const int nrr2, const int nrr3, const int nmumin, const int nmumax, std::array<double, 4> & pos, std::array<double, 4> & kmu);
    */
    //void findturns();
    void solvere();
    //void testturns() const;
    //double deltatau(double mu, const double & h) const;
    void caltaur_r2r3(); 
    double caltaur_r1(const double r) const;
    double caltaur_r2(const double r) const;
    //double caltaur_r3(const double r) const;
    double caltaur_inf_negq(const double r) const;
    double caltaur_negq(const double r1, const double r2) const;
    double caltaue(const double r) const;
    double calximu(const double mu) const;
    //double calximu_lowen(const double mu) const;
    double caltaur_all(const double re, const double rf, const double kre, const int nrr2, const int nrr3);
    double calximu_all(const double mue, const double muf, const double kmue, const int nmuminus, const int nmuplus);
    double caldlambda_muall(const double mue, const double muf, const double kmue,
        const int nmuminus, const int nmuplus, const double ximue, const double inte) const;
    double caldlambda_muplus(const double mu) const;
    //double caldlambda_muplus_lowen(const double mu) const;
    double caldlambda_r1(const double r) const;
    double caldlambda_r2(const double r) const;
    double caldlambda_negq(const double r1, const double r2) const;
    double caldlambda_rall(const double re, const double rf, const double kre, const int nrr2, const int nrr3, const double dlambdare, const double inte);
    double caldlambda_re(const double re) const;
    double solver_r1(const double tau) const;
    double solver_r2(const double tau) const;
    double solver_inf_negq(const double tau) const;
    double solver_inte(const double inte, const double re, const double kre, double & krfsign, int & nrr2, int & nrr3);
    double solver_inte(const double inte, const double re, const double kre, const double taue, double & krfsign, int & nrr2, int & nrr3);
    double solver_selfirr(const double re, const double kr, double & krfsign);
    double solvemu_inte(const double inte, const double mue, const double kmue, double & kmufsign, int & nmuminus, int & nmuplus);
    double solvemu_inte(const double inte, const double mue, const double kmue, const double xie, double & kmufsign, int & nmuminus, int & nmuplus);
    double solvemu(const double tau) const;
    void calmuobs(const double r, const double mu, const double kr, const double ktheta);
    //void calmuobs_disk(const double kr, const double re);
    int findxtrid(const tridgeo & trid, const double r, const double mu, const double kr,
        const double kmu, const double rin, double & rnow, double & munow, std::array<double, 4> & kmunow);
    int findxtrid_old(const tridgeo & trid, const double r, const double mu, const double kr,
    const double kmu, double & rnow, double & munow, std::array<double, 4> & kmunow);
    int findxtrid_sim5(const double rin, const tridgeo & trid, std::array<double, 4> & pos, std::array<double, 4> & kmu);
    bool locatere(const double r) const;
    double calfate(const double r, const double mu, const double kr, const double kmu, const double rin);
    double calfate_sim5(const double rin, std::array<double, 4> & pos, std::array<double, 4> & kmu);
private:
    void init();
    //bool turnsfound = false;
};

/*
class slabgeo {
public:
    geoinf * geo;
    double h;
    slabgeo(geoinf & geoo, double hh) : geo(&geoo), h(hh){};
};
std::ostream & operator << (std::ostream & out, const geoinf & geo);
double calslabrad(geoinf & geo, const double h);
double calslabrad_sim5(const geoinf & geo1, const double h, const double smax, double & rsign, double & tsign);
double deltatau(double mu, const slabgeo & slab);
*/
#endif
