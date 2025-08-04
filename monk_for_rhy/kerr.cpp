//! \file kerr.cpp
//! A few functions for various calculations in Kerr spacetime.
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <functional>
#include <random>
#include <complex>
#include "kerr.h"
#include "sim5lib.h"
#include "const.h"
#include "utils.h"

/*!
\brief Calculates the four velocity of particle doing Keplerian motion
@param a BH spin.
@param r radius in \f$[\rm GM~c^{-2}]\f$.
*/
std::vector<double> keplerumu(const double & a, const double & r) {
  std::vector<double> umu(4);
  umu[0] = (r * r + a * std::sqrt(r)) / r / std::sqrt(r * r - 3. * r + 2. * a * std::sqrt(r));
  umu[1] = 0.;
  umu[2] = 0.;
  umu[3] = 1. / std::sqrt(r) / std::sqrt(r * r - 3. * r + 2. * a * std::sqrt(r));
  return umu;
}

/*!
Reference: Dovciak+2004 (http://arxiv.org/abs/astro-ph/0407330), Equations A53-A56.
@param a BH spin
@param r radius in \f$[\rm GM~c^{-2}]\f$
@param l constant of motion
@param q constant of motion
@param rplus whether \f$k^r\f$ is positive
@param thetaplus whether \f$k^\theta\f$ is positive
*/
std::vector<double> fourp_disk(const double & a, const double & r, const double & l, const double & q,
  const bool rplus, const bool thetaplus) {
    std::vector<double> pmu(4);
    double delta = (r * r - 2. * r + a * a);
    pmu[0] = (a * (l - a) + (r * r + a * a) * (r * r + a * a - a * l) / delta) / r / r;
    double temp1 = r * r + a * a - a * l;
    temp1 *= temp1;
    pmu[1] = std::sqrt(temp1 - delta * ((l - a) * (l - a) + q)) / r / r;
    pmu[2] = std::sqrt(q) / r / r;
    pmu[3] = 2. * a / r / delta;
    pmu[1] *= (rplus ? 1. : -1.);
    pmu[2] *= (thetaplus ? 1 : -1);

    return pmu;
}

/*!
@param a BH spin
@param h height above the equatorial plane, in \f$[\rm GM~c^{-2}\f$
@param q constant of motion
@param rplus sign of \f$k^r\f$
@param thetaplus sign of \f$k^\theta\f$
*/
std::vector<double> fourp_axis(const double & a, const double &h, const double & q,
  const bool rplus, const bool thetaplus) {
    std::vector<double> pmu(4);
    double rho2 = h * h + a * a;
    double delta = h * h - 2. * h + a * a;
    double g00, g11, g22;
    g00 = -1. + 2. * h / rho2;
    g11 = rho2 / delta;
    g22 = rho2;
    pmu[0] = -1. / g00;
    pmu[2] = std::sqrt((q + a * a)) / g22;
    pmu[1] = std::sqrt(-1. * (1. / g00 + g22 * pmu[2] * pmu[2]) / g11);
    pmu[3] = 0.;
    pmu[1] *= (rplus ? 1. : -1.);
    pmu[2] *= (thetaplus ? 1 : -1);

    return pmu;
  }

/*
Return \f$\mu_i \equiv {\rm cos} \delta_i\f$ as measure in the disc fluid rest frame, where \f$\delta_i\f$ is the incident angle
of the photon arriving the disc on the equatorial plane, given constants of motion.
Reference: 
@param a bh spin
@param l constant of motion
@param q constant of motion
@param r radius where the photon arrives at the disc on the equatorial plane
double calmui(const double & a, const double & l, const double & q, const double & r) {
  double nom, denom;
  std::vector<double> umu, pmu;
  // four velocity of fluid
  umu = keplerumu(a, r);
  // four velocity of photon
  pmu = fourp_disk(a, r, l, q, true, true);

  sim5metric met;
  kerr_metric(a, r, 0., &met);
  // normal vector
  std::vector<double> nmu = {0., 0., 1. / std::sqrt(met.g22), 0.};

  nom = dotprod(pmu.data(), nmu.data(), &met);
  denom = dotprod(pmu.data(), umu.data(), &met);
  return -1. * nom / denom;
}
*/

/*!
Reference: Eq. 9 of Wilkins & Fabian 2012 (http://adsabs.harvard.edu/abs/2012MNRAS.424.1284W)
@param r radius in \f$[\rm GM~c^{-2}]\f$
@param p pointer to the parameter object. To be compatible with integration libraries,
    the params are better to be passed by a void pointer and then reinterpreted as 
    a \ref ds_params object within the function.
*/
double ds(const double r, void *p) {
  ds_params & params = *reinterpret_cast<ds_params *>(p);
  double a = params.a;
  double delta;
  delta = r * r - 2. * r + a * a;
  return 2. * M_PI * r * std::sqrt((r * r + a * a + 2. * a * a / r) / delta);
}

/*
//! calculate the proper area of a ring with inner and outer radii of
//! r0 and r1 in Kerr spacetime.
//! @param a bh spin
//! @param r0 inner radius
//! @param r1 outer radius
double calds_gsl(const double a, const double r0, const double r1) {
    double result, error;
    size_t neval;
    ds_params params;
    params.a = a;
    gsl_function F;
    F.function = &ds;
    F.params = reinterpret_cast<void *> (&params);

    const double epsabs=1e-4;
    const double epsrel=1e-4;

    int code = gsl_integration_qng(&F, r0, r1, epsabs, epsrel, &result, &error, &neval);
    if (code) {
        std::cerr << "Error in integration!" << std::endl;
        return 0.;
    } else {
        return result;
    }
}
*/

/*!
\f$\mathcal{S} = \int_{r_0}^{r_1} s(r)dr\f$, where \f$s(r)\f$ is defined in \ref ds. 
The integral is evaluated with the Gauss-Legendre method.
@param a BH spin
@param r0 lower boundary of the radial bin, in \f$[\rm GM~c^{-2}]\f$
@param r1 upper boundary of the radial bin, in \f$[\rm GM~c^{-2}]\f$
*/
double calds(const double a, const double r0, const double r1) {
    double weights[NGAUSS], abissca[NGAUSS];
    double inte = 0.;
    gauleg(r0, r1, abissca, weights, NGAUSS);
    
    ds_params params;
    params.a = a;
    void *p = reinterpret_cast<void *> (&params);
    
    for (long i = 0; i < NGAUSS; ++i)
        inte += ds(abissca[i], p) * weights[i];
    
    return inte;
}

/*!
@param tbb photon temperature with same unit as en
@param en input energy array
@param flux output array of specific intensity. Its unit is \f$\rm en^3\f$. To obtain
    the intensity in proper unit (e.g., \f$[\rm erg~s^{-1}~cm^{-2}~Hz^{-1}~sr^{-1}]\f$, one needs to 
    do further operation.
*/
void planck_kev(const double tbb, const std::vector<double> & en, std::vector<double> & flux) {
    double ra, ex;
    for (size_t i = 0; i < en.size(); ++i) {
        ra = en[i] / tbb;
        if (ra < 1e-3) {
            flux[i] = 2. * en[i] * en[i] * tbb;
        } else {
            ex = std::exp(-1. * en[i] / tbb);
            flux[i] = 2. * en[i] * en[i] * en[i] * ex / (1. - ex);
        }
    }
        //flux[i] = 2. * one_over_c2h3 * en[i] * en[i] * en[i] / (std::exp(en[i] / tbb) - 1.);
}

/*
Calculate polarization degree from pure scattering disc atomsphere with
infinite optical depth; following Eqs. 3-4 of Li 2009
@param mu cosine of polar emission angle
double semiplane_poldeg(const double mu) {
    return 0.;
}
*/

/*!
@param mu cosine of polar emission angle
*/
double phil(const double mu) {
    return 0.75 * (1. - mu * mu);
}

/*!
@param mu cosine of polar emission angle
*/
double phir(const double mu) {
    return 0.375 * (1. - mu * mu);
}

/*! 
One recursive form of \f$H(\mu)\f$ is
\f$1/H(\mu) = \left[1-2\int_0^1 \Phi(\mu^\prime)d\mu^\prime\right]^{1/2} +
\int_0^1[\mu^\prime\Phi(\mu^\prime)H(\mu^\prime)]/(\mu+\mu^\prime)d\mu^\prime\f$
(Eq. 13, Chandrasekhar 1960, http://adsabs.harvard.edu/abs/1960ratr.book.....C),
where \f$\Phi(\mu)\f$ is the characteristic function.
The evaluation follows a fast algorithm by
Bosma & de Rooij 1983 (http://adsabs.harvard.edu/abs/1983A%26A...126..283B).
@param mu array of cosine polar emission angle
@param res output; the values of H function
@param phifunc the characteristic function
*/
void hfunc(const std::vector<double> & mu, std::vector<double> & res, std::function<double (const double)> phifunc) {
    
    double phi0, iterconst, inte, deltamax, deltanow;
    
    //mu.push_back(0.);
    
    std::vector<double> weights(NGAUSS), abscissa(NGAUSS), absres(NGAUSS), absres1(NGAUSS), phimu(NGAUSS);
    std::vector<double> res1(mu.size());
    std::fill(res.begin(), res.end(), 1.);
    std::fill(absres.begin(), absres.end(), 1.5);
    
    gauleg(0., 1., abscissa.data(), weights.data(), NGAUSS);
    
    for (unsigned i = 0; i < NGAUSS; ++i)
        phimu[i] = phifunc(abscissa[i]);
    
    phi0 = 0.;
    for (unsigned long i = 0; i < NGAUSS; ++i)
        phi0 += weights[i] * phimu[i];
    iterconst = std::sqrt(1. - 2. * phi0);
    
    // find values at abscissa
    deltamax = 1.;
    for (unsigned long i = 0; i < MAXITER; ++i) {
        
        // Update abscissa values
        for (unsigned long j = 0; j < NGAUSS; ++j) {
            inte = 0.;
            for (unsigned long k = 0; k < NGAUSS; ++k)
                inte += (weights[k] * abscissa[k] * phimu[k] * absres[k] / (abscissa[j] + abscissa[k]));
            absres1[j] = 1. / (iterconst + inte);   // update
            deltanow = std::abs(absres1[j] - absres[j]);
            if (deltanow <= deltamax)
                deltamax = deltanow;
            absres[j] = absres1[j];
        }
        if (deltamax <= CRITERION) {
            break;
        }
    }
    
        // find value of desired mu arr
    for (unsigned long j = 0; j < mu.size(); ++j) {
        deltanow = 1.;
        for (unsigned long i = 0; i < MAXITER; ++i) {
            inte = 0.;
            for (unsigned long k = 0; k < NGAUSS; ++k)
                inte += (weights[k] * abscissa[k] * phimu[k] * absres[k] / (std::abs(mu[j]) + abscissa[k]));
            res1[j] = 1. / (iterconst + inte);   // update
            deltanow = std::abs(res1[j] - res[j]);
            res[j] = res1[j];
            if (deltanow <= CRITERION)
                break;
        }
    }
}

/*!
Calculate two Chandrasekhar H functions \f$H_l(\mu)\f$ and \f$H_r(\mu)\f$, corresponding to the characteristic
functions \ref phil() and \ref phir().
@param muarr the array of cosine polar emission angle
@param ilarr output; array of \f$H_l(\mu)\f$
@param irarr output; array of \f$H_r(\mu)\f$
*/
void evaluate_hfunc(const std::vector<double> & muarr, std::vector<double> & ilarr, std::vector<double> & irarr) {
    size_t nmu = muarr.size();
    
    std::vector<double> hl(nmu), hr(nmu);
    hfunc(muarr, hl, phil);
    hfunc(muarr, hr, phir);
    for (size_t i = 0; i < nmu; ++i) {
        ilarr[i] = hl[i] * POL_Q;
        irarr[i] = hr[i] * (std::abs(muarr[i]) + POL_C);
    }
}

/*!
@param muarr the array of cosine polar emission angle
@param poldeg output; array of pol. deg.
*/
void calpoldeg(const std::vector<double> & muarr, std::vector<double> & poldeg) {
    size_t nmu = muarr.size();
    std::vector<double> ilarr(nmu), irarr(nmu);
    evaluate_hfunc(muarr, ilarr, irarr);
    for (size_t i = 0; i < nmu; ++i)
        poldeg[i] = std::abs(irarr[i] - ilarr[i]) / (irarr[i] + ilarr[i]);
}

/*! 
Reference: Eq X.82, Chandrasekhar 1960.
@param muarr the array of cosine polar emission angle
@param dark output; the limb-darkening factor
*/
void calchandra_darkening(const std::vector<double> & muarr, std::vector<double> & dark) {
    double fac = 3. / 8. / std::sqrt(2.);
    size_t nmu = muarr.size();
    std::vector<double> ilarr(nmu), irarr(nmu);
    evaluate_hfunc(muarr, ilarr, irarr);
    for (size_t i = 0; i < nmu; ++i)
        dark[i] = (irarr[i] + ilarr[i]) * fac;
}

/*!
The Hamilton-Jacobi equation can be separated. After the separation, we can related two integrals of \f$r\f$ and
\f$\mu\equiv {\rm cos}\theta\f$:
\f{equation}
\int \frac{dr}{\sqrt{R(r)}} = \int \frac{d\mu}{\sqrt{M(\mu)}},
\f}
where
\f{eqnarray}
 R(r) &=& r^4 + (a^2 - l^2 - Q) r^2 + 2\left[Q + (l - a)^2\right] r - a^2 Q,\\
 M(\mu) &=& Q + (a^2 - l^2 - Q)\mu^2  - a^2 \mu^4.
\f}
This function calculates the \f$r-\f$integrand, i.e., \f$1/\sqrt{R(r)}\f$.
@param r radius in \f$[\rm GM~c^{-2}]\f$
@param params parameter array, \f$[a, l, q]\f$
*/
double kerr_rintegrand(const double r, const std::array<double, 3> & params) {
    double a = params[0], l = params[1], q = params[2];
    double a2 = a * a;
    double R, r2 = r * r;
    R = r2 * r2 + (a2 - l * l - q) * r2 + 2. * (q + (l - a) * (l - a)) * r - a2 * q;
    
    return 1. / std::sqrt(R);
}

/*! \brief Evaluates \f$r/\sqrt{R(r)}\f$, the integrand of \f$dt_{r, p0}\f$; see \ref geoinf_tphi::caldt_r1_p2().
@param r radius
@param params parameter array, \f$[a, l, q]\f$
 */
double kerr_dtintegrand_r_p0(const double r, const std::array<double, 3> & params) {
    return r * kerr_rintegrand(r, params);
}

/*!
The \f$\mu-\f$ integrand is \f$1/\sqrt{M(\mu)}\f$ (see Eq. 1 in \ref kerr_rintegrand()).
@param mu \f$\equiv {\rm cos}\theta\f$, where \f$\theta\f$ is one of the Boyer-Lindquist coordinates
@param params parameter array, \f$[a, l, q]\f$
 */
double kerr_muintegrand(const double mu, const std::array<double, 3> & params) { 
    double a = params[0], l = params[1], q = params[2];
    double T, mu2 = mu * mu;
    T = q + (a * a - l * l - q) * mu2 - a * a * mu2 * mu2;
    return std::sqrt(1. / T);
}
    
/*!
The affine parameter: 
\f{equation}
\label{eq:lambda}
 \lambda = \int \frac{r^2dr}{\sqrt{R(r)}} + a^2 \int \frac{\mu^2 d\mu}{\sqrt{M(\mu)}}.
\f}
This function calcuates the \f$\mu-\f$integrand, i.e., \f$\frac{\mu^2 d\mu}{\sqrt{M(\mu)}}\f$.
@param mu \f$\equiv {\rm cos}\theta\f$, where \f$\theta\f$ is one of the Boyer-Lindquist coordinates
@param params parameter array, \f$[a, l, q]\f$
 */
double kerr_lambdaintegrand_mu(const double mu, const std::array<double, 3> & params) {
    double a = params[0], l = params[1], q = params[2];
    double T, mu2 = mu * mu;
    T = q + (a * a - l * l - q) * mu2 - a * a * mu2 * mu2;
    return a * a * mu2 * std::sqrt(1. / T);
}

/*!
This function calcuates the \f$r-\f$integrand in the equation in \ref kerr_lambdaintegrand_mu, i.e., \f$\frac{r^2dr}{\sqrt{R(r)}}\f$.
@param r radius in \f$[\rm GM~c^{-2}]\f$
@param params parameter array, \f$[a, l, q]\f$
*/
double kerr_lambdaintegrand_r(const double r, const std::array<double, 3> & params) {
    double a = params[0], l = params[1], q = params[2];
    double R, r2 = r * r;
    R = r2 * r2 + (a * a - l * l - q) * r2 + 2. * (q + (l - a) * (l - a)) * r - a * a * q;
    
    return r2 / std::sqrt(R);
}

/*! \brief Evaluates \f$l/(1-\mu^2)\sqrt{M(\mu)}\f$, the integrand of \f$d\phi_\mu\f$; see \ref geoinf_tphi::caldphimu.
@param mu \f$\equiv {\rm cos}\theta\f$, where \f$\theta\f$ is one of the Boyer-Lindquist coordinates
@param params parameter array, \f$[a, l, q]\f$
 */
double kerr_phiintegrand_mu(const double mu, const std::array<double, 3> & params) {
    double a = params[0], l = params[1], q = params[2];
    double T, mu2 = mu * mu;
    T = q + (a * a - l * l - q) * mu2 - a * a * mu2 * mu2;
    return l / (1. - mu2) / std::sqrt(T);
}

/*! \brief Evaluates \f$a(2r-al)/\Delta \sqrt{R(r)}\f$, the integrand of \f$d\phi_r\f$; see \ref geoinf_tphi::caldphimu
@param r radius
@param params parameter array, \f$[a, l, q]\f$
*/
double kerr_phiintegrand_r(const double r, const std::array<double, 3> & params) {
    double a = params[0], l = params[1], q = params[2];
    double R, r2 = r * r, Delta;
    R = r2 * r2 + (a * a - l * l - q) * r2 + 2. * (q + (l - a) * (l - a)) * r - a * a * q;
    Delta = r2 + a * a - 2. * r;
    return a * (2. * r - a * l) / Delta / std::sqrt(R);
}

/*!
Reference: Page & Thorne 1974, http://adsabs.harvard.edu/abs/1974ApJ...191..499P.
@param r radius in \f$[\rm GM~c^{-2}]\f$
@param a BH spin
 */
double nt_ffunc(const double r, const double a) {
    double x, xm1, xm2, xm3, xm4, a2;
    double varc, varcderiv, xvarf, xvarfderiv, item1, dLdx, dLdr;
    
    x = std::sqrt(r);
    xm1 = 1. / x;
    xm2 = xm1 * xm1;
    xm3 = xm1 * xm2;
    xm4 = xm3 * xm1;
    a2 = a * a;
    
    varc = 1. - 3. * xm2 + 2. * a * xm3;
    xvarf = x - 2. * a * xm2 + a2 * xm3;
    xvarfderiv = 1. + 4. * a * xm3 - 3. * a2 * xm4;
    varcderiv = 6. * (xm3 - a * xm4);
    
    item1 = std::sqrt(varc) / (1. + a * xm3);
    
    dLdx = (xvarfderiv * varc - 0.5 * xvarf * varcderiv) / varc / std::sqrt(varc);
    dLdr = 0.5 * dLdx * xm1;
    
    return dLdr * item1;
}

/*!
@param a BH spin
@param re emission radius
@param theta photon polar emission angle measured by disc fluid
@param phi photon azimuthal emission angle measured by disc fluid
@param bl_k output; photon wave vector in Boyer-Lindquist frame
 */ 
void kmu_disk_angle(const double a, const double re, const double theta, const double phi, std::array<double, 4> & bl_k) {
    double x, omega, sint, knorm;
    sim5metric met;
    sim5tetrad tet;
    std::array<double, 4> on_k;
    
    x = std::sqrt(re);
    omega = 1. / (x * x * x + a);
    
    kerr_metric(a, re, 0., &met);
    tetrad_azimuthal(&met, omega, &tet);
        
    sint = std::sin(theta);
    on_k[0] = 1.;
    on_k[2] = std::cos(theta);
    on_k[1] = sint * std::cos(phi);
    on_k[3] = sint * std::sin(phi);
    
    on2bl(on_k.data(), bl_k.data(), &tet);
                // normalize wave vector such that k_t = -1;
    knorm = -1. / (met.g00 * bl_k[0] + met.g03 * bl_k[3]);
    std::transform(bl_k.cbegin(), bl_k.cend(), bl_k.begin(), std::bind2nd(std::multiplies<double>(), knorm));
}

/*
void zamokmu_angle(const double a, const double re, const double mu, const double theta, const double phi, std::array<double, 4> & bl_k) {
    double sint, knorm;
    sim5metric met;
    sim5tetrad tet;
    std::array<double, 4> on_k;
    
    sint = std::sin(theta);
    kerr_metric(a, re, mu, &met);
    
    if (mu < 1.) {
        tetrad_zamo(&met, &tet);
        
        on_k[0] = 1.;
        on_k[2] = std::cos(theta);
        on_k[1] = sint * std::cos(phi);
        on_k[3] = sint * std::sin(phi);
        
        on2bl(on_k.data(), bl_k.data(), &tet);
                    // normalize wave vector such that k_t = -1;
        knorm = -1. / (met.g00 * bl_k[0] + met.g03 * bl_k[3]);
        std::transform(bl_k.cbegin(), bl_k.cend(), bl_k.begin(), std::bind2nd(std::multiplies<double>(), knorm));
    } else {
        knorm = 1. / std::sqrt(-1. * met.g00);
        bl_k[0] = knorm / std::sqrt(-1. * met.g00);
        bl_k[1] = knorm * std::cos(theta) / std::sqrt(met.g11);
        bl_k[2] = knorm * sint / std::sqrt(met.g22);
        bl_k[3] = 0.;
        //std::cout << "k * k = " << dotprod(bl_k.data(), bl_k.data(), &met) << std::endl;
    }
}
*/

/*!
\brief  A wrapper for `photon_momentum()` in \f$\textsc{sim5}\f$.

This function also offers special treatment for photon
on the rotation axis
@param a BH spin
@param r radius in \f$[\rm GM~c^{-2}]\f$.
@param mu \f$\mu\equiv {\rm cos}\theta\f$
@param l constant of motion
@param q constant of motion
@param krsign sign of \f$k^r\f$
@param kthetasign sign of \f$k^\theta\f$
@param kmu output; photon wave vector
*/
void photon_momentum1(const double a, const double r, const double mu, const double l, const double q, 
    const double krsign, const double kthetasign, std::array<double, 4> & kmu) {
    
    if (mu < 1.) {
        photon_momentum(a, r, mu, l, q, krsign, kthetasign, kmu.data());
    } else {
        sim5metric met;
        //std::cout << "inside photon_momentum1()" << std::endl;
        kerr_metric(a, r, mu, &met);
        kmu[0] = -1. / met.g00;
        kmu[2] = kthetasign * std::sqrt(q + a * a) / met.g22;
        kmu[1] = krsign * std::sqrt(-1. * (met.g00 * kmu[0] * kmu[0] + met.g22 * kmu[2] * kmu[2]) / met.g11);
        kmu[3] = 0.;
    }
}
        
/*!
Following Pozdnyakov+1983 (http://adsabs.harvard.edu/abs/1983ASPRv...2..189P).
@param tbb photon temperature.
@param gen 64-bit MT random number generator
 */
double sample_bb(const double tbb, std::mt19937_64 & gen) {
    std::uniform_real_distribution<> dis(0., 1.);
    double xi1, xi2, xi3, xi4, alpha, m = 2., en;
    double series0, series1;
    
    series0 = 1.;
    series1 = 9. / 8.;
    
    xi1 = dis(gen);
    xi2 = dis(gen);
    xi3 = dis(gen);
    xi4 = dis(gen);
    
    if (xi1 < 1. / APERY_CONST) {
        alpha = 1.;
    } else {
        while(true) {
            if (xi1 >= series0 / APERY_CONST && xi1 < series1 / APERY_CONST) {
                alpha = m;
                break;
            } else {
                ++m;
                series0 = series1;
                series1 += 1. / m / m / m;
            }
        }
    }
    
    en = -1. * tbb * std::log(xi2 * xi3 * xi4) / alpha;
    
    return en;
}

/*!
Following Pozdnyakov+1983 (http://adsabs.harvard.edu/abs/1983ASPRv...2..189P).
To increase the statistics at high energy, we increase the probability of photons above `ecut` and decrease their weight correspondingly.
We use a highly inefficient rejection method (1/ratio of \ref sample_bb).
@param tbb photon temperature.
@param gen 64-bit MT random number generator
@param ecut the critical energy above which the probability will be increased by a factor of `ratio`, and correspondingly weight
    dereased by the same factor
@param ratio the factor by which the probability increased
@param en energy of the sampled photon
@param weight statistical weight of the sampled photon
*/
void sample_bb1(const double tbb, std::mt19937_64 & gen, const double ecut, const double ratio, double & en, double & weight) {
    std::uniform_real_distribution<> dis(0, 1);
    bool found = false;
    double r1;
    while (!found) {
        en = sample_bb(tbb, gen);
        if (en >= ecut) {
            found = true;
            weight = 1. / ratio;
        } else {
            r1 = dis(gen);
            if (r1 < 1. / ratio) {
                found = true;
                weight = 1.;
            }
        }
    }
}

/*!
The probability density distribution
\f$p(x) \propto x^\gamma\f$ between \f$x_0\f$ and \f$x_1\f$.
@param gamma powerlaw index \f$\gamma\f$
@param x0 lower limit
@param x1 upper limit
@param gen 64-bit MT19937 random number generator
*/
double sample_pl(const double gamma, const double x0, const double x1, std::mt19937_64 & gen) {
    double x, y, temp;
    std::uniform_real_distribution<> dis(0, 1);
    y = dis(gen);
    
    temp = (std::pow(x1, gamma + 1.) - std::pow(x0, gamma + 1.)) * y + std::pow(x0, gamma + 1.);
    x = std::pow(temp, 1. / (gamma + 1.));
    return x;
}
    
/*!
Reference: Li+2009 (http://adsabs.harvard.edu/abs/2009ApJ...691..847L).
@param muobs cosine of observer inclination
@param ktheta sign of \f$k_\theta\f$ at infinity
@param a BH spin
@param l motion constant
@param q motion constant
@param poldeg polarization degree
@param wp Walker-Penrose constant
@param qfac output; Q/I, Q and I are two of the four Stokers parameters
@param ufac output; U/I, U and I are two of the four Stokers parameters
*/
void cal_pol_inf(const double muobs, const double ktheta, const double a, const double l, const double q, 
    const double poldeg, const std::complex<double> wp, double & qfac, double & ufac) {
    // calculate two motion constants S and T; LNM09, Eqs. B12-13
    double signktheta, sintheta, S, T, s2plust2, ftheta, fphi, sin2phi, cos2phi;
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
    ftheta = -1. * (wp.imag() * S + wp.real() * T) / s2plust2;
    fphi = (wp.imag() * T - wp.real() * S) / s2plust2;
    
    sin2phi = 2. * ftheta * fphi;
    cos2phi = 2. * ftheta * ftheta - 1.;
    qfac = poldeg * cos2phi;
    ufac = poldeg * sin2phi;
}

/*!
 */
double sample_cutoffpl(const double alpha, const double emin, const double ecut, std::mt19937_64 & gen) {
    double emax = 10. * ecut;
    std::uniform_real_distribution<> dis(0., 1.);
    double x, xi;
    bool find = false;
    while (!find) {
        x = sample_pl(alpha, emin, emax, gen);
        xi = dis(gen);
        if (xi < std::exp(-1. * x / ecut))
            find = true;
    }
    return x;
}

double sample_expo(const double lambda, std::mt19937_64 & gen) {
    std::uniform_real_distribution<> dis(0., 1.);
    double xi = dis(gen);
    return -1. * lambda * std::log(1. - xi);
}
    
