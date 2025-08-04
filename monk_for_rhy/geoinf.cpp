//! \file geoinf.cpp
#include <cmath>
#include <array>
#include <iostream>
#include <algorithm>
#include "geoinf.h"
#include "quadroots.h"
#include "sim5lib.h"
#include "rootsearch.h"
#include "utils.h"
#include "tridgeo.h"
#include "const.h"

/*! \brief Constructor with impact factors
@param aa BH spin
@param mu cosine of inclination
@param alphaa impact factor
@param betaa impact factor
*/
geoinf::geoinf(const double aa, const double mu, const double alphaa, const double betaa) {
    a = aa;
    a2 = a * a;

    spin0 = (std::abs(a) <= 1e-6);

    muobs = mu;
    alpha = alphaa;
    beta = betaa;
    l = -1. * alpha * std::sqrt(1. - muobs * muobs);
    q = beta * beta + (alpha * alpha - a2) * muobs * muobs;

    if (alpha == 0.) {
        muplus2 = 1.;
    } else if (beta == 0.) {
        muplus2 = muobs * muobs;
    } else {
        muplus2 = spin0 ? q / (l * l + q) : (std::sqrt((l * l + q - a2) * (l * l + q - a2) + 4. * a2 * q) - (l * l + q - a2)) / 2. / a2;
    }

    if (beta == 0) {
        muplus = muobs;
    } else {
        muplus = std::sqrt(muplus2);
    }

    ismuobsset=true;
    init();
}

/*! \brief Constructor; with position and wave vector
@param aa BH spin
@param r current Boyer-Lindquist coordinate r
@param mu cosine of current Boyer-Lindquist coordinate \f$\theta\f$
@param kmu wave vector in Boyer-Lindquist frame
*/
geoinf::geoinf(const double aa, const double r, const double mu, std::array<double, 4> & kmu) {
    a = aa;
    a2 = a * a;
    if (mu < 1.) {
        photon_motion_constants(a, r, mu, kmu.data(), &l, &q);
    } else {
        l = 0.;
        double g22 = r * r + a2;
        double kon_theta = g22 * kmu[2];
        q = kon_theta * kon_theta - a2;
        //std::cout << "kon_theta = " << kon_theta << std::endl;
        //std::cout << "q = " << q << std::endl;
    }

    spin0 = (std::abs(a) <= 1e-6);

    muplus2 = spin0 ? q / (l * l + q) : (std::sqrt((l * l + q - a2) * (l * l + q - a2) + 4. * a2 * q) - (l * l + q - a2)) / 2. / a2;
    muplus = std::sqrt(muplus2);

    ismuobsset=false;
    init();
}

/*! \brief Constructor; from constants of motion
@param aa BH spin
@param ll \f$l\equiv L_z/E_\infty\f$, where \f$L_z\f$ is phonton angular momentum about spin axis, and \f$E_\infty\f$ is photon energy at infinity
@param qq \f$Q\equiv \mathcal{Q}/E_\infty^2\f$, where \f$\mathcal{Q}\f$ is another constant of motion
*/
geoinf::geoinf(const double aa, const double ll, const double qq) {
    a = aa;
    a2 = a * a;

    spin0 = (std::abs(a) <= 1e-6);
    // the case for E_infty \approx 0

    l = ll;
    q = qq;

    muplus2 = spin0 ? q / (l * l + q) : (std::sqrt((l * l + q - a2) * (l * l + q - a2) + 4. * a2 * q) - (l * l + q - a2)) / 2. / a2;
    muplus = std::sqrt(muplus2);
    init();
}

/*! \brief Pre-calculate a few values 
 */
void geoinf::init() {
    double temp3;
    //lowen = (l >= 1e4 && q >= 1e4);

    rh = 1. + std::sqrt(1. - a2);

    muminus2 = std::abs(q) / a2 / muplus2;
    if (spin0) {
        taumunorm = 1. / std::sqrt(l * l + q);
        taumu0 = taumunorm * 0.5 * M_PI;
        dlambda_mu0 = 0.;
    } else {
        if (q >= 0.) {
            mmu = muplus2 / (muplus2 + muminus2);
            taumunorm = 1. / std::abs(a) / std::sqrt(muplus2 + muminus2);
        } else {
            muminus = std::sqrt(muminus2);
            mmu = 1. - muminus2 / muplus2;
            taumunorm = 1. / std::abs(a) / muplus;
        }
        taumu0 = elliptic_k(mmu) * taumunorm;
        dlambda_mu0 = elliptic_e_sin(1., mmu) / taumunorm;
    }


    if (ismuobsset) {
        if (beta == 0.) {
            taumuinf = 0.;
        } else {
            taumuinf = spin0 ? calximu(muobs) : jacobi_icn(muobs / muplus, mmu) * taumunorm;
        }

        // if beta < 0, no turning point in mu
        if (beta < 0.) {
            nmuturn = 0;
            taumu = taumu0 - taumuinf;
        } else {
            nmuturn = 1;
            taumu = taumu0 + taumuinf;
        }
    }

    // then comes to radius
    quadroots roots(a, l, q);
    roots.sort();
    nrr = roots.nrr;

    if (nrr == 0) {
        nrturn = 0;
        rr1 = roots.roots[0].real();
        rr2 = roots.roots[0].imag();
        rr3 = roots.roots[2].real();
        rr4 = roots.roots[2].imag();
        u = rr3;
        v = rr4;

    } else {
        rr1 = roots.roots[0].real();
        rr2 = roots.roots[1].real();

        if (nrr == 2) {
            u  = roots.roots[2].real();
            v  = roots.roots[2].imag();
            if (q >= 0.) {
                rr3 = sqrt(sqr(rr1-u) + v * v);
                rr4 = sqrt(sqr(rr2-u) + v * v);
            } else {
                rr3 = u;
                rr4 = v;
            }

            if (q >= 0.) {
                m4 = (sqr(rr3+rr4) - sqr(rr1-rr2)) / (4.*rr3*rr4);
                taurnorm = 1. / std::sqrt(rr3 * rr4);
                taurinf = taurnorm * jacobi_icn((rr3 - rr4) / (rr3 + rr4), m4);
            }
        } else if (nrr == 4) {
            rr3 = roots.roots[2].real();
            rr4 = roots.roots[3].real();
            temp = sqrt((rr1 - rr3) * (rr2 - rr4));
            taurnorm = 2. / temp;
            temp3 = (rr2 - rr4) / (rr1 - rr4);
            m4 = (rr2 - rr3) / (rr1 - rr3) / temp3;
            taurinf = taurnorm * jacobi_isn(std::sqrt(temp3), m4);
            if (rr2 > 0.) {
                taur2r3 = taurnorm * elliptic_k(m4);
                dlambda_r2r3 = caldlambda_r2(rr3);
            }
        }

        nrturn = ((taurinf > taumu) ? 0 : 1);
    }
}

/*!
\brief Obsolete; solve the radius where the photon reach the disc on the equatorial plane. Now this work is done by \ref calfate.
 */
void geoinf::solvere() {
    if (rr1 <= rh && nrturn == 1) {
        fallinbh = true;
    } else {
        re = solver_r1(std::abs(taumu - taurinf));
        if (re <= rh) {
            re = 0.;
            fallinbh = true;
        }
    }
}

/*!
\brief Invert R-integral \f$\int_{r_1}^r \frac{1}{\sqrt{R(r)}} =\tau_r\f$ to solve r, where \f$r_1\f$ is the largest root
of \f$R(r)=0\f$. For definition of \f$R(r)\f$ see \ref kerr_rintegrand.

This works when \f$Q \geq 0\f$ and \f$r\geq r_1\f$;

Reference: Li+2005 (http://adsabs.harvard.edu/abs/2005ApJS..157..335L)
@param tau the integral \f$\tau_r\f$.
*/
double geoinf::solver_r1(const double tau) const {
    if (q < 0.) {
        throw("geoinf::solver_r1: q < 0");
    }

    if (tau > taurinf) {
        throw("geoinf::solver_r1: intrtoolarge");
    }

    double xi2, cn, xi4, sn2, res;
    if (nrr == 2) {
        xi2 = tau * std::sqrt(rr3 * rr4);
        cn = jacobi_cn(xi2, m4);
        res = (rr2 * rr3 - rr1 * rr4 - (rr2 * rr3 + rr1 * rr4) * cn) /
            (rr3 - rr4 - (rr3 + rr4) * cn);
    } else {
        xi4 = tau * temp / 2.;
        sn2 = jacobi_sn(xi4, m4);
        sn2 *= sn2;
        res = (rr1 * (rr2 - rr4) - rr2 * (rr1 - rr4) * sn2) /
            (rr2 - rr4 - (rr1 - rr4) * sn2);
    }
    return res;
}

/*!
\brief Invert R-integral \f$\int_{r}^{r_2} \frac{1}{R(r)}=\tau_r\f$ to solve r. Here \f$r_2\f$ is the second largest root of \f$R(r)=0\f$.

This works only for the case when R(r)=0 has four real roots, and \f$r\leq r_2\f$.
Reference: Li+2005.
@param tau the integral \f$\tau_r\f$.
*/
double geoinf::solver_r2(const double tau) const {
    double sn2, xi4, res;
    if (q < 0.) {
        throw("geoinf::solver_r2: q < 0");
    }
    if (nrr != 4) {
        throw("geoinf::solver_r2: nrr != 4");
    }

    if (tau > taur2r3) {
        throw("geoinf::solver_r2: intrtoolarge");
    }

    xi4 = tau * temp / 2.;
    sn2 = jacobi_sn(xi4, m4);
    sn2 *= sn2;

    res = (rr2 * (rr1 - rr3) - rr1 * (rr2 - rr3) * sn2) /
        (rr1 - rr3 - (rr2 - rr3) * sn2);

    return res;
}

/*!
\brief Invert R-integral \f$\int_{r}^\infty \frac{1}{R(r)}=\tau_r\f$ to solve r. 

This function works for the case of \f$q < 0\f$.
In the case of 0 real root, following Dexter & Agol 2009 (http://adsabs.harvard.edu/abs/2009ApJ...696.1616D), Table 2; otherwise following Li+2005 Eq. A14.
@param tau the integral \f$\tau_r\f$.
*/
double geoinf::solver_inf_negq(const double tau) const {
    if (q >= 0.)
        throw("geoinf::solver_inf_negq: qge0");

    double e, m, n, c1, m1, c4, c5, temp3, jval, uf;
    temp3 = rr3 * rr3 + rr4 * rr4;
    e = -1. * a2 * q;

    // Two pairs of conjugate complex roots
    // Case 6 of Dexter & Agol 2009 table 2
    if (nrr == 0) {
        double c2, c3, temp1, scj, ic3;
        double p, r;
        double temp2 = rr1 * rr1 + rr2 * rr2;
        m = rr1 / temp2;
        n = std::abs(rr2 / temp2);
        p = rr3 / temp3;
        r = std::abs(rr4 / temp3);

        if (m < p) {
            temp2 = m;
            m = p;
            p = temp2;
            temp3 = n;
            n = r;
            r = temp3;
        }

        temp1 = (m - p) * (m - p);
        c4 = std::sqrt(temp1 + (n + r) * (n + r));
        c5 = std::sqrt(temp1 + (n - r) * (n - r));
        c1 = std::sqrt(e) * (c4 + c5) / 2.;
        c2 = std::sqrt((4. * n * n - (c4 - c5) * (c4 - c5)) / ((c4 + c5) * (c4 + c5) - 4. * n * n));
        c3 = m + c2 * n;
        ic3 = caltaur_inf_negq(1. / c3);
        m1 = 1. - (c4 - c5) * (c4 - c5) / (c4 + c5) / (c4 + c5);
        jval = tau - ic3;
        scj = jacobi_sn(c1 * jval, m1) / jacobi_cn(c1 * jval, m1);

        uf = c3 + (n * (1. + c2 * c2) * scj) / (1. - c2 * scj);
    } else {
        // two real roots and one pair of complex conjugate roots
        // Case 5 of Dexter & Agol 2009 table 2
        double cnj, u1, u4, inteb;

        u1 = 1. / rr1;
        u4 = 1. / rr2;

        //std::cout << "u1 = " << u1 << ", u4 = " << u4 << std::endl;
        m = rr3 / temp3;
        n = std::abs(rr4) / temp3;
        c4 = std::sqrt((m - u4) * (m - u4) + n * n);
        c5 = std::sqrt((m - u1) * (m - u1) + n * n);

        m1 = ((c4 + c5) * (c4 + c5) - (u4 - u1) * (u4 - u1)) / 4. / c4 / c5;
        c1 = std::sqrt(e * c4 * c5);

        inteb = jacobi_icn((u4 * c5 - u1 * c4) / (-1. * u4 * c5 - u1 * c4), m1) / c1;

        jval = c1 * (tau + inteb);
        cnj = jacobi_cn(jval, m1);
        uf = (u4 * c5 - u1 * c4 + (u4 * c5 + u1 * c4) * cnj) / ((c4 + c5) * cnj - c4 + c5);
    }

    return 1. / uf;
}

/*!
\brief Calculate \f$\int_{r_1}^r \frac{1}{R(r)}\f$. The inverse function of \ref solver_r1.

Reference see Li+2005, Eqs. A12 and A14.
@param r radius
*/
double geoinf::caltaur_r1(const double r) const {
    double tau;
    if (std::abs(r - rr1) <= DEPS)
        return 0.;
    if (r < rr1)
        throw("geoinf::caltaur_r1(): r < r1");
    if (q < 0.)
        throw("geoinf::caltaur_r1(): q < 0");
    if (nrr == 2) {
        tau = taurnorm * jacobi_icn(((rr3 - rr4) * r + rr1 * rr4 - rr2 * rr3) /
            ((rr3 + rr4) * r - rr1 * rr4 - rr2 * rr3), m4);
    } else {
        tau = taurnorm * jacobi_isn(std::sqrt((rr2 - rr4) * (r - rr1) / (rr1 - rr4) / (r - rr2)), m4);
    }
    return tau;
}

/*!
\brief Calculate \f$\int_{r_3}^{r_2}\frac{dr}{\sqrt{R(r)}}\f$,
where \f$r_2\f$ and \f$r_3\f$ are the second and third largest roots of \f$R(r)=0\f$, respectively.

This works only for the case when R(r)=0 has four real roots.
*/
void geoinf::caltaur_r2r3() {
    if (nrr != 4)
        throw("geoinf::caltaur_r2r3(): nrr != 4");
    if (q < 0.)
        throw("q < 0");

    taur2r3 = taurnorm * elliptic_k(m4);
}

/*!
\brief Calculate \f$\int_{r}^{r_2}\frac{dr}{\sqrt{R(r)}}\f$. The inverse function of \ref solver_r2.

Reference: Li+2005.
@param r radius
*/
double geoinf::caltaur_r2(const double r) const {
    if (std::abs(r - rr2) <= DEPS)
        return 0.;
    if (std::abs(r - rr3) <= DEPS) {
        return taur2r3;
    }
    if (nrr != 4)
        throw("geoinf::caltaur_r2(): nrr != 2");
    if (r > rr2) {
        std::cout << "nrr = " << nrr << " rr1 = " << rr1 << ", rr2 = " << rr2 << std::endl;
        std::cout << "r = " << r << std::endl;
        throw("geoinf::caltaur_r2(): r > r2");
    }
    if (r < rr3)
        throw("geoinf::caltaur_r2(): r < r3");

    return taurnorm * jacobi_isn(std::sqrt((rr1 - rr3) * (rr2 - r) / (rr2 - rr3) / (rr1 - r)), m4);
}

/*!
\brief Evaluate the integral \f$\int_r^{\infty} \frac{dr}{\sqrt{R(r)}}\f$, in the case of Carter's constant \f$q < 0\f$. Inverse function of
\ref solver_inf_negq.

In the case of 0 real root, the integral is evaluated following Carlson 1991 eq. 2.14, otherwise following Carlson 1992 eq. 2.36, .
//! @param r radius
*/
double geoinf::caltaur_inf_negq(const double r) const {
    double f2, g2;
    double eta2, M2;
    double Delta, Delta_plus, Delta_minus, L2_plus, L2_minus, delta12sqr, delta11sqr, delta22sqr;
    double r2 = r * r;

    if (q >= 0.) {
        throw("geoinf::caltaur_inf_negq: q>=0");
    }

    f2 = rr3 * rr3 + rr4 * rr4;
    g2 = -2. * rr3;
    eta2 = std::sqrt(f2 + g2 * r + r2);

    if (nrr == 0) {
        double f1, g1;
        double eta1;

        f1 = rr1 * rr1 + rr2 * rr2;
        g1 = -2. * rr1;

        eta1 = std::sqrt(f1 + g1 * r + r2);
        M2 = (2. * eta1 + g1 + 2. * r) * (2. * eta2 + g2 + 2. * r);

        delta12sqr = 2. * (f1 + f2) - g1 * g2;
        delta11sqr = 4. * f1 - g1 * g1;
        delta22sqr = 4. * f2 - g2 * g2;
        Delta = std::sqrt(delta12sqr * delta12sqr - delta11sqr * delta22sqr);
        Delta_plus = delta12sqr + Delta;
        Delta_minus = delta12sqr - Delta;
        L2_plus = M2 + Delta_plus;
        L2_minus = M2 + Delta_minus;
    } else {
        double Y1, Y4;
        double c14sqr, c11, c44;
        c14sqr = 2. * f2 + g2 * (rr1 + rr2) + 2. * rr1 * rr2;
        c11 = std::sqrt(2. * (f2 + g2 * rr1 + rr1 * rr1));
        c44 = std::sqrt(2. * (f2 + g2 * rr2 + rr2 * rr2));

        Y1 = std::sqrt(r - rr1);
        Y4 = std::sqrt(r - rr2);
        M2 = (Y1 + Y4) * (Y1 + Y4) * (2. * eta2 + g2 + 2. * r);

        L2_plus = M2 + c14sqr + c11 * c44;
        L2_minus = M2 + c14sqr - c11 * c44;

    }

    return 4. * rf(M2, L2_minus, L2_plus);
}

/*!
\brief Evaluate the integral \f$\int_{r_a}^{r_b} \frac{dr}{\sqrt{R(r)}}\f$, in the case of Carter's constant \f$q < 0\f$.

Reference: Carlson 1992 (http://adsabs.harvard.edu/abs/1992MaCom..59..165C)
@param r1 lower bound of integral \f$r_a\f$
@param r2 upper bound of integral \f$r_b\f$
*/
double geoinf::caltaur_negq(const double r1, const double r2) const {
    double f2, g2;
    double xi2, eta2, theta1, theta2, zeta1, zeta2, M, M2;
    double Delta, Delta_plus, Delta_minus, L2_plus, L2_minus, delta12sqr, delta11sqr, delta22sqr;
    double temp1, r1sqr = r1 * r1, r2sqr = r2 * r2;

    if (q >= 0.) {
        throw("geoinf::caltaur_inf_negq: q>=0");
    }

    f2 = rr3 * rr3 + rr4 * rr4;
    g2 = -2. * rr3;
    xi2 = std::sqrt(f2 + g2 * r2 + r2sqr);
    eta2 = std::sqrt(f2 + g2 * r1 + r1sqr);
    temp1 = (r2 - r1) * (r2 - r1);

    if (nrr == 0) {
        double f1, g1;
        double xi1, eta1;
        f1 = rr1 * rr1 + rr2 * rr2;
        g1 = -2. * rr1;


        eta1 = std::sqrt(f1 + g1 * r1 + r1sqr);
        xi1 = std::sqrt(f1 + g1 * r2 + r2sqr);

        theta1 = xi1 * xi1 + eta1 * eta1 - temp1;
        theta2 = xi2 * xi2 + eta2 * eta2 - temp1;
        zeta1 = std::sqrt(2. * xi1 * eta1 + theta1);
        zeta2 = std::sqrt(2. * xi2 * eta2 + theta2);
        M = zeta1 * zeta2 / (r2 - r1);
        M2 = M * M;

        delta12sqr = 2. * (f1 + f2) - g1 * g2;
        delta11sqr = 4. * f1 - g1 * g1;
        delta22sqr = 4. * f2 - g2 * g2;
        Delta = std::sqrt(delta12sqr * delta12sqr - delta11sqr * delta22sqr);
        Delta_plus = delta12sqr + Delta;
        Delta_minus = delta12sqr - Delta;
        L2_plus = M2 + Delta_plus;
        L2_minus = M2 + Delta_minus;
    } else {
        double X1, Y1, X4, Y4;
        double c14sqr, c11, c44;

        X1 = std::sqrt(r2 - rr1);
        Y1 = std::sqrt(r1 - rr1);
        X4 = std::sqrt(r2 - rr2);
        Y4 = std::sqrt(r1 - rr2);

        M2 = (X1 * Y4 + X4 * Y1) * (X1 * Y4 + X4 * Y1) * (2. * xi2 * eta2 + 2. * f2 + g2 * (r1 + r2) + 2. * r1 * r2) / temp1;
        std::cout << "f2 = " << f2 << ", g2 = " << g2 << std::endl;

        c14sqr = 2. * f2 + g2 * (rr1 + rr2) + 2. * rr1 * rr2;
        c11 = std::sqrt(2. * (f2 + g2 * rr1 + rr1 * rr1));
        c44 = std::sqrt(2. * (f2 + g2 * rr2 + rr2 * rr2));

        L2_plus = M2 + c14sqr + c11 * c44;
        L2_minus = M2 + c14sqr - c11 * c44;
    }

    return 4. * rf(M2, L2_minus, L2_plus);
}

/*!
\brief Calculate the integral over radius along the whole geodesic, given position, sign of \f$k^r\f$ at the starting point, the position of end point,
and number of encounting turning points.

@param re starting radius
@param rf end radius
@param kre sign of \f$k^r\f$ at the beginning
@param nrr2 number of the photon encounting the turning point at the lower boundary
@param nrr3 number of the photon encounting the turning point at the upper boundary
*/
double geoinf::caltaur_all(const double re, const double rf, const double kre, const int nrr2, const int nrr3) {
    double taur, taurer2, taurer3, taurfr2, taurfr3;
    if (q >= 0.) {
        if (locatere(re)) {
            // in this case at most one turning point
            taur = (nrr3 > 0) ? caltaur_r1(re) + caltaur_r1(rf) : std::abs(caltaur_r1(re) - caltaur_r1(rf));
        } else {
            // the photon may bounce several times between r2 and r3, provided that r3 is outside of the BH event horizon; otherwise there
            // is at most one turning point at r2
            if (rr3 > rh) {
                taurer2 = caltaur_r2(re);
                taurer3 = taur2r3 - taurer2;
                taurfr2 = caltaur_r2(re);
                taurfr3 = taur2r3 - taurfr2;

                if (kre >= 0.) {
                    taur = (nrr2 == nrr3) ? taurer2 - taurfr2 + (double)(2 * nrr3) * taur2r3 : taurer2 + taurfr2 + (double)(2 * nrr3) * taur2r3;
                } else {
                    taur = (nrr2 == nrr3) ? taurer3 - taurfr3 + (double)(2 * nrr2) * taur2r3 : taurer3 + taurfr3 + (double)(2 * nrr2) * taur2r3;
                }
            } else {
                taur = (nrr2 > 0) ? caltaur_r2(re) + caltaur_r2(rf) : std::abs(caltaur_r2(re) - caltaur_r2(rf));
            }
        }
    } else {
        taur = std::abs(caltaur_inf_negq(re) - caltaur_inf_negq(rf));
    }

    return taur;
}

/*!
\brief Calculate the integral along r between \f$r_e\f$ and the corresponding turning point.

The integral is:
- if \f$q \geq 0\f$, and \f$R(r)=0\f$ has four real roots and \f$r_e \geq r_1\f$: \f$\int_{r_1}^{r_e} dr/\sqrt{R(r)}\f$
- if \f$q \geq 0\f$, and \f$R(r)=0\f$ has four real roots and \f$r_3 \leq r_e \leq r_2\f$: \f$\int_{r_e}^{r_2} dr/\sqrt{R(r)}\f$
- if \f$q \geq 0\f$, and \f$R(r)=0\f$ has two real roots: \f$\int_{r_1}^{r_e} dr/\sqrt{R(r)}\f$
- if \f$q < 0\f$: \f$\int_{r_e}^{\infty} dr/\sqrt{R(r)}\f$
@param re radius
*/
double geoinf::caltaue(const double re) const {
    double taue;
    if (q >= 0.) {
        if (locatere(re)) {
            taue = caltaur_r1(re);
        } else {
            taue = caltaur_r2(re);
        }
    } else {
        taue = caltaur_inf_negq(re);
    }
    return taue;
}

/*!
\brief Evaluate r-part in the integral of the affine parameter \f$\lambda\f$ between \f$r_e\f$ and the corresponding turning point.
The turning point is the same as \ref caltaue.

The affine parameter: 
\f{equation}
\label{eq:lambda}
 \lambda = \int \frac{r^2dr}{\sqrt{R(r)}} + a^2 \int \frac{\mu^2 d\mu}{\sqrt{M(\mu)}}.
\f}
@param re radius
*/
double geoinf::caldlambda_re(const double re) const {
    double dlambdae;
    if (q >= 0.) {
        if (locatere(re)) {
            dlambdae = caldlambda_r1(re);
        } else {
            dlambdae = caldlambda_r2(re);
        }
    } else {
        dlambdae = 0.;
    }
    return dlambdae;
}

/*!
\brief Calculate the integral over \f$\mu\f$ along the whole geodesic, given position, sign of \f$k^\theta\f$ at the starting point,
the position of end point, and number of encounting turning points.

@param mue starting radius
@param muf end radius
@param kmue sign of \f$k^r\f$ at the beginning
@param nmuminus number of the photon encounting the turning point at the lower boundary
@param nmuplus number of the photon encounting the turning point at the upper boundary
*/
double geoinf::calximu_all(const double mue, const double muf, const double kmue, const int nmuminus, const int nmuplus) {
    double ximu, xiall;
    double ximue_muplus, ximue_muminus, ximuf_muplus, ximuf_muminus;
    
    bool lower = false;
    if (q < 0.) {
        lower = (mue < 0.);
    } else if (q == 0.) {
        lower = (mue < 0.) || ((mue == 0.) && (kmue > 0.));
    }
    
    if (lower) 
        return calximu_all(-1. * mue, -1. * muf, -1. * kmue, nmuminus, nmuplus);

    ximue_muplus = calximu(mue);
    ximuf_muplus = calximu(muf);

    xiall = (q > 0.) ? 2. * taumu0 : taumu0;

    ximue_muminus = xiall - ximue_muplus;
    ximuf_muminus = xiall - ximuf_muplus;

    if (kmue <= 0.) {
        ximu = (nmuplus == nmuminus) ?
            ximue_muplus - ximuf_muplus + (double)(2 * nmuminus) * xiall :
            ximue_muplus + ximuf_muplus + (double)(2 * nmuminus) * xiall;
    } else {
        ximu = (nmuplus == nmuminus) ?
            ximue_muminus - ximuf_muminus + (double)(2 * nmuplus) * xiall:
            ximue_muminus + ximuf_muminus + (double)(2 * nmuplus) * xiall;
    }
    return ximu;
}

/*!
\brief Evaluate \f$\int_\mu^{\mu_+}d\mu/\sqrt{M(\mu)}\f$. For \f$M(\mu)\f$ see \ref kerr_muintegrand.

Reference: Li+2005, Dexter & Agol 2009.
@param mu cosine of polar angle.
*/
double geoinf::calximu(const double mu) const {
    double x, tau;
    if ((mu > muplus) && (std::abs(mu - muplus) <= DEPS)) {
        return 0.;
    }
    
    if (q <= 0. && mu <= 0.)
        return calximu(-1. * mu);
    
    if (q < 0.) {
        if (std::abs(mu - muminus) <= DEPS) {
            return taumu0;
        }
    }

    if (mu > muplus) {
        std::cout << "mu = " << mu << ", muplus = " << muplus << ", minus = " << muminus <<", geo.q = " << q << std::endl;
        throw("geoinf::calximu(): mu > muplus");
    }
    if (q < 0. && mu < muminus) {
        std::cout << "mu = " << mu << ", muplus = " << muplus << ", minus = " << muminus <<", geo.q = " << q << std::endl;
        throw("geoinf::calximu(): mu < muminus");
    }

    if (spin0) {
        tau = (M_PI * 0.5 - std::asin(mu / muplus)) * taumunorm;
    } else {
        if (q >= 0.) {
            x = mu / muplus;
        } else {
            x = std::sqrt((mu * mu - muminus2) / (muplus2 - muminus2));
        }
        if (q == 0.) {
            double inttemp = std::sqrt(1. - x * x);
            tau = taumunorm * std::log((1. + inttemp) / x);
        } else {
            tau = taumunorm * jacobi_icn(x, mmu);
        }
    }

    return tau;
}

/*
double geoinf::calximu_lowen(const double mu) const {
    double muplusnow = std::sqrt(q / (l * l + q));
    if (!lowen) {
        std::cerr << "Inside calximu_lowen()" << std::endl;
        throw("Not lowen case!");
    }

    return M_PI - std::asin(mu / muplusnow);
}
*/

/*!
\brief Evaluate \f$\mu\f$-part in the integral of the affine parameter \f$\lambda\f$ between \f$\mu_e\f$ and \f$\mu_+\f$.

For expression of \f$\lambda\f$ see \ref caldlambda_re. 
This function evaluates the integral \f$\int_{\mu_e}^{\mu_+} \mu^2 d\mu /M(\mu)\f$, where \f$\mu_+\f$.

@param mu \f$\mu_e\f$
*/
double geoinf::caldlambda_muplus(const double mu) const {
    double dlambda, x;
    if (spin0)
        return a * a * 0.5 * std::sqrt(q * q / (l * l + q) / (l * l + q) / (l * l + q)) *
        (M_PI * 0.5 - std::asin(mu / muplus) + mu / muplus * std::sqrt(1. - mu * mu / muplus / muplus));
    //std::cout << "mu = " << mu << std::endl;
    
    if (std::abs(mu - muplus) <= DEPS) {
        return 0.;
    }
    
    if (q <= 0. && mu <= 0.)
        return caldlambda_muplus(-1. * mu);

    if (q < 0.) {
        if (std::abs(mu - muminus) <= DEPS) {
            return dlambda_mu0;
        }
    }

    if (mu > muplus)
        throw("geoinf::caldlambda_muplus(): mu > muplus");
    
    if (q < 0. && mu < muminus) {
        std::cout << "mu < muminus" << std::endl;
        throw("geoinf::caldlambda_muplus(): mu < muminus");
    }

    if (q != 0.) {
        if (q > 0.) {
            x = mu / muplus;
        } else {
            x = std::sqrt((mu * mu - muminus2) / (muplus2 - muminus2));
        }
        dlambda = elliptic_e_cos(x, mmu) / taumunorm;
    } else {
        dlambda = a * std::sqrt(muplus2 - mu * mu);
    }

    return dlambda;
}

/*
double geoinf::caldlambda_muplus_lowen(const double mu) const {
    double muplusnow = std::sqrt(q / (l * l + q));
    double t = std::abs(mu) / muplus;
    double inte;
    if (!lowen) {
        std::cerr << "Inside calximu_lowen()" << std::endl;
        throw("Not lowen case!");
    }

    inte = 0.5 * (0.5 * M_PI - std::asin(t) + t * std::sqrt(1. - t * t));
    if (mu < 0)
        inte = M_PI * 0.5 - inte;
    return muplusnow * muplusnow * inte * a * a;
}
*/

/*!
\brief Similar with \ref calximu_all, but for \f$\lambda\f$. 

@param mue starting point
@param muf end point
@param kmue sign of \f$k^\theta\f$ at the beginning
@param nmuminus number of the photon encounting the turning point at the lower boundary
@param nmuplus number of the photon encounting the turning point at the upper boundary
@param ximue the \f$\mu-\f$integral of \f$\lambda\f$ between \f$\mu_e\f$ and \f$\mu_+\f$
@param inte \f$\int_{\mu_e}^{\mu_+} dr/\sqrt{R(r)}\f$.
*/
double geoinf::caldlambda_muall(const double mue, const double muf, const double kmue, 
    const int nmuminus, const int nmuplus, const double ximue, const double inte) const {
        
    double ximu, xiall;
    double ximue_muminus, ximuf_muplus, ximuf_muminus;

    if (spin0)
        return 0.;
    
    bool lower = false;
    if (q < 0.) {
        lower = (mue < 0.);
    } else if (q == 0.) {
        lower = (mue < 0.) || ((mue == 0.) && (kmue > 0.));
    }
    
    if (lower)
        return caldlambda_muall(-1. * mue, -1. * muf, -1. * kmue, nmuminus, nmuplus,
            ximue, inte);
    //std::cout << "mue = " << mue << ", muf = " << muf << ", kmue = " << kmue << std::endl;

    ximuf_muplus = caldlambda_muplus(muf);
    
    xiall = (q > 0.) ? 2. * dlambda_mu0 : dlambda_mu0;
    ximue_muminus = xiall - ximue;
    ximuf_muminus = xiall - ximuf_muplus;
    
    if (kmue <= 0.) {
        ximu = (nmuplus == nmuminus) ?
            ximue - ximuf_muplus + (double)(2 * nmuminus) * xiall :
            ximue + ximuf_muplus + (double)(2 * nmuminus) * xiall;
    } else {
        ximu = (nmuplus == nmuminus) ?
            ximue_muminus - ximuf_muminus + (double)(2 * nmuplus) * xiall:
            ximue_muminus + ximuf_muminus + (double)(2 * nmuplus) * xiall;
    }


    if (q > 0.)
        ximu -= (a2 * muminus2 * inte);
    return ximu;
}

/*! \brief Evaluate \f$r-\f$ part of \f$\lambda\f$, between \f$r_1\f$ and \f$r\f$ where \f$r_1\f$ is the largest root of \f$R(r)=0\f$.
 This function works only when \f$q\geq 0\f$ and \f$r\geq r_1\f$.

For expression of \f$\lambda\f$, see \ref kerr_lambdaintegrand_r.

Reference: 
    - Four real roots: Eq. 2.33 of Carlson 1988 (http://www.ams.org/mcom/1988-51-183/S0025-5718-1988-0942154-7/);
    - Two real roots: Carlson 1991 (http://adsabs.harvard.edu/abs/1991MaCom..56..267C)
@param r radius
*/
double geoinf::caldlambda_r1(const double r) const {
    double tau, I2, A;
    if (std::abs(r - rr1) <= DEPS)
        return 0.;
    if (r < rr1)
        throw("geoinf::caltaur_r1(): r < r1");
    if (q < 0.)
        throw("geoinf::caltaur_r1(): q < 0");

    if (nrr == 4) {
        double one_over_rmrr1, U12sqr, U13sqr, U14sqr;
        one_over_rmrr1 = 1. / (r - rr1);
        U12sqr = (r - rr2) * (rr1 - rr3) * (rr1 - rr4) * one_over_rmrr1;
        U13sqr = (r - rr3) * (rr1 - rr2) * (rr1 - rr4) * one_over_rmrr1;
        U14sqr = (r - rr4) * (rr1 - rr2) * (rr1 - rr3) * one_over_rmrr1;
        I2 = 2. * (rr1 - rr2) * (rr1 - rr3) * rd(U12sqr, U13sqr, U14sqr) / 3.;
        A = std::sqrt((r - rr1) * (r - rr2) * (r - rr3) / (r - rr4));
        tau = 0.5 * (rr2 - rr4) * (rr3 - rr4) * I2 + A;

    } else {
        double f, g, xi;
        double c14sqr, M2, L2_plus, L2_minus, U;

        f = u * u + v * v;
        g = -2. * u;

        xi = std::sqrt(f + g * r + r * r);
        c14sqr = 2. * f + g * (rr1 + rr2) + 2. * rr1 * rr2;

        U = std::sqrt((r - rr2) / (r - rr1)) * rr3;
        M2 = (rr1 - rr2) * (2. * xi * rr3 + 2. * f + g * (r + rr1) + 2. * r * rr1) / (r - rr1);
        L2_plus = M2 + c14sqr + 2. * rr3 * rr4;
        L2_minus = M2 + c14sqr - 2. * rr3 * rr4;
        I2 = (2. * rr3 / 3. / rr4) * (4. * (c14sqr + 2. * rr3 * rr4) * rd(M2, L2_minus, L2_plus) -
            6. * rf(M2, L2_minus, L2_plus) + 3. / U);
        A = std::sqrt((r - rr1) / (r - rr2)) * xi;
        tau = 0.5 * rr4 * rr4 * I2 + A;
    }

    return tau;
}

/*! \brief Evaluate \f$r-\f$ part of \f$\lambda\f$, between \f$r\f$ and \f$r_2\f$ where \f$r_2\f$ is the second largest root of \f$R(r)=0\f$.
 This function works only when \f$R(r)=0\f$ has four real roots and and \f$r_3 \leq r\leq r_2\f$.
 
Reference: Carlson 1988.
@param r radius
*/
double geoinf::caldlambda_r2(const double r) const {
    double tau, I2, A;
    if (q < 0.)
        throw("geoinf::caldlambda_r2(): q < 0");
    if (nrr != 4)
        throw("geoinf::caldlambda_r2(): nrr != 4");
    if (std::abs(r - rr2) <= DEPS)
        return 0.;
    if (rr3 - r > DEPS)
        throw("geoinf::caltaur_r1(): r < r3");
    if (r > rr2) {
        std::cout << "rr1 = " << rr1 << ", rr2 = " << rr2 << std::endl;
        std::cout << "r = " << r << std::endl;
        throw("geoinf::caltaur_r1(): r > r2");
    }

    double one_over_rmrr2, U12sqr, U13sqr, U14sqr;
    one_over_rmrr2 = 1. / (rr2 - r);
    U12sqr = (rr1 - r) * (rr2 - rr3) * (rr2 - rr4) * one_over_rmrr2;
    U13sqr = (r - rr3) * (rr1 - rr2) * (rr2 - rr4) * one_over_rmrr2;
    U14sqr = (r - rr4) * (rr1 - rr2) * (rr2 - rr3) * one_over_rmrr2;
    I2 = 2. * (rr1 - rr2) * (rr2 - rr3) * rd(U12sqr, U13sqr, U14sqr) / 3.;
    A = std::sqrt((rr1 - r) * (rr2 - r) * (r - rr3) / (r - rr4));
    tau = -0.5 * (rr1 - rr4) * (rr3 - rr4) * I2 - A;

    return tau;
}

/*!
\brief Calculate r-part of lambda integral between \f$r_a\f$ and \f$r_b\f$, in the case of q < 0.

References: Carlson 1991; Eq. 2.44 of Carlson 1992 (http://adsabs.harvard.edu/abs/1992MaCom..59..165C)

@param r1 the lower bound of the integral \f$r_a\f$
@param r2 the upper bound of the integral \f$r_b\f$
*/
double geoinf::caldlambda_negq(const double r1, const double r2) const {
    double M2, L2_plus, L2_minus, U;
    double I1, I2, A;
    if (q >= 0.)
        throw("geoinf::caldlambda_inf_negq(): q < 0");
    double f2, g2, xi2, eta2, tau, temp1;
    f2 = u * u + v * v;
    g2 = -2. * u;
    xi2 = std::sqrt(f2 + g2 * r2 + r2 * r2);
    eta2 = std::sqrt(f2 + g2 * r1 + r1 * r1);
    temp1 = (r2 - r1) * (r2 - r1);

    if (nrr == 2) {
        double c11, c44, c14sqr, Y1, Y4, rfval;

        double X1, X4;

        X1 = std::sqrt(r2 - rr1);
        Y1 = std::sqrt(r1 - rr1);
        X4 = std::sqrt(r2 - rr2);
        Y4 = std::sqrt(r1 - rr2);

        M2 = (X1 * Y4 + X4 * Y1) * (X1 * Y4 + X4 * Y1) * (2. * xi2 * eta2 + 2. * f2 + g2 * (r1 + r2) + 2. * r1 * r2) / temp1;

        c14sqr = 2. * f2 + g2 * (rr1 + rr2) + 2. * rr1 * rr2;
        c11 = std::sqrt(2. * (f2 + g2 * rr1 + rr1 * rr1));
        c44 = std::sqrt(2. * (f2 + g2 * rr2 + rr2 * rr2));

        L2_plus = M2 + c14sqr + c11 * c44;
        L2_minus = M2 + c14sqr - c11 * c44;

        U = (X1 * X4 * eta2 + Y1 * Y4 * xi2) / (r2 - r1);

        rfval = rf(M2, L2_minus, L2_plus);

        I1 = 4. * rfval;
        I2 = (2. * c11 / 3. / c44) * (4. * (c14sqr + c11 * c44) * rd(M2, L2_minus, L2_plus) -
            6. * rfval + 3. / U) + 2. * X1 * Y1 / X4 / Y4 / U;
        A = X1 * xi2 / X4 - Y1 * eta2 / Y4;

        tau = 0.5 * ((rr2 - u) * (rr2 - u) + v * v) * I2 + (rr1 * rr1 - 0.5 * ((rr1 - u) * (rr1 - u) + v * v)) * I1 +
            A;
    } else {
        double f1, g1;
        double xi1, eta1, xi1_prime, eta1_prime;
        double M, r1sqr = r1 * r1, r2sqr = r2 * r2, theta1, theta2, zeta1, zeta2, delta12sqr, delta11sqr, delta22sqr, Delta, Delta_minus, Delta_plus;
        double U, G, B, rfval, Sigma, phi0, mu0, T0, V0sqr, a0, b0sqr, X0, H0, Omega0sqr, S;

        f1 = rr1 * rr1 + rr2 * rr2;
        g1 = -2. * rr1;

        xi1 = std::sqrt(f1 + g1 * r2 + r2sqr);
        eta1 = std::sqrt(f1 + g1 * r1 + r1sqr);

        xi1_prime = (g1 + 2. * r2) / 2. / xi1;
        eta1_prime = (g1 + 2. * r1) / 2. / eta1;

        phi0 = g1 - g2;
        X0 = (xi1_prime * xi2 + eta1_prime * eta2) / (r1 - r2);
        mu0 = 1. / (xi1 * eta1);

        theta1 = xi1 * xi1 + eta1 * eta1 - temp1;
        theta2 = xi2 * xi2 + eta2 * eta2 - temp1;
        zeta1 = std::sqrt(2. * xi1 * eta1 + theta1);
        zeta2 = std::sqrt(2. * xi2 * eta2 + theta2);
        M = zeta1 * zeta2 / (r2 - r1);
        M2 = M * M;

        delta12sqr = 2. * (f1 + f2) - g1 * g2;
        delta11sqr = 4. * f1 - g1 * g1;
        delta22sqr = 4. * f2 - g2 * g2;
        Delta = std::sqrt(delta12sqr * delta12sqr - delta11sqr * delta22sqr);
        Delta_plus = delta12sqr + Delta;
        Delta_minus = delta12sqr - Delta;
        L2_plus = M2 + Delta_plus;
        L2_minus = M2 + Delta_minus;

        U = (xi1 * eta2 + xi2 * eta1) / (r2 - r1);
        S = (M2 + delta12sqr) / 2. - U * U;
        T0 = mu0 * S + 2.;
        V0sqr = mu0 * mu0 * (S * S + delta11sqr * U * U);
        Omega0sqr = M2 + delta11sqr;
        a0 = S * Omega0sqr / U + 2. * delta11sqr * U;
        b0sqr = (S * S / (U * U) + delta11sqr) * (Omega0sqr * Omega0sqr);

        rfval = rf(M2, L2_minus, L2_plus);
        G = 2. * Delta * Delta_plus * rd(M2, L2_minus, L2_plus) / 3. + Delta / 2. / U +
            (delta12sqr * theta1 - delta11sqr * theta2) / 4. / xi1 / eta1 / U;
        B = xi1_prime * xi2 - eta1_prime * eta2;
        Sigma = G - Delta_plus * rfval + B;
        H0 = delta11sqr * phi0 *(rj(M2, L2_minus, L2_plus, Omega0sqr) / 3. + rc(a0 * a0, b0sqr) / 2.) -
            X0 * rc(T0 * T0, V0sqr);

        tau = Sigma + g1 * g1 * rfval - (g1 + g2) * H0;
    }

    return tau;
}

/*!
\brief similar with \ref caltaur_all, but for \f$\lambda\f$.
@param re starting radius
@param rf end radius
@param kre sign of \f$k^r\f$ at the beginning
@param nrr2 number of the photon encounting the turning point at the lower boundary
@param nrr3 number of the photon encounting the turning point at the upper boundary
@param dlambdare \f$r-\f$ part of \f$\lambda\f$ integral between \f$r_e\f$ and the corresponding turning point; see \ref caldlambda_re
@param inte \f$r-\f$ integral between \f$r_e\f$ and the corresponding turning point; see \ref caltaue
 */
double geoinf::caldlambda_rall(const double re, const double rf, const double kre, const int nrr2, 
    const int nrr3, const double dlambdare, const double inte) {
    
    double taur, taurf_r1, I1part;
    if (q >= 0.) {
        if (locatere(re)) {
            // in this case at most one turning point
            taurf_r1 = caldlambda_r1(rf);
            if (kre >= 0.)
                taur = taurf_r1 - dlambdare;
            else
                taur = (nrr3 > 0) ? dlambdare + taurf_r1 : dlambdare - taurf_r1;

            if (nrr == 4) {
                I1part = (rr1 * rr1 - 0.5 * (rr1 - rr2) * (rr1 - rr3)) * inte;
            } else {
                I1part = (rr1 * rr1 - 0.5 * rr3 * rr3) * inte;
            }
            taur += I1part;

        } else {
            double taurer3, taurfr2, taurfr3;
            taurer3 = dlambda_r2r3 - dlambdare;
            taurfr2 = caldlambda_r2(rf);
            taurfr3 = dlambda_r2r3 - taurfr2;

            if (rr3 > rh) {
                if (kre >= 0.) {
                    taur = (nrr2 == nrr3) ? dlambdare - taurfr2 + (double)(2 * nrr3) * taur2r3 : dlambdare + taurfr2 + (double)(2 * nrr3) * taur2r3;
                } else {
                    taur = (nrr2 == nrr3) ? taurer3 - taurfr3 + (double)(2 * nrr2) * taur2r3 : taurer3 + taurfr3 + (double)(2 * nrr2) * taur2r3;
                }
            } else {
                if (kre >= 0.)
                    taur = (nrr2 > 0) ? dlambdare + taurfr2 : dlambdare - taurfr2;
                else
                    taur = taurfr2 - dlambdare;
            }
            I1part = (rr2 * rr2 + 0.5 * (rr1 - rr2) * (rr2 - rr3)) * inte;

            taur += I1part;
        }
    } else {
        taur = (kre >= 0.) ? caldlambda_negq(re, rf) : caldlambda_negq(rf, re);
    }
    return taur;
}

/*!
\brief invert \f$\mu-\f$integral to solve for \f$\mu\f$.

Reference: Dexter & Agol 2009
//! @param tau \f$\mu-\f$integral
*/
double geoinf::solvemu(const double tau) const {
//    std::cout << "taumu0 = " << taumu0 << ", muplus = " << muplus << std::endl;
    double x, mu, xiall;
    xiall = (q > 0. ? 2. * taumu0 : taumu0);
    if (tau > xiall) {
        throw("muinttoolarge");
    }

    if (spin0) {
        mu = muplus * std::cos(tau / taumunorm);
    } else {
        if (q == 0.) {
            double tempnow = std::exp(tau / taumunorm);
            x = 2. * tempnow / (tempnow * tempnow + 1.);
        } else {
            x = jacobi_cn(tau / taumunorm, mmu);
        }

        if (q >= 0.) {
            mu = x * muplus;
        } else {
            mu = std::sqrt(x * x * (muplus2 - muminus2) + muminus2);
        }
    }

    return mu;
}

/*
//! find turning points
void geoinf::findturns() {
    double tau = std::abs(taumuinf - taurinf);
    if (nmuturn == 0) {
        if (nrturn == 1)
            muturn = solvemu(taumuinf + taurinf);
    } else {
        rturn = ((tau <= 1e-6) ? rr1 : solver_r1(tau));
        if (nrturn == 1)
            muturn = ((tau <= 1e-6) ? muplus : solvemu(tau));
    }
    turnsfound = true;
}
*/

/*
//! calculate the integration of geodesic along
double deltatau(double mu, const slabgeo & slab) {
    double taumul, taurl;
    //! starting with the simplest case with no turning point
    if (slab.geo->nmuturn == 0 && slab.geo->nrturn == 0) {
        taumul = slab.geo->calximu(mu) - slab.geo->taumuinf;
        taurl = slab.geo->taurinf - slab.geo->caltaur_r1(slab.h / mu);
    } else if (slab.geo->nmuturn == 1 && slab.geo->nrturn == 0) {
        taumul = slab.geo->calximu(mu);
        taurl = slab.geo->caltaur_r1(slab.geo->rturn) - slab.geo->caltaur_r1(slab.h / mu);
    } else if (slab.geo->nmuturn == 1 && slab.geo->nrturn == 1) {
        if (slab.geo->taumuinf > slab.geo->taurinf) {
            taumul = slab.geo->calximu(slab.geo->muturn) - slab.geo->calximu(mu);
            taurl = slab.geo->caltaur_r1(slab.h / mu);
        } else {
            taumul = slab.geo->calximu(mu);
            taurl = slab.geo->caltaur_r1(slab.geo->rturn) - slab.geo->caltaur_r1(slab.h / mu);
        }
    } else {
        throw("nosuchcase");
    }
    return taumul - taurl;
}
*/

/*
//! calculate the position where the geodesic reaches a slab at height h
//! @param h height of the slab
double calslabrad(geoinf & geo, const double h) {
    // starting with the simplest case with no turning point
    double mu = 0., ROUT = 1000.;
    double mumin, mumax;
    double rmax, munow, intmu, intr;
    double region1hmin;
    double hturn1, hturn2, muturn2, rturn2;
    long res;
    std::vector<double> x(4), kmu(4);
    // relative positoin of r, mu turning points with respect to
    double kmunow;
    raytrace_data rtd;
    double precision = 0.01, step = 0.01;
    bool found;

    if (!geo.areturnsfound()) {
        std::cerr << "Turning points not found yet!" << std::endl;
        return 0.;
    }
    slabgeo slab(geo, h);

    if (geo.nmuturn == 0 && geo.nrturn == 0) {
        //std::cout << "case1" << std::endl;
        mumin = h / ROUT;
        mumax = std::min(geo.muobs, h / geo.rr1 - 1e-6);
        res = rtbis<slabgeo>(mumin, mumax, 1e-4, &deltatau, slab, mu);
        if (res == 0) {
            throw("case1notfound");
        }
        // we do not check event horizon as mu <= mumax < 1, r = h * mu > h.
        // As long as $h > r_h$, we are guaranteed that r > r_h.
    }
    // above checked

    if (geo.nmuturn == 1 && geo.nrturn == 0) {
        //std::cout << "case2" << std::endl;
        // In this case the intersection must be before the turning point.
        found = false;
        // should be checked that geo.rturn is solved without problem
        if (geo.rturn * geo.muobs < h) {
            // for this part we are calculating integral with respective
            // to infinity; and rmax > r_h
            rmax = h / geo.muplus;
            intr = geo.caltaur_r1(rmax);
            intmu = geo.taumuinf + geo.taurinf - intr;
            munow = geo.solvemu(intmu);

            // calculate four momentum
            x = {0., rmax, munow, 0.};
            photon_momentum(geo.a, rmax, munow, geo.l, geo.q, -1., -1., kmu.data());

            // find the intersection point using step-by-step raytracing
            // again x[1] >= h > rh
            raytrace_prepare(geo.a, x.data(), kmu.data(), NULL, precision, 0, &rtd);
            kmunow = kmu[2];
            while(kmunow * kmu[2] > 0) {
                //raytrace1(x.data(), kmu.data(), &step, &rtd);
				raytrace(x.data(), kmu.data(), NULL, &step, &rtd);
                if (x[1] * x[2] <= h) {
                    found = true;
                    break;
                }
                kmunow = kmu[2];
            }
        }

        if (found) {
            mu = x[2];
            // in the future refining the result
        } else {
            if (geo.muplus * geo.rturn < h) {
                throw("case2notfound");
            } else if(geo.muplus * geo.rturn == h) {
                mu = geo.muplus;
            } else {
                // if muplus is properly calculated, then it it guaranteed that
                // h / mu > r_h.
                mumin = h / geo.rturn;
                geo.solvere();
                mumax = std::min(geo.muplus, h / geo.re);
                res = rtbis<slabgeo>(mumin, mumax, 1e-4, &deltatau, slab, mu);
            }
        }
    }

    // whether rturn is calculated properly
    if (geo.nmuturn == 0 && geo.nrturn == 1) {
        //std::cout << "case3" << std::endl;
        if (geo.muturn * geo.rr1 < h) {
            mumin = geo.muturn;
            mumax = geo.muobs;
            res = rtbis<slabgeo>(mumin, mumax, 1e-4, &deltatau, slab, mu);
            if (res == 0) throw("case3notfound");
            // r > h / muobs > rh
        } else if (geo.muturn * geo.rr1 == h) {
            mu = geo.muturn;
        } else {
            found = false;
            if (geo.rr1 < geo.rh)
                throw("withinrh");
            x = {0., geo.rr1, geo.muturn, 0.};
            photon_momentum(geo.a, geo.rr1, geo.muturn, geo.l, geo.q, 1., 1., kmu.data());
            kmu[1] = 0.;
            raytrace_prepare(geo.a, x.data(), kmu.data(), NULL, precision, 0, &rtd);
            while(x[1] <= geo.re) {
                //raytrace1(x.data(), kmu.data(), &step, &rtd);
                raytrace(x.data(), kmu.data(), NULL, &step, &rtd);
                if (x[1] * x[2] <= h) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                throw("case3notfound");
            } else {
                mu = x[2];
            }
        }
    }

    if (geo.nmuturn == 1 && geo.nrturn == 1) {
        //std::cout << "case4" << std::endl;
        if (geo.taurinf >= geo.taumuinf) {
            region1hmin = geo.rturn * geo.muplus;
            hturn1 = geo.muplus * geo.rturn;
            hturn2 = geo.rr1 * geo.muturn;
            rturn2 = geo.rr1;
            muturn2 = geo.muturn;
        } else {
            region1hmin = geo.rr1 * geo.muobs;
            hturn1 = geo.muturn * geo.rr1;
            hturn2 = geo.muplus * geo.rturn;
            muturn2 = geo.muplus;
            rturn2 = geo.rturn;
        }
        found = false;
        if (region1hmin < h) {
            // rmax > rh;
            rmax = h / geo.muobs;
            intr = geo.taurinf - geo.caltaur_r1(rmax);
            intmu = intr + geo.taumuinf;
            munow = geo.solvemu(intmu);
            x = {0., rmax, munow, 0.};
            photon_momentum(geo.a, rmax, munow, geo.l, geo.q, -1., -1., kmu.data());
            raytrace_prepare(geo.a, x.data(), kmu.data(), NULL, precision, 0, &rtd);
            found = false;
            kmunow = kmu[2];
            while(kmunow * kmu[2] > 0) {
                //raytrace1(x.data(), kmu.data(), &step, &rtd);
                raytrace(x.data(), kmu.data(), NULL, &step, &rtd);
                if (x[1] * x[2] <= h) {
                    found = true;
                    break;
                }
            }
            // r > rh now;

            if (found) {
                mu = x[2];
            } else {
                if (hturn1 < h) {
                    throw("case4notfound");
                } else if (hturn1 == h) {
                    mu = geo.muplus;
                    found = true;
                }
            }
        }

        if (!found) {
            if (hturn2 == h) {
                mu = muturn2;
                found = true;
            } else if (h > geo.rr1 * geo.muturn && h < geo.rturn * geo.muplus) {
                res = rtbis<slabgeo>(geo.muturn, std::min(geo.muplus, h / geo.rr1 - 1e-6), 1e-4, &deltatau, slab, mu);
                if (res == 0) {
                    throw("rootnotfound");
                }
                found = true;
            } else {
                x = {0., rturn2, muturn2, 0.};
                photon_momentum(geo.a, rturn2, muturn2, geo.l, geo.q, 1., 1., kmu.data());
                // p^r = 0 at r turning point;
                if (geo.taurinf >= geo.taumuinf)
                    kmu[1] = 0.;
                else
                    kmu[2] = 0.;
                raytrace_prepare(geo.a, x.data(), kmu.data(), NULL, precision, 0, &rtd);
                while(x[2] > 0.) {
                    //raytrace1(x.data(), kmu.data(), &step, &rtd);
                    raytrace(x.data(), kmu.data(), NULL, &step, &rtd);
                    //std::cout << "x = " << x << std::endl;
                    //std::cout << "kmu = " << kmu << std::endl;
                    if (x[1] * x[2] <= h) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    throw("case4notfound");
                }
                mu = x[2];
            }
        }
    }
    return mu;
};
*/

/*
double calslabrad_sim5(const geoinf & geo1, const double h, const double smax, double & rsign, double & tsign) {
    double intr;
    double rmax = std::sqrt(h * h + smax * smax);
    double signtheta, xnow, munow, ximu, r = 0.;
    std::vector<double> x(4), kmu(4);
    raytrace_data rtd;
    double step;
    bool found;

    if (rmax <= geo1.rr1) {
        throw("rmax <= r1");
    } else {
        intr = geo1.taurinf - geo1.caltaur_r1(rmax);
        if (geo1.beta >= 0.) {
            if (intr > geo1.taumuinf) {
                ximu = intr - geo1.taumuinf;
                signtheta = 1.;
            } else {
                ximu = geo1.taumuinf - intr;
                signtheta = -1.;
            }
        } else {
            ximu = geo1.taumuinf + intr;
            signtheta = 1.;
        }

        munow = geo1.solvemu(ximu);

        x = {0., rmax, munow, 0.};
        photon_momentum(geo1.a, rmax, munow, geo1.l, geo1.q, -1., signtheta, kmu.data());
        raytrace_prepare(geo1.a, x.data(), kmu.data(), NULL, 0.01, 0, &rtd);
        while (x[2] > 0){
            //raytrace1(x.data(), kmu.data(), &step, &rtd);
            raytrace(x.data(), kmu.data(), NULL, &step, &rtd);
            if (x[1] * x[2] < h) {
                //std::cout << "x = " << x << std::endl;
                //std::cout << "k = " << kmu << std::endl;
                kmu[1] *= -1.;
                kmu[2] *= -1.;
                rsign = (kmu[1] > 0 ? 1. : -1.);
                tsign = (kmu[2] > 0 ? 1. : -1.);
                photon_momentum(geo1.a, x[1], x[2], geo1.l, geo1.q, rsign, tsign, kmu.data());
                while (x[1] * x[2] < h) {
                    xnow = x[1];
                    munow = x[2];
                    //raytrace1(x.data(), kmu.data(), &step, &rtd);
                    raytrace(x.data(), kmu.data(), NULL, &step, &rtd);
                }
                r = (xnow + x[1]) / 2.;
                //std::cout << "r * mu = " << r * (munow + x[2]) / 2. << std::endl;
                found = true;
                break;
            }
        }
        if (!found) {
            r = 0.;
            rsign = 0.;
            tsign = 0.;
            throw("not found by sim5");
        }
    }
    return r;
}
*/

/*!
\brief To see if the photon is located in the allowed region.
@param r radius
@param mu \f$\mu\f$
*/
bool geoinf::isallowed(const double r, const double mu) const {
    bool allowed = true;
    if (mu > muplus) {
        allowed = false;
    }
    if (q >= 0.) {
        if (nrr == 4 && ((r < rr1 && r > rr2) || r < rr3)) {
            /*
            if (r < rr3) {
                std::cout << "Isallowed: r < rr3" << std::endl;
            } else {
                std::cout << "Isallowed: rr2 < r < rr1" << std::endl;
            }
            */
            allowed = false;
        }
        if (nrr == 2 && r < rr1) {
            allowed = false;
        }
    } else {
        // for q < 0: no turning point in r; \sqrt{M_} <= mu <= \sqrt{M_+}
        if (mu > muplus || mu < muminus) {
            std::cout << "Isallowed: mu > muplus || mu < muminus" << std::endl;
            allowed = false;
        }
    }
    return allowed;
}

/*! \brief See if the photon can reach the corona.
@param trid \ref tridgeo object; corona info
@param r \f$r\f$
@param mu \f$\mu\f$
@param kr sign of \f$k^r\f$
@param kmu sign of \f$k^\theta\f$
 */
bool geoinf::ispossible(const tridgeo & trid, const double r, const double mu, const double kr, const double kmu) {
    bool status = true;
    if (q >= 0.) {
        if ((trid.mumin > muplus) || (nrr == 2 && trid.rmax < rr1) ||
            (nrr == 4 && (((trid.rmax < rr1 && locatere(r)) || (trid.rmin > rr2 && !locatere(r))) || trid.rmax < rr3))) {
            status = false;
        }
        if (kr >= 0. && r >= rr1 && r > trid.rmax) {
            status = false;
        }
        if (kmu >= 0. && mu < trid.mumin) {
            status = false;
        }
        if (kr <= 0. && r < trid.rmin) {
            if ((nrr == 2 || (nrr == 4 && locatere(r))) && rr1 <= rh)
                status = false;
            if (nrr == 4 && !locatere(r) && rr3 <= rh)
                status = false;
        }
    } else {
        if ((trid.mumin > muplus) || (trid.mumax < muminus))
            status = false;
        if (kr >= 0. && r > trid.rmax)
            status = false;
        // even if there are two real roots, both of them will be negative
        if (kr <= 0. && r < trid.rmin)
            status = false;
    }

    return status;
}

/*
void geoinf::testturns() const {
    double taurturn, taurtoinf;
    double taumuturn, taumutoinf;
    if (!areturnsfound())
        std::cerr << "findturns() must be applied first!" << std::endl;
    if ((nmuturn + nrturn) > 0)
        std::cout << "Testing turnings points found:" << std::endl;
    if (nmuturn ==1 && nrturn == 0) {
        taurturn = caltaur_r1(rturn);
        taurtoinf = taurinf - taurturn;
        std::cout << taurtoinf << " == " << taumuinf << std::endl;
    }
    if (nmuturn == 0 && nrturn == 1) {
        taumuturn = calximu(muturn);
        taumutoinf = taumuturn - taumuinf;
        std::cout << taumutoinf << " == " << taurinf << std::endl;
    }
}
*/

/*
void geoinf::calmuobs_disk(const double kr, const double re) {
    double taur, taue, ximuobs;

    if (q <= 0.) {
        status = 0;
        return;
    }
    if (re < rr1) {
        status = 0;
        return;
    }

    taue = caltaur_r1(re);

    if (kr > 0.)
        taur = taurinf - taue;
    else if (kr == 0.)
        taur = taurinf;
    else
        taur = taurinf + taue;

    if (taur <= 2. * taumu0) {
        ximuobs = std::abs(taur - taumu0);
        muobs = solvemu(ximuobs);
        nmuturn = (taur > taumu0 ? 1 : 0);
        ismuobsset = true;
        status = 1;
    } else {
        status = 0;
    }
}
*/

/*! 
\brief When \f$q \geq 0\f$, tells which region is the photon in.

Return:
    - true: \f$r\geq r_1\f$
    - false: \f$r_3 \leq r \leq r_2\f$

@param r radius
*/
bool geoinf::locatere(const double r) const {
    return (q >= 0. && (r >= rr1 || std::abs(r - rr1) <= DEPS));
}

/*!
\brief Obsolete; calculate the inclination angle at infinity muobs. Now this work is done with \ref calfate.

@param kr r- component of photon wave vector.
@param r radius
@param mu mu = cos(theta)
@param ktheta theta- component of photon wave vector
*/
void geoinf::calmuobs(const double r, const double mu, const double kr, const double ktheta) {
    double taur, taue, ximuobs, ximu, kmufsign;

    if (q >= 0.) {
        if (!locatere(r)) {
            status = 0;
            ismuobsset = true;
            return;
        }

        taue = caltaur_r1(r);
        // special case; locating at mu turning point
        ximu = calximu(mu);

        if (kr > 0.)
            taur = taurinf - taue;
        else if (kr == 0.)
            taur = taurinf;
        else
            taur = taurinf + taue;

        if (ktheta >= 0.)
            ximuobs = ximu + taur;
        else
            ximuobs = std::abs(taur - ximu);
        //    std::cout << "taumu0 = " << taumu0 << std::endl;

        if (ximuobs < taumu0) {
            muobs = solvemu(ximuobs);
            status = 1;
            if (ktheta >= 0 || (ktheta < 0. && taur >= ximu))
                nmuturn = 1;
            else
                nmuturn = 0;
        } else {
            status = 0;
        }
    } else {
        // When q < 0, there is no turning point in r; hence one photon
        // with kr <= 0. in vacuum will never reach infinity
        if (kr <= 0.) {
            status = 0;
            ismuobsset = true;
            return;
        } else {
            int n1, n2;
            status = 1;
            taur = caltaur_inf_negq(r);
            // subtract integration times 2 * \int_mumin^mumax dmu / M(mu);
            muobs = solvemu_inte(taur, mu, ktheta, kmufsign, n1, n2);
            nmuturn = (kmufsign >= 0.) ? 1 : 0;
        }

    }
        // calculate r integral
    ismuobsset = true;
}

/*!
\brief Find out the final destination of the photon and save the result to the class variable *status*. 

Status:
   - 0: enters BH event horizon
   - 1: escapes to infinity
   - 2: strikes at the disc on the equatorial plane

Return:
    - 0: -1
    - 1: \f$\mu_\infty\f$
    - 2: \f$r_{\rm disc}\f$.
    
@param r starting point radius
@param mu starting point mu
@param kr sign of starting point \f$k^r\f$
@param kmu sign of starting point \f$k^\theta\f$ 
@param rin inner edge radius of the disc on the equatorial plane
*/
double geoinf::calfate(const double r, const double mu, const double kr, const double kmu, const double rin) {
    
    double taue, muinf, kmusign, krsign, ximudisc, rdisc, kmusigndisc;
    int nmuminus, nmuplus, nrr2, nrr3;
    
    //! when q < 0. the photon cannot reach the equatorial plane. 
    if (q < 0.) {
        //! no turning point in r; 
        if (kr >= 0.) {
            taue = caltaur_inf_negq(r);
            muinf = solvemu_inte(taue, mu, kmu, kmusign, nmuminus, nmuplus);
            status = 1;
        } else {
            status = 0;
        }
    } else if (q == 0.) {
        if (locatere(r)) {
            if (rr1 > rh) {
                status = 1;
            } else {
                status = (kr >= 0.) ? 1 : 0;
            }
            taue = (kr >= 0.) ? taurinf - caltaur_r1(r) : taurinf + caltaur_r1(r);
            if (status == 1) 
                muinf = solvemu_inte(taue, mu, kmu, kmusign, nmuminus, nmuplus);
        } else {
            status = 0;
        }
    } else {
        if (locatere(r)) {
            if (rr1 > rh) {
                taue = (kr >= 0.) ? taurinf - caltaur_r1(r) : taurinf + caltaur_r1(r);
                status = 1;
            } else {
                if (kr > 0.) {
                    taue = taurinf - caltaur_r1(r);
                    status = 1;
                } else {
                    taue = caltaur_r1(r) - caltaur_r1(rh);
                    status = 0;
                }
            }
        } else {
            if (rr3 > rh) {
                taue = 1e3;
                status = 1;
            }  else {
                taue = kr > 0. ? caltaur_r2(rh) + caltaur_r2(r) : caltaur_r2(rh) - caltaur_r2(r);
                status = 0;
            }
        }
        if (status == 1) 
            muinf = solvemu_inte(taue, mu, kmu, kmusign, nmuminus, nmuplus);
            
        ximudisc = ((mu >= 0. && kmu <= 0.) || (mu <= 0. && kmu >= 0.)) ?
            taumu0 + calximu(std::abs(mu)) : taumu0 - calximu(std::abs(mu));
        kmusigndisc = (mu >= 0.) ? 1. : -1.;
        
        while(ximudisc <= taue) {
            rdisc = solver_inte(ximudisc, r, kr, krsign, nrr2, nrr3);
            if (rdisc >= rin) {
                status = 2;
                //std::cout << "krsign = " << krsign << std::endl;
                break;
            }
            ximudisc += 2. * taumu0;
            kmusigndisc *= -1.;
        }
    }
    ismuobsset = true;
    switch (status) {
        case 0: 
            return -2.;
        case 1: {
            muobs = muinf;
            kmuobssign = kmusign;
            return muinf;
        }
        case 2: {
            re = rdisc;
            kmuobssign = kmusigndisc;
            krobssign = krsign;
            return rdisc;
        }
    }
}

/*! \brief The same with \ref calfate, but the ray-tracing part is done with \f$\textsc{sim5}\f$ package.
@param rin inner edge radius of the disc on the equatorial plane
@param pos photon four-position 
@param kmu photon wave vector
 */
double geoinf::calfate_sim5(const double rin, std::array<double, 4> & pos, std::array<double, 4> & kmu) {
    double precision = 1e-2, step = 0.01, rnow, munow, rdisc, muinf;
    sim5metric met;
    raytrace_data rtd;
    rnow = pos[1], munow = pos[2];
    raytrace_prepare(a, pos.data(), kmu.data(), NULL, precision, 0, &rtd);
    
    while(true) {
        //std::cout << "pos = " << pos << std::endl;
        //std::cout << "kmu = " << kmu << std::endl;
        raytrace1(pos.data(), kmu.data(), &step, &rtd);
        
        if (pos[1] <= 1.01 * rh) {
            status = 0;
            break;
        }
        
        if (pos[2] * munow <= 0. && rnow >= rin) {
            status = 2;
            rdisc = rnow;
            break;
        }
        
        if (pos[1] > 1e3) {
            status = 1;
            muinf = pos[2];
            break;
        }
        rnow = pos[1];
        munow = pos[2];
        //std::cout << "rnow = " << rnow << ", munow = " << munow << std::endl;
    }
    
    ismuobsset = true;
    switch (status) {
        case 0: 
            return -2.;
        case 1: {
            muobs = muinf;
            kmuobssign = dsign(kmu[2]);
            return muinf;
        }
        case 2: {
            re = rdisc;
            kmuobssign = dsign(kmu[2]);
            krobssign = dsign(kmu[1]);
            return rdisc;
        }
    }
    
}
            

/*
//! correct kmu; paying special attention to the case where the photon is not in the allowed region (in practice this is not possible
//! if the position is solved by evaluating the integration, but this may results from step-by-step raytracing routine (e.g., raytrace() inside SIM5)).
void geoinf::correctkmu(const double re, const double mue, const double kre, const double kmue, const int nrr2, const int nrr3,
    const int nmumin, const int nmumax, std::array<double, 4> & pos, std::array<double, 4> & kmu) {
    double krsign, kthetasign;
    double taue, taur2, taur3, ximu, xie;
    double disrr2, disrr3;
    krsign = (kmu[1] > 0. ? 1. : -1.);
    kthetasign = (kmu[2] > 0. ? 1. : -1.);
    int n1, n2;

    if (mue == 0.) {
        xie = taumu0;
    } else {
        xie = calximu(mue);
    }

    if (!isallowed(pos[1], pos[2])) {
        if (q >= 0.) {
        // mu turning point
            if (pos[2] > muplus) {
                pos[2] = muplus;
                // try to find the r-coordinate of mu turning point
                // note that taumu0 is always calculated
                pos[1] = solver_inte(xie, re, kre, krsign, n1, n2);
            } else {
                // r turning point
                // first find the turning point
                disrr2 = std::abs(pos[1] - rr2);
                disrr3 = std::abs(pos[1] - rr3);

                // if the turning point is r1
                if (re >= rr1) {
                    pos[1] = rr1;
                    taue = caltaur_r1(re);

                } else {
                    taur2 = caltaur_r2(re);
                    taur3 = taur2r3 - taur2;

                    if (disrr2 < disrr3) {
                        pos[1] = rr2;
                        taue = (kre >= 0.) ? taur2 + 2. * (double)(nrr3) * taur2r3 : taur3 + (double)(2. * nrr3 - 1) * taur2r3;
                    } else {
                        pos[1] = rr3;
                        taue = (kre >= 0.) ? caltaur_r2(re) + (double)(2. * nrr2 - 1) * taur2r3 : taur3 * (double)(2. * nrr2) * taur2r3;
                    }
                }

                pos[2] = solvemu_inte(taue, mue, kmue, kthetasign, n1, n2);
            }
        } else {
            // q < 0 case; now there are two turning points in mu but no turning point in r

            double ximue_mumax, ximue_mumin;
            ximue_mumax = calximu(mue);
            ximue_mumin = taumu0 - calximu(mue);
            taue = caltaur_inf_negq(re);

            if (pos[2] > muplus) {
                pos[2] = muplus;

                if (kmue >= 0.) {
                    ximu = ximue_mumin + (double)(2. * nmumin - 1) * taumu0;
                } else {
                    ximu = ximue_mumax + (double)(2. * nmumin) * taumu0;
                }
            } else {
                pos[2] = muminus;
                if (kmue >= 0.) {
                    ximu = ximue_mumin + (double)(2. * nmumax) * taumu0;
                } else {
                    ximu = ximue_mumax + (double)(2. * nmumax - 1) * taumu0;
                }
            }

            pos[1] = solver_inte(ximu, re, kre, krsign, n1, n2);
        }
    }

    photon_momentum(a, pos[1], pos[2], l, q, krsign, kthetasign, kmu.data());
}
*/

/*! \brief Inverse function of \ref calximu_all. Solves for \f$\mu\f$ given \f$\mu-\f$integral along the whole geodesic,
initial position and momentum. 

@param inte \f$\mu-\f$integral along the whole geodesic
@param mue \f$\mu\f$ at the starting point
@param kmue sign of \f$k^\theta\f$ at the starting point
@param kmufsign output; sign of \f$k^\theta\f$ at the end point
@param nmuminus output; number of encounting the lower turning point
@param nmuplus output; number of encounting the upper turning point
*/
double geoinf::solvemu_inte(const double inte, const double mue, const double kmue, double & kmufsign, int & nmuminus, int & nmuplus) {
    double xie, muf;
    
    bool lower = false;
    if (q < 0.) {
        lower = (mue < 0.);
    } else if (q == 0.) {
        lower = (mue < 0.) || ((mue == 0.) && (kmue > 0.));
    }
    
    if (lower) {
        muf = -1. * solvemu_inte(inte, -1. * mue, -1. * kmue, kmufsign, nmuminus, nmuplus);
        kmufsign *= -1.;
        return muf;
    }
    
    xie = calximu(mue);

    muf = solvemu_inte(inte, mue, kmue, xie, kmufsign, nmuminus, nmuplus);
    return muf;
}

/*! 
\brief Same with \ref solvemu_inte(const double, const double, const double, double &, int &, int &), but with pre-calculated
\f$\mu-\f$integral between \f$\mu_e\f$ and \f$\mu_+\f$.

Return: \f$\mu_f\f$, value of \f$\mu\f$ at the end point of the geodesic

@param inte \f$\mu-\f$integral along the whole geodesic
@param mue \f$\mu\f$ at the starting point
@param kmue sign of \f$k^\theta\f$ at the starting point
@param xie \f$\mu-\f$integral between \f$\mu_e\f$ and \f$\mu_+\f$; see \ref calximu.
@param kmufsign output; sign of \f$k^\theta\f$ at the end point
@param nmuminus output; number of encounting the lower turning point
@param nmuplus output; number of encounting the upper turning point
*/
double geoinf::solvemu_inte(const double inte, const double mue, const double kmue, 
    const double xie, double & kmufsign, int & nmuminus, int & nmuplus) {
    double inte1, muf, ximuf, ximuf_minus, xiall;
    int nbounce;
    xiall = (q > 0. ? 2. * taumu0 : taumu0);
    
    bool lower = false;
    if (q < 0.) {
        lower = (mue < 0.);
    } else if (q == 0.) {
        lower = (mue < 0.) || ((mue == 0.) && (kmue > 0.));
    }

    if (lower) {
        muf = -1. * solvemu_inte(inte, -1. * mue, -1. * kmue, xie, kmufsign, nmuminus, nmuplus);
        kmufsign *= -1.;
        return muf;
    }

    double ximue_mumin;
    ximue_mumin = xiall - xie;
    nbounce = (int)(std::floor(inte / 2. / xiall));
    inte1 = inte - 2. * (double)(nbounce) * xiall;
    nmuplus = nbounce;
    nmuminus = nbounce;

    if (kmue <= 0.) {
        //ximuf = (nmumin == nmumax) ? ximue_mumax - inte + (double)(2 * nmumin) * taumu0 : inte - ximue_mumax - (double)(2 * nmumin) * taumu0;
        if (inte1 < xie) {
            ximuf = xie - inte1;
            kmufsign = -1.;
        } else if (inte1 < xie + xiall) {
            ximuf = inte1 - xie;
            kmufsign = 1.;
            nmuplus += 1;
        } else {
            ximuf = xiall - (inte1 - xie - xiall);
            kmufsign = -1.;
            nmuplus += 1;
            nmuminus += 1;
        }
    } else {
        //ximuf_minus = (nmumin == nmumax) ? ximue_mumin - inte + (double)(2 * nmumax) * taumu0 : inte - ximue_mumin - (double)(2 * nmumax) * taumu0;
        if (inte1 < ximue_mumin) {
            ximuf_minus = ximue_mumin - inte1;
            kmufsign = 1.;
        } else if (inte1 < ximue_mumin + xiall) {
            ximuf_minus = inte1 - ximue_mumin;
            kmufsign = -1.;
            nmuminus += 1;
        } else {
            ximuf_minus = xiall - (inte1 - ximue_mumin - xiall);
            kmufsign = 1.;
            nmuminus += 1;
            nmuplus += 1;
        }
        ximuf = xiall - ximuf_minus;
    }

    muf = solvemu(ximuf);
    return muf;
}

/*! 
\brief Inverse function of \ref caltaur_all. Solve for \f$r\f$ at the end of the geodesic given \f$r-\f$ integral along the whole 
geodesic, initial position and momentum.

Return: \f$r_f\f$, \f$r\f$ at the end of the geodesic

@param inte \f$r-\f$integral along the whole geodesic
@param re \f$r\f$ at the starting point
@param kre sign of \f$k^r\f$ at the starting point
@param krfsign output; sign of \f$k^r\f$ at the end
@param nrr2 output; number of encounting the lower turning point
@param nrr3 output; number of encounting the upper turning point
*/
double geoinf::solver_inte(const double inte, const double re, const double kre, double & krfsign, int & nrr2, int & nrr3) {
    double taue = caltaue(re);
    double rf = solver_inte(inte, re, kre, taue, krfsign, nrr2, nrr3);
    return rf;
}

/*! 
\brief Same with \ref solver_inte(const double, const double, const double, double &, int &, int &), but
with pre-calculated \f$r-\f$ integral between \f$r_e\f$ and the corresponding turning point (see \ref caltaue).

Return: \f$r_f\f$, \f$r\f$ at the end of the geodesic

@param inte \f$r-\f$integral along the whole geodesic
@param re \f$r\f$ at the starting point
@param kre sign of \f$k^r\f$ at the starting point
@param taue \f$r-\f$ integral between \f$r_e\f$ and the corresponding turning point (see \ref caltaue).
@param krfsign output; sign of \f$k^r\f$ at the end
@param nrr2 output; number of encounting the lower turning point
@param nrr3 output; number of encounting the upper turning point
 */
double geoinf::solver_inte(const double inte, const double re, const double kre, const double taue, double & krfsign, int & nrr2, int & nrr3) {
    double inte1, taue3, tauf, tauf3;
    double rf = 0.;
    double nbounce;

    if (q >= 0.) {
        if (locatere(re)) {
            if (rr1 > rh) {
                tauf = (kre >= 0.) ? taue + inte : std::abs(taue - inte);
                rf = solver_r1(tauf);
                krfsign = (kre < 0. && inte < taue) ? -1. : 1.;
                nrr3 = (kre < 0. && inte > taue) ? 1 : 0;
                nrr2 = 0;
            } else {
                tauf = (kre >= 0.) ? taue + inte : taue - inte;
                rf = solver_r1(tauf);
                nrr2 = 0;
                nrr3 = 0;
                krfsign = kre >= 0. ? 1. : -1.;
            }
                
        } else {
            taue3 = taur2r3 - taue;
            nbounce = (double)(std::floor(inte / 2. / taur2r3));
            inte1 = inte - nbounce * taur2r3;
            nrr2 = nbounce;
            nrr3 = nbounce;
            if (rr3 > rh) {
                if (kre >= 0.) {
                    //tauf = (nrr2 == nrr3) ? taue - inte + (double)(2 * nrr3) * taur2r3 : inte - taue - (double)(2 * nrr3) * taur2r3;
                    if (inte1 < taue) {
                        tauf = taue - inte1;
                        krfsign = 1.;
                    } else if (inte1 < taue + taur2r3) {
                        tauf = inte1 - taue;
                        krfsign = -1.;
                        nrr2 += 1;
                    } else {
                        tauf = taur2r3 - (inte1 - taue - taur2r3);
                        krfsign = 1.;
                        nrr2 += 1;
                        nrr3 += 1;
                    }
                } else {
                    //tauf3 = (nrr2 == nrr3) ? taue - inte + (double)(2 * nrr2) * taur2r3 : inte - taue3 -(double)(2 * nrr2) * taur2r3;
                    if (inte1 < taue3) {
                        tauf3 = taue3 - inte1;
                        krfsign = -1.;
                    } else if (inte1 < taue3 + taur2r3) {
                        tauf3 = inte1 - taue3;
                        krfsign = 1.;
                        nrr3 += 1;
                    } else {
                        tauf3 = taur2r3 - (inte1 - taue3 - taur2r3);
                        krfsign = -1.;
                        nrr3 += 1;
                        nrr2 += 1;
                    }
                    tauf = taur2r3 - tauf3;
                }
            } else {
                tauf = (kre <= 0.) ? taue + inte : std::abs(taue - inte);
                krfsign = (kre > 0. && inte < taue) ? 1. : -1.;
                nrr2 = (kre > 0 && inte > taue) ? 1 : 0;
            }

            rf = solver_r2(tauf);
        }

    } else {
        // q < 0 case; now there are two turning points in mu but no turning point in r
        nrr2 = 0;
        nrr3 = 0;
        tauf = (kre > 0.) ? taue - inte : taue + inte;

        rf = solver_inf_negq(tauf);
        krfsign = (kre > 0.) ? 1. : -1.;
    }

    return rf;
}

/*!
\brief Obsolete; solve the radius for self-irradiation. Now this can be done with \ref calfate.
*/
double geoinf::solver_selfirr(const double re, const double kr, double & krfsign) {
    // Photons with negative q cannot reach the disc
    if (q < 0.)
        throw("geoinf::solver_selfirr(): q < 0");
    double intr, rf;
    int nrr2, nrr3;
    intr = 2. * taumu0;
    rf = solver_inte(intr, re, kr, krfsign, nrr2, nrr3);
    if (rf <= rh)
        rf = -1.;

    if (locatere(re)) {
        if (nrr3 > 0 && rr1 <= rh)
            rf = -1.;
    } else {
        if (rr2 <= rh || (rr3 <= rh && nrr3 > 0))
            rf = -1.;
    }

    return rf;
}

/*!
\brief Obsolete. Just an older version of \ref findxtrid.
*/
int geoinf::findxtrid_old(const tridgeo & trid, const double r, const double mu, const double kr,
    const double kmu, double & rnow, double & munow, std::array<double, 4> & kmunow) {

    // we will choose step size adaptively
    //double dinte1 = 1e-2, dinte2 = 1e-4;
    double inte = 0;
    int status;
    double krsign, kmusign;
    double taue, xie;
    int nrr2, nrr3, nmuplus = 0, nmuminus = 0;
    double rfunc, rfunc0, rfunc2, dr, rmin;
    std::array<double, 5> rfunc_coeff = {1., 0., a2 - l * l - q, 2. * (q + (l -a) * (l - a)), -1. * a2 * q};

    taue = caltaue(r);
    xie = calximu(mu);

    rnow = r;
    munow = mu;
    krsign = kr;
    kmusign = kmu;

    //std::cout << "kmu = " << kmu << ", kr = " << kr << std::endl;

    while(true) {
        //std::cout << "inte = " << inte << ", rnow = " << rnow << ", munow = " << munow <<
        // ", nmuminus = " << nmuminus << ", xie = " << xie << std::endl;
        //std::cout << "inte = " << inte << ", rnow = " << rnow << ", munow = " << munow <<
        //  ", nmuminus = " << nmuminus << ", xie = " << xie << std::endl;
        //std::cout << "mumin = " << trid.mumin << ", mumax = " << trid.mumax << std::endl;

        // munow < 0.: below equatorial plane
        // nmuninus > 0: above the equatorial plane, but has returned from the turning point below the equatorial plane
        if (munow < -1e-8 || (q > 0. && nmuminus > 0)) {
            status = 2;
         //   std::cout << "Return radiation:";
         //   std::cout << "\n";
            break;
        }

        // hit
        if (trid.inside(rnow, munow)) {
            status = 1;
            photon_momentum(a, rnow, munow, l, q, krsign, kmusign, kmunow.data());
            break;
        }

        // not possible
        if (!ispossible(trid, rnow, munow, krsign, kmusign)) {
         //   std::cout << "Not possible: ";
         //   std::cout << "\n";
            status = 0;
            break;
        }

        if (rnow < 2.0) {
            dr = 0.01;
        } else {
            if (rnow < trid.rmin - 0.1) {
                dr = trid.rmin - 0.1 - rnow;
            } else if (rnow > trid.rmax + 0.1) {
                dr = rnow - trid.rmax - 0.1;
            } else {
                dr = 0.01;
            }
        }

        rmin = (rnow - dr < 1.1 * rh) ? 1.1 * rh : rnow - dr;
        if (rnow < 2.0)
            rmin = rnow - dr;
        std::cout << "rmin = " << rmin << std::endl;

        rfunc = polyval(rfunc_coeff, rnow);
        std::cout << "rfunc = " << rfunc << std::endl;
        rfunc0 = polyval(rfunc_coeff, rnow + dr);
        std::cout << "rfunc = " << rfunc0 << std::endl;
        rfunc2 = polyval(rfunc_coeff, rmin);
        std::cout << "rfunc = " << rfunc2 << std::endl;
        if (rfunc0 > rfunc) rfunc = rfunc0;
        if (rfunc2 > rfunc) rfunc = rfunc2;


        inte += dr * std::sqrt(rfunc);
        //std::cout << "dr = " << dr << ", inte = " << inte << std::endl;

        try {
            rnow = solver_inte(inte, r, kr, taue, krsign, nrr2, nrr3);
        }
        catch(const char * error) {
            std::cout << error << std::endl;
            std::cout << "l = " << l << std::endl;
            std::cout << "q = " << q << std::endl;
            std::cout << "rnow = " << rnow << std::endl;
            std::cout << "nrr = " << nrr << std::endl;
            std::cout << "rr1 = " << rr1 << ", rr2 = " << rr2 << ", rr3 = " << rr3 << std::endl;
            return 2;
        }

        if (rnow <= 1.01 * rh) {
            status = 3;
            std::cout << "BH:";
            //std::cout << "\n";
            break;
        }

        munow = solvemu_inte(inte, mu, kmu, xie, kmusign, nmuminus, nmuplus);

    }

    return status;
}

/*!
\brief Find the position where the photon hits the corona.

Return: an integer indicating the status.
    - 0: it is impossible for the photon to reach the corona
    - 1: the photon enters the corona
    - 2: the photon strikes at the disc on the equatorial plane
    - 3: the photon enters BH event horizon
    
@param trid a \ref tridgeo object for different kinds of corona
@param r starting \f$r\f$
@param mu starting \f$\mu\f$
@param kr sign of \f$k^r\f$
@param kmu sign of \f$k^\theta\f$
@param rin inner edge radius of the disc on the equatorial plane
@param rnow output; \f$r\f$ where the photon enters the corona
@param munow output; \f$\mu\f$ where the photon enters the corona
@param kmunow output; photon wave vector while entering the corona
*/
int geoinf::findxtrid(const tridgeo & trid, const double r, const double mu, const double kr,
    const double kmu, const double rin, double & rnow, double & munow, std::array<double, 4> & kmunow) {

    // we will choose step size adaptively
    //double dinte1 = 1e-2, dinte2 = 1e-4;
    double inte, dinte, dr = 1e-3, muprev;
    int status;
    double krsign, krsigndisc, kmusign;
    double taue, xie;
    int nrr2, nrr3, nmuplus = 0, nmuminus = 0, nmuminusdisc, nmuplusdisc;
    double rfunc, rfunc0, rfunc2;
    std::array<double, 5> rfunc_coeff = {1., 0., a2 - l * l - q, 2. * (q + (l -a) * (l - a)), -1. * a2 * q};
    double rstart, krstart, mustart, kmustart;
    double rdisc, intedisc;
    bool precal = false;

    taue = caltaue(r);
    xie = calximu(mu);

    rnow = r;
    munow = mu;
    muprev = mu + 1e-3 * kmu;
    krsign = kr;
    kmusign = kmu;

    // first check the possibility
    if (!ispossible(trid, rnow, munow, krsign, kmusign)) {
        //std::cout << "rnow = " << rnow << std::endl;
        status = 0;
        //std::cout << "now possible " << std::endl;
        return status;
    }

    // first step: propagate the photon close to the corona
    if (q >= 0.) {
        if ((nrr == 4 && locatere(rnow)) || nrr == 2) {
            if (krsign >= 0. && rnow < trid.rmin - 0.1) {
                inte = caltaur_r1(trid.rmin - 0.1) - caltaur_r1(rnow);
                krsign = 1.;
                munow = solvemu_inte(inte, mu, kmu, kmusign, nmuminus, nmuplus);
                rnow = trid.rmin - 0.1;
                precal = true;
            } else if (krsign <= 0. && rnow > trid.rmax + 0.1) {
                inte = caltaur_r1(rnow) - caltaur_r1(trid.rmax + 0.1);
                krsign = -1.;
                munow = solvemu_inte(inte, mu, kmu, kmusign, nmuminus, nmuplus);
                rnow = trid.rmax + 0.1;
                precal = true;
            } else {
                inte = 0.;
            }
                
            if (q > 0.) {
                // find out the radius when crossing the equatorial plane
                intedisc = ((mu >= 0. && kmu <= 0.) || (mu <= 0. && kmu >= 0.)) ?
                    taumu0 + calximu(std::abs(mu)) : taumu0 - calximu(std::abs(mu));
                while(intedisc < inte) {
                    rdisc = solver_inte(inte, r, kr, krsigndisc, nmuminusdisc, nmuplusdisc);
                    if (rdisc >= rin) {
                        status = 2;
                        return status;
                    }
                    intedisc += 2. * taumu0;
                }
            }
        }
            
    } else {
        if (krsign > 0. && rnow < trid.rmin - 0.1) {
            inte = caltaur_inf_negq(rnow) - caltaur_inf_negq(trid.rmin - 0.1);
            munow = solvemu_inte(inte, mu, kmu, kmusign, nmuminus, nmuplus);
            krsign = 1.;
            rnow = trid.rmin - 0.1;
        } else if (krsign < 0. && rnow > trid.rmax + 0.1) {
            inte = caltaur_inf_negq(trid.rmax + 0.1) - caltaur_inf_negq(rnow);
            krsign = -1.;
            munow = solvemu_inte(inte, mu, kmu, kmusign, nmuminus, nmuplus);
            rnow = trid.rmax + 0.1;
        }
    }

    // then start
    // If the photon crossed the disc during the steps above, we need to find out 
    // the radius
    /*
    if (munow < -1e-8 || (q > 0. && nmuminus > 0)) {
        status = 2;
        return status;
    }
    */

    // hit
    if (trid.inside(rnow, munow)) {
        status = 1;
        photon_momentum(a, rnow, munow, l, q, krsign, kmusign, kmunow.data());
        return status;
    }

    // not possible
    if (!ispossible(trid, rnow, munow, krsign, kmusign)) {
        status = 0;
        return status;
    }
    rstart = rnow;
    mustart = munow;
    krstart = krsign;
    kmustart = kmusign;

    inte = 0.;
    long count = 0;
    while (true) {
        rfunc = polyval(rfunc_coeff, rnow);
        rfunc0 = polyval(rfunc_coeff, rnow + dr);
        rfunc2 = polyval(rfunc_coeff, rnow - dr);
        if (rfunc0 > rfunc) rfunc = rfunc0;
        if (rfunc2 > rfunc) rfunc = rfunc2;
        
        dinte = dr / std::sqrt(rfunc);
        if (dinte > 1e-3)
            dinte = 1e-3;

        inte += dinte;
        muprev = (count == 0 && mu == 0 && !precal) ? kmu * -1e-3 : munow;
        try {
            rnow = solver_inte(inte, rstart, krstart, taue, krsign, nrr2, nrr3);
        }
        catch(const char * error) {
            std::cout << error << std::endl;
            std::cout << "l = " << l << std::endl;
            std::cout << "q = " << q << std::endl;
            std::cout << "r = " << r << std::endl;
            std::cout << "rstart = " << rstart << ", kr = " << kr << std::endl;
            std::cout << "nrr = " << nrr << std::endl;
            std::cout << "rr1 = " << rr1 << ", rr2 = " << rr2 << ", rr3 = " << rr3 << std::endl;
            std::cout << "trid.rmin = " << trid.rmin << std::endl;
            std::cout << "ispossible = " << ispossible(trid, r, mu, kr, kmu) << std::endl;
            std::cout << "int_r2r3 = " << taur2r3 << std::endl;
            std::cout << "dr * std::sqrt(rfunc) = " << dr / std::sqrt(rfunc) << std::endl << std::endl;
            return 2;
        }

        if (rnow <= 1.01 * rh) {
            status = 3;
            break;
        }

        munow = solvemu_inte(inte, mustart, kmustart, xie, kmusign, nmuminus, nmuplus);
        //std::cout << "rnow = " << rnow << ", munow = " << munow << std::endl;

        if (munow * muprev <= 0. && rnow >= rin) {
            status = 2;
            break;
        }

        // hit
        if (trid.inside(rnow, munow)) {
            status = 1;
            photon_momentum(a, rnow, munow, l, q, krsign, kmusign, kmunow.data());
            break;
        }

        // not possible
        if (!ispossible(trid, rnow, munow, krsign, kmusign)) {
            status = 0;
            break;
        }
        count++;
    }

    //std::cout << "status = " << status << std::endl;
    return status;
}

/*!
\brief Same with \ref findxtrid, but the ray-tracing is done with \f$\textsc{sim5}\f$ package.

@param rin the inner edge of the disc
@param trid the \ref tridgeo object
@param pos initial photon four position
@param kmu initial photon wave vector
 */
int geoinf::findxtrid_sim5(const double rin, const tridgeo & trid, std::array<double, 4> & pos, std::array<double, 4> & kmu) {
    double precision = 1e-2, step = 0.01, kmusign, muprev;
    kmusign = (kmu[2] >= 0.) ? -1. : 1.;
    long count = 0;
    raytrace_data rtd;
    raytrace_prepare(a, pos.data(), kmu.data(), NULL, precision, 0, &rtd);

    while(true) {
        //std::cout << "pos = " << pos << std::endl;
        //std::cout << "kmu = " << kmu << std::endl;
        muprev = (count == 0) ? pos[2] + kmusign * 1e-3 : pos[2];
        raytrace1(pos.data(), kmu.data(), &step, &rtd);
        if (pos[1] <= 1.01 * rh) {
            status = 3;
            break;
        }
        
        if (muprev * pos[2] <= 0. && pos[1] >= rin) {
            status = 2;
            break;
        }

        if (trid.inside(pos[1], pos[2])) {
            status = 1;
            break;
        }
        
        if (!ispossible(trid, pos[1], pos[2], kmu[1], kmu[2])) {
            status = 0;
            break;
        }
        ++count;
    }
    return status;
}

/*!
Overloaded << operator; printing object info
 */
std::ostream & operator << (std::ostream & cout, const geoinf & geo) {
    cout << "======================================" << std::endl;
    cout << "Physical pamameters:" << std::endl;
    cout << "a = " << geo.a << ", muobs = " << geo.muobs << ", alpha = " << geo.alpha
        <<", beta = " << geo.beta << std::endl;
    cout << "--------------------------------------" << std::endl;
    cout << "Integrations:" << std::endl;
    cout << "taumuinf = " << geo.taumuinf << ", taumu0 = " << geo.taumu0 << ", taurinf = " << geo.taurinf
        << std::endl;
    cout << "Re = " << geo.re <<", r1 = " << geo.rr1 << ", muplus = " << geo.muplus << std::endl;
    /*
    if (geo.areturnsfound()){
        cout << "--------------------------------------" << std::endl;
        cout << "Turning points:" << std::endl;
        cout << geo.nmuturn << " turning points in mu";
        if (geo.nmuturn > 0)
            std::cout << ": [" << geo.rturn << ", " << geo.muplus << "]" << std::endl;
        else
            cout << "." << std::endl;
        cout << geo.nrturn << " turning points in r";
        if (geo.nrturn > 0)
            cout << ": [" << geo.rr1 << ", " << geo.muturn<< "]" << std::endl;
        else
            cout << "." << std::endl;
    }
    */
    cout << "======================================" << std::endl;
    return cout;
}

/*
//! To be documentated
//! @param a spin
void compare_raytrace(const double a, const std::array<double, 4> & pos0, const std::array<double, 4> & kmu0) {
    double step;
    raytrace_data rtd;
    std::vector<double> rarrray, muarrray, rarrint, muarrint;
    std::array<double, 4> pos = pos0, kmu = kmu0;

    geoinf geo(a, pos[1], pos[2], kmu);

    raytrace_prepare(a, pos.data(), kmu.data(), NULL, 0.01, 0, &rtd);
    while(pos[2] >= 0.) {
        step = 1e-3;
        raytrace(pos.data(), kmu.data(), NULL, &step, &rtd);
        rarrray.push_back(pos[1]);
        muarrray.push_back(pos[2]);
        //std::cout << "pos = " << pos << std::endl;
        //std::cout << "kmu = " << kmu << std::endl;
        if (pos[1] > 1e4) {
            std::cout << "geo.l = " << geo.l << ", geo.muplus = " << geo.muplus << ", geo.rr1 = " 
                << geo.rr1 << ", geo.rr2 = " << geo.rr2 << ", geo.rr3 = " << geo.rr3 << std::endl;
            std::cout << "geo.q = " << geo.q << ", geo.nrr = " << geo.nrr << ", geo.rr1 = " << geo.rr1 
                << ", geo.rr2 = " << geo.rr2 << ", geo.rr3 = " << geo.rr3 << std::endl;
            std::cout << "pos[1] = " << pos[1] << std::endl;
        }
        if (pos[1] > 1e2 || pos[1] != pos[1]) {
            break;
        }
    }

    double inte, dinte = 1e-4, rnow, munow;
    int nmuminus, nmuplus, nrr2, nrr3;
    double mufsign, rfsign;

    for (long i = 0; i < 10000; ++i) {
        inte = (double)(i) * dinte;
        try {
            rnow = geo.solver_inte(inte, pos0[1], kmu0[1], rfsign, nrr2, nrr3);
            munow = geo.solvemu_inte(inte, 0., -1., mufsign, nmuminus, nmuplus);
        }
        catch(const char * error) {
            std::cout << error << std::endl;
        }
        rarrint.push_back(rnow);
        muarrint.push_back(munow);
        if (munow < 0.) break;
    }
    wdoublevec(rarrray, "rarr.dat");
    wdoublevec(muarrray, "muarr.dat");
    wdoublevec(rarrint, "rarrint.dat");
    wdoublevec(muarrint, "muarrint.dat");
}
*/
