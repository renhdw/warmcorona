/*! \file calhotcross.cpp 
 * This file contains several functions for calculating the ``hot`` Compton cross section, i.e., the scattering 
cross section of hot electrons with Maxwell-Jutter distribution:
\f{equation}
\sigma_{h} = \frac{1}{n_e} \int \frac{dn_e}{d\gamma} (1-\mu \beta) \sigma_{\rm KN}(x^\prime),
\f}
where \f$n_e\f$ is electron number density, \f$dn_e/d\gamma\f$ is electron velocity distribution, \f$beta\f$ is electron dimensionless velocity,
\f$\gamma\f$ is electron Lorentz factor, \f$\sigma_{\rm KN}\f$ is the total Klein-Nishina cross section, \f$x^\prime \equiv \gamma x (1-\mu \beta)\f$,
and \f$x\f$ is photon dimensionless energy. The maxwell Juttern distribution:
\f{equation}
\frac{dn_e}{d\gamma} = \frac{n_e\gamma^2 \beta}{\theta_T K_2(\theta_T)} e^{-\gamma/\theta_T},
\f}
where \f$\theta_T\f$ is electron dimensionless temperature, and \f$K_2\f$ is the modified Bessel function of the second kind. Wring the hot cross section
explicitely:
\f{equation}
\sigma_{h} = \frac{1}{2\theta_T K_2(\theta_T)} \int_0^1 d\beta \int_{-1}^1 d\mu\ \gamma^5 
\beta^2 e^{-\gamma/\theta_T}(1-\mu \beta) \sigma_{\rm KN}(x^\prime).
\f}

We pre-calculate the hot cross section for a grid of \f$\theta_T\f$ and \f$x\f$, and write the data into a binary file.
This file also contains \ref readhotx, an interface to read the binary file and perform interpolation.

*/
#include <iostream>
#include <cmath>
#include <string>
#include <experimental/filesystem>
#include "scatter.h"
#include "utils.h"
#include "const.h"
#include "calhotcross.h"

//const long NGAUSS1 = 11;
//const long NBIN = 10000;
//! Maximum order in trapezoid method
const int NMAX = 23;
//double hotcross_calmuinte(const double beta, const double x);
double hotcross_integrand_mu(const double mu, const double beta, const double x);
double hotcross_integrand_beta(const double thetat, const double beta);
double hotcross_calbetainte(const double thetat, const double x);
double int_trape1(double (*pt2func)(const double, void *), const double a, const double b, const double acc, void * params);
void int_trapezoid(double(*pt2func)(const double, void *), const double a, const double b, const int n, void * params, double & res);
double hotcross_integrand_beta1(const double beta, void * params);
double hotcross_integrand_mu1(const double mu, void * params);
//double x2_integrand(const double x, void * param);

namespace fs = std::experimental::filesystem;
/*
void init_hotx();
void test();

int main() {
    test();
}
*/
//! Global vector used in reading data from disc; 
std::vector<double> HOTX::LOGXARR;
//! Global vector used in reading data from disc; 
std::vector<double> HOTX::LOGTHETATARR;
//! Global vector used in reading data from disc; 
std::vector<double> HOTX::HOTXARR;

/*! \brief Write data to disk.
 */
void writehotx() {
    std::vector<double> logthetat(HOTX::NTHETAT), logx(HOTX::NX), hotcross(HOTX::NX * HOTX::NTHETAT);
    
    for (long i = 0; i < HOTX::NTHETAT; ++i) {
        logthetat[i] = HOTX::LOGTHETATMIN + (double)(i) * HOTX::DLOGTHETAT;
    }
    for (long i = 0; i < HOTX::NX; ++i) {
        logx[i] = HOTX::LOGXMIN + (double)(i) * HOTX::DLOGX;
    }
    
    for (long ix = 0; ix < HOTX::NX; ++ix) {
        for (long it = 0; it < HOTX::NTHETAT; ++it) {
            hotcross[ix * HOTX::NTHETAT + it] = hotcross_calbetainte(std::pow(10., logthetat[it]), std::pow(10., logx[ix]));
            showprogress((double)(ix * HOTX::NTHETAT + it) / (double)(hotcross.size()));
        }
    }
    wdoublevec(logthetat, "logthetat.dat");
    wdoublevec(logx, "logx.dat");
    wdoublevec(hotcross, "hotx.dat");
}

/*! \brief Read data into memory.
 */
void init_hotx() {
    init_hotx_local(HOTX::LOGTHETATARR, HOTX::LOGXARR, HOTX::HOTXARR);
}

void init_hotx_local(std::vector<double> & logthetatarr, std::vector<double> & logxarr, std::vector<double> & hotxarr) {
    std::string tpath = HOTX::HOTXDIR + "logthetat.dat", xpath = HOTX::HOTXDIR + "logx.dat", hotpath = HOTX::HOTXDIR + "hotx.dat";
    if (!fs::exists(tpath)) {
		std::cout << "tpath = " << tpath << std::endl;
        throw(tpath + " does not exist!");
        return;
    }
    if (!fs::exists(xpath)) {
        throw(xpath + " does not exist!");
        return;
    }
    if (!fs::exists(hotpath)) {
        throw(hotpath + " does not exist!");
        return;
    }
    rdoublevec(logthetatarr, tpath);
    rdoublevec(logxarr, xpath);
    rdoublevec(hotxarr, hotpath);
}

double readhotx(const double thetat, const double x) {
    if (HOTX::LOGTHETATARR.size() == 0 || HOTX::LOGXARR.size() == 0 || HOTX::HOTXARR.size() == 0) {
        std::cerr << "Hotcross section not read!" << std::endl;
        throw("Hotcross not read");
    }
    return readhotx_local(thetat, x, HOTX::LOGTHETATARR, HOTX::LOGXARR, HOTX::HOTXARR);
}

/*
double readhotx_mpi(const double thetat, const double x) {
    return readhotx_local(thetat, x, HOTX::LOGTHETATARR, HOTX::LOGXARR, HOTX::HOTXARR);
}
*/
    
/*! \brief Given \f$\theta_T\f$ and \f$x\f$, calculate the hot cross section by bilinear interpolation.

For cold plasma and low energy photon, we just return Thompson scattering cross section. For cold plasma, we return
total Klein Nishina cross section.
@param thetat dimensionless electron temperature
@param x dimensionless photon energy
 */
double readhotx_local(const double thetat, const double x, const std::vector<double> & logthetatarr, 
    const std::vector<double> & logxarr, const std::vector<double> & hotxarr) {
        
    double logx, logthetat, hotx, hotx_t0, hotx_t1;
    long it0, it1, ix0, ix1;
    long index00, index01, index10, index11;
    
    if (thetat * x < 1.e-6)
        return 1.;
    
    if (thetat < -4.)
        return KN_xsection(x);
    
    if ((x > HOTX::XMIN && x < HOTX::XMAX) && (thetat > HOTX::THETATMIN && thetat < HOTX::THETATMAX)) {
        logx = std::log10(x);
        logthetat = std::log10(thetat);
        ix0 = (long)(std::floor((logx - std::log10(HOTX::XMIN)) / HOTX::DLOGX));
        ix1 = ix0 + 1;
        it0 = (long)(std::floor((logthetat - std::log10(HOTX::THETATMIN)) / HOTX::DLOGTHETAT));
        it1 = it0 + 1;
        //std::cout << ix0 << " " << ix1 << " " << it0 << " " << it1 << std::endl;
        
        index00 = ix0 * HOTX::NTHETAT + it0;
        index01 = ix0 * HOTX::NTHETAT + it1;
        index10 = ix1 * HOTX::NTHETAT + it0;
        index11 = ix1 * HOTX::NTHETAT + it1;
        
        hotx_t0 = (logxarr[ix1] - logx) * hotxarr[index00] + (logx - logxarr[ix0]) * hotxarr[index10];
        hotx_t1 = (logxarr[ix1] - logx) * hotxarr[index01] + (logx - logxarr[ix0]) * hotxarr[index11];
        hotx = ((logthetatarr[it1] - logthetat) * hotx_t0 + (logthetat - logthetatarr[it0]) * hotx_t1) / HOTX::DLOGTHETAT / HOTX::DLOGX;
    } else {
        hotx = hotcross_calbetainte(thetat, x);
    }
    
    return hotx;
}

//! A test function
void test() {
    //void *param;
    
    //double thetat;
    double te = 1000.;
    std::vector<double> ephoton = {1., 10., 100., 500., 1000., 5000.};
    //std::vector<double> ephoton = {std::pow(10., 2.) * ME_KEV};
    double hotx, x0;
    std::vector<double> thetatarr = {te / ME_KEV};
    //std::vector<double> thetatarr = {-2., -1., 0., 1., 2., 3.};
    std::cout << "log(w) = " << std::log10(ephoton[0] / ME_KEV) << std::endl;
    
    for (unsigned int i = 0; i < ephoton.size(); ++i) {
        for (unsigned int j = 0; j < thetatarr.size(); ++j) {
            x0 = ephoton[i] / ME_KEV;
            //thetat = std::pow(10., thetatarr[j]);
            //thetat = thetatarr[j];
            hotx = readhotx(thetatarr[j], x0) * 0.665245873e-24;
            std::cout << "hotcross = " << hotx / 1e-24 << " barns" << std::endl;
            //std::cout << "log hotcross / (cm^{-2}) = " << std::log10(hotx) << std::endl;
        }
    }
}

/*! \brief Evaluating \f$\beta\f$ integration. We first evaluate \f$\mu\f$ integral given \f$beta\f$.
@param thetat dimensionless electron temperature
@param x dimensionless photon energy
 */
double hotcross_calbetainte(const double thetat, const double x) {
    //double abissca[NGAUSS1], weights[NGAUSS1];
    //double integauss = 0.;
    //gauleg(0., 1., abissca, weights, NGAUSS1);
    //
    //for (int i = 0; i < NGAUSS1; ++i) {
    //    integauss += weights[i] * hotcross_integrand_beta(thetat, abissca[i]) * hotcross_calmuinte(abissca[i], x);
    //}
    //return 0.5 * integauss / thetat / bessel_k2(1. / thetat);
    std::array<double, 2> paramarr = {thetat, x};
    void * param = reinterpret_cast<void *>(&paramarr);
    double inte = int_trape1(&hotcross_integrand_beta1, 0., 1., 1e-3, param);
    
    return 0.5 * inte / thetat / bessel_k2(1. / thetat);
}
    
/*! \brief The \f$\beta\f$ integrand
@param beta the dimensionless electron velocity
@param params a reinterpreted pointer to std:array<double, 2> which contains \f$[\theta_T, x]\f$
 */
double hotcross_integrand_beta1(const double beta, void * params) {
    std::array<double, 2> * paramarr = reinterpret_cast<std::array<double, 2> *> (params);
    std::array<double, 2> muparams = {beta, (*paramarr)[1]};
    void * intparams = reinterpret_cast<void *>(&muparams);
    double muint = int_trape1(&hotcross_integrand_mu1, -1., 1., 1e-3, intparams);
    return hotcross_integrand_beta((*paramarr)[0], beta) * muint;
}

/*! \brief The part of the integrand which is function of \f$\beta\f$ only; i.e., \f$\gamma^5 \beta^2 e^{-\gamma/\theta_T}\f$
@param thetat dimensionless electron temperature
@param beta the dimensionless electron velocity
*/
double hotcross_integrand_beta(const double thetat, const double beta) {
    double gamma = 1. / std::sqrt(1. - beta * beta);
    double gamma2 = gamma * gamma;
    return gamma2 * gamma2 * gamma * beta * beta * std::exp(-1. * gamma / thetat);
}

/*! \brief \f$\mu\f$ integrand: \f$(1-\mu\beta)\sigma_{\rm KN}(x^\prime)\f$.
@param mu cosine of angle made by photon and electron
@param beta the dimensionless electron velocity
@param x dimensionless photon energy
*/
double hotcross_integrand_mu(const double mu, const double beta, const double x) {
    //std::cout << "beta = " << beta << std::endl;
    double xp, gamma;
    gamma = 1. / std::sqrt(1. - beta * beta);
    xp = gamma * (1. - mu * beta) * x;
    return (1. - mu * beta) * KN_xsection(xp);
}

/*! \brief the same with \ref hotcross_integrand_mu, but the arguments are passed with a void pointer,
to be compatible with \ref int_trape1.
@param mu cosine of angle made by photon and electron
@param params a reinterpreted pointer to std:array<double, 2> which contains \f$[\beta, x]\f$
 */
double hotcross_integrand_mu1(const double mu, void * params) {
    std::array<double, 2> * paramarr = reinterpret_cast<std::array<double, 2> *> (params);
    return hotcross_integrand_mu(mu, (*paramarr)[0], (*paramarr)[1]);
}

/*! \brief Trapezoid integral of a given order. Taken from `sim5` package; re-written to be able to pass additional arguments with a void pointer.
@param pt2func pointer to function to be integrated
@param a lower boundary of integral
@param b upper boundary of integral
@param n order
@param params void pointer to pass additional parameters
@param res result
 */
void int_trapezoid(double(*pt2func)(const double, void *), const double a, const double b, const int n, void * params, double & res) {
    double x, tnm, sum, del;
    int it, j;

    if(n==1){
       res = 0.5 * (b-a) * (pt2func(a, params) + pt2func(a, params));
    } else {
        for(it=1, j=1; j < n-1; j++)
            it <<= 1;
        tnm = (double) it;
        del = (b - a) / tnm;
        x = a + 0.5 * del;
        for(sum=0.0, j=1; j<=it; j++, x+=del) {
            sum += pt2func(x, params); 
        }
        res = 0.5 * (res + del * sum);
    }
}

/*! \brief Trapezoid integral. The function will return the result if 1) the given accurancy is reached; 2) the maximum
order is reached. Taken from `sim5` package; re-written to be able to pass additional arguments with a void pointer.
@param pt2func pointer to function to be integrated
@param a lower boundary of integral
@param b upper boundary of integral
@param acc the desired accurancy
@param params void pointer to pass additional parameters
 */
double int_trape1(double (*pt2func)(const double, void *), const double a, const double b, const double acc, void * params) {
    int n;
    double s = 0.0;
    double olds = -1000.;
    
    for (n = 0; n <= NMAX; ++n) {
        int_trapezoid(pt2func, a, b, n, params, s);
        
        if (n > 3) {
            if (std::abs(s - olds) < acc * std::abs(olds) || ((s == 0.) && (olds == 0.))) 
                break;
        }
        olds = s;
    }
    
    //std::cout << "n = " << n << std::endl;
    return s;
}

/*
double x2_integrand(const double x, void * param) {
    return x * x;
}
*/
