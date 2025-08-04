//! \file scatter.h
#ifndef _SLABSCATTER_H
#define _SLABSCATTER_H
#include <vector>
#include <array>
#include <random>
#include "detector.h"
#include "tridgeo.h"
#include "diskobj.h"
#include "superphoton.h"
#include "photon_dist.h"

class Ggammaparam {
public:
    double thetat, result;
};

class xsectphicdfparams {
public:
    double mu;
    double poldeg;
    double result;
};

class pol_phi_params {
public:
    double par1;
    double rndnum;
};

double KN_pdfx(const double x, const double xp);
void sample_omega(std::mt19937_64 & gen, double & theta, double & ephi);
double sample_maxjutt(const double te, std::mt19937_64 & gen);
double sample_maxjutt1(const double te, std::mt19937_64 & gen);
void sample_thompsonx(std::mt19937_64 & gen, const double poldeg, double & mu, double & phi);
//void sample_thompsonx_fixmu(std::mt19937_64 & gen, const double poldeg, const double mu, double & phi);
double sample_scatter_mu(const double beta, std::mt19937_64 & gen);
double sample_kn_phi(std::mt19937_64 & gen, const double poldeg, const double mu, const double x1_over_x);
double sample_thompson_phi(std::mt19937_64 & gen, const double poldeg, const double mu);
double loboost(const double nx, const double ny, const double nz, const double beta, const double gamma,
    const std::array<double, 4> & kmu_f, std::array<double, 4> & kmu_e);
double sample_KN_xp(const double x, std::mt19937_64 & gen);

void scatter(const bool pol, const bool recoil, const bool KN, const double en, const double poldeg, const std::array<double, 4> & pe,
    std::array<double, 4> & kmu_f, std::array<double, 4> & fmu_f, double & en_p, double & poldeg_p,
    std::array<double, 4> & kmu_f_p, std::array<double, 4> & fmu_f_p, std::mt19937_64 & gen);

void scatter_restframe(const bool KN, const double * pk, const double * pf, const double coskt,
    const double kphi, const double poldeg, double & poldeg_p, const double z, double * pk_e_p, double * pf_e_p);
void tridsca_sp(const bool pol, const double a, const double rin, const double mfp, const double te, const double weightnorm, 
    const double qfac0, const double ufac0, const double muinf_n0scatter, const double dr, superphoton sp, 
    const tridgeo & trid, std::vector<std::vector<double>> & vecarr_inf, std::vector<int> & nsca_inf, 
    std::vector<std::vector<double>> & vecarr_disc, std::vector<int> & nsca_disc, std::mt19937_64 & gen);
//void tridsca_sp_simple(const double a, const double mfp, const double te, const double dr,
//    superphoton sp, const tridgeo & trid, std::vector<std::ofstream> & ofiles);
void forward_rotation(const double * itx, const double * ity, const double * itz, const double * it, double * it1);
void backward_rotation(const double * itx, const double * ity, const double * itz, const double * it, double * it1);
//std::vector<double> receive_scatter_pack(const double tbb, const std::vector<double> & en,
//    const std::vector<double> & loboostarr, const std::vector<double> & scattermuarr, 
//    const std::vector<double> & inloboostarr, const std::vector<double> gfacarr);
void scatter_bias_photon_KN_pol(const bool pol, const bool recoil, const int stype, const double te, 
    const double dl, const double minweight, const double mfp, const double weightnorm,
    const tridgeo_sr & trid, superphoton_sr_pol & sp, std::mt19937_64 & gen, std::vector<std::ofstream> & ofiles);
void scatter_sr_slab(const bool pol, const bool recoil, const int stype, const double te, const double dl, const double minweight, const double mfp,
    const slab3d_sr & trid, superphoton_sr_pol & sp, std::vector<superphoton_sr_pol> & sps_escaped, std::mt19937_64 & gen);
double KN_xsection(const double x);
double klein_nishina(double a, double ap);
//double hc_klein_nishina(double we);
void sphere_KN(const bool progress, const bool recoil, const int type, const int pattern, const long npack, const double tau, const double te, const double dl, 
    const double weightnorm, const photon_dist & pphoton, std::vector<std::ofstream> & ofiles);
double bessel_k0(const double x);
double bessel_k1(const double x);
double bessel_k2(const double x);
void dummy_vec(const double * itz, double * itx);
void sample_electron_KN(const double te, const double x, std::mt19937_64 & gen, double & gamma, double & mu, double & phi);
void calpe(const double emu, const double ephi, const double * kmu, std::array<double, 4> & pe);

const std::array<double, 5> k0pi = {1.0, 2.346487949187396e-1, 1.187082088663404e-2, 2.150707366040937e-4, 1.425433617130587e-6};
const std::array<double, 3> k0qi = {9.847324170755358e-1, 1.518396076767770e-2, 8.362215678646257e-5};
const std::array<double, 5> k0p = {1.159315156584126e-1, 2.770731240515333e-1, 2.066458134619875e-2, 4.574734709978264e-4, 3.454715527986737e-6};
const std::array<double, 3> k0q = {9.836249671709183e-1, 1.627693622304549e-2, 9.809660603621949e-5};
const std::array<double, 8> k0pp = {1.253314137315499, 1.475731032429900e1,  6.123767403223466e1, 1.121012633939949e2,
    9.285288485892228e1, 3.198289277679660e1, 3.595376024148513, 6.160228690102976e-2};
const std::array<double, 8> k0qq = {1.0, 1.189963006673403e1, 5.027773590829784e1, 9.496513373427093e1,
    8.318077493230258e1, 3.181399777449301e1, 4.443672926432041, 1.408295601966600e-1};
const std::array<double, 5> k1pi = {0.5, 5.598072040178741e-2, 1.818666382168295e-3, 2.397509908859959e-5, 1.239567816344855e-7};
const std::array<double, 3> k1qi = {9.870202601341150e-1, 1.292092053534579e-2, 5.881933053917096e-5};
const std::array<double,5> k1p = {-3.079657578292062e-1, -8.109417631822442e-2, -3.477550948593604e-3, -5.385594871975406e-5, -3.110372465429008e-7};
const std::array<double, 3> k1q = {9.861813171751389e-1, 1.375094061153160e-2, 6.774221332947002e-5};
const std::array<double, 8> k1pp = {1.253314137315502, 1.457171340220454e1, 6.063161173098803e1, 1.147386690867892e2,
    1.040442011439181e2,  4.356596656837691e1, 7.265230396353690, 3.144418558991021e-1};
const std::array<double, 8> k1qq = {1.0, 1.125154514806458e1, 4.427488496597630e1, 7.616113213117645e1, 5.863377227890893e1,
    1.850303673841586e1, 1.857244676566022, 2.538540887654872e-2};

#endif
