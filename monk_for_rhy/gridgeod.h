//! \file gridgeod.h
#ifndef _GRIDGEOD_H
#define _GRIDGEOD_H
#include "diskobj.h"
#include "sim5lib.h"
#include "tridgeo.h"
#include "detector.h"

/*
void onaxisgrid(const long & ntheta, diskobj & disk, const double & a, const double & height, const double & vr,
    std::vector<double> & rearr, std::vector<double> & garr, std::vector<double> & domegaarr, bool todisk);
void lpoffaxis(const double & a, const double & height, const double & s, sim5metric & met,
	sim5tetrad & tet, std::vector<double> & four_v);
void lpoffaxis(const double & a, const double & height, const double & s, const double & omega, sim5metric & met,
	sim5tetrad & tet, std::vector<double> & four_v);

void offaxisgrid(const diskobj & disk, sim5tetrad & tet, long & ntheta, const long & nphi,
    std::vector<double> & rearr, std::vector<double> & garr, std::vector<double> & domegaarr);
void calmuarr_disktoslab(const double a, const double re, const double h, const double smax,
    const twodgeo & twodg, const long ntheta, const long nphi, std::vector<double> & muarr,
    std::vector<double> & garr, std::vector<double> & theta, std::vector<double> & domegaarr);
void calmuarr_inftoslab(const double a, const double h, const double smax, const double muobs,
    const twodgeo & twodg, const long nalpha, const long nbeta, const double alphamax,
    double & domega, std::vector<double> & muarr, std::vector<double> & garr);
void griddisk_nt(const double a, const double rout, const double m, const double mdot, 
    const long nr, const long ntheta, const long nphi, detector & det, const bool pol);
void grid3dcorona(const double a, const double rout, const double m, const double mdot, const double te, 
    const double mfp, const tridgeo & trid, const long nr, const long ntheta, 
    const long nphi, const long npack, detector & det);
void grid3dcorona_sp(const bool pol, const double a, const double rout, const double m, const double mdot, const double fcol, const double te, 
    const double mfp, const tridgeo & trid, const long nr, const long ntheta, const long nphi, 
    const long nphoton, std::ofstream & enfile, std::ofstream & weightfile, std::ofstream & muinffile, std::ofstream & nscafile,
    std::ofstream & qweightfile, std::ofstream & uweightfile);
void griddisk_nt_sp_pol(const double a, const double rout, const double m, const double mdot, const double fcol,
    const double thetamax, const long nr, const long ntheta, const long nphi, const long nphoton, std::ofstream & enfile,
    std::ofstream & weightfile, std::ofstream & muinffile, std::ofstream & qweightfile, std::ofstream & uweightfile);
void selfirr(const double a, const double m, const double mdot, const double fcol, 
    const double rtr, const long nr, const long ntheta, const long nphi, const long ne, 
    const double emin, const double emax);
void selfirr_pol(const double a, const double rout, const long nr, const long ntheta, const long nphi);
    */
void calmuobsarr(const double a, const double rin, const double rout, const double thetamax, const long nr, 
    const long ntheta, const long nphi, std::vector<double> & muobsarr, std::vector<double> & radius, std::vector<double> & theta,
    std::vector<double> & phi, std::vector<double> & rr1arr, std::vector<double> & muplusarr, std::vector<double> & larr,
    std::vector<double> & qarr, std::vector<double> & kmuobssignarr);
void griddisk_nt_sp_pol_writeparam(const bool pol, const double a, const double rout, const double thetamax,
    const long nr, const long ntheta, const long nphi, std::string paramfile);
void grid3dcorona_sp_writeparam(const double rin, const double rout, const double thetamax, const tridgeo & trid, const long nr, const long ntheta, const long nphi);
void register_sp_nt(const bool pol, const bool wsample, const std::string parafile, const double m, const double mdot, const double fcol,
    const long nphoton, std::mt19937_64 & gen, std::vector<std::ofstream> & ofiles);
void register_sp_tridsca(const bool pol, const bool wsample, const bool progress, const std::string parafile,
    const double m, const double mdot, const double fcol,
    const double te, const double mfp, const double dr, const long nphoton, std::vector<std::ofstream> & ofiles_inf,
    std::ofstream & nscafile_inf, std::vector<std::ofstream> & ofiles_disc, std::ofstream & nscafile_disc);
#endif
