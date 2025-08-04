//! \file calhotcross.h
#ifndef CALHOTCROSS_H
#define CALHOTCROSS_H
double readhotx(const double thetat, const double x);
double readhotx_local(const double thetat, const double x, const std::vector<double> & logthetatarr, 
    const std::vector<double> & logxarr, const std::vector<double> & hotxarr);
void init_hotx_local(std::vector<double> & logthetatarr, std::vector<double> & logxarr, std::vector<double> & hotxarr);
void init_hotx();

namespace HOTX {
    //! Directory containing data files
    const std::string HOTXDIR = "/home/hdw/data/monk/monk_for_rhy/data/";
    //! Lower boundary of \f$\theta_T\f$ in pre-calculation
    const double THETATMIN = 1e-4;
    //! Upper boundary of \f$\theta_T\f$ in pre-calculation
    const double THETATMAX = 10.;
    //! Lower boundary of \f$x\f$ in pre-calculation
    const double XMIN = 1e-12;
    //! Upper boundary of \f$x\f$ in pre-calculation
    const double XMAX = 1e6;
    //! Step size in \f${\rm log}(\theta_T)\f$
    const double DLOGTHETAT = 0.1;
    //! Step size in \f${\rm log}(x)\f$
    const double DLOGX = 0.1;
    const double LOGTHETATMIN = std::log10(THETATMIN), LOGTHETATMAX = std::log10(THETATMAX), LOGXMIN = std::log10(XMIN), LOGXMAX = std::log10(XMAX);
    //! Number of \f$\theta_T\f$
    const long NTHETAT = (LOGTHETATMAX - LOGTHETATMIN) / DLOGTHETAT + 1;
    //! Number of \f$x\f$
    const long NX = (LOGXMAX - LOGXMIN) / DLOGX + 1;
    
    //! Global vector used in reading data from disc; 
    extern std::vector<double> LOGXARR;
    //! Global vector used in reading data from disc; 
    extern std::vector<double> LOGTHETATARR;
    //! Global vector used in reading data from disc; 
    extern std::vector<double> HOTXARR;
};
#endif
