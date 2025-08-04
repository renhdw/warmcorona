/*! \file sphere.cpp

This is the source file for executable binary file `sphere`. With sphere one can calculate the spectrum of a spherical plasma cloud Comptonising
primary radiation. While the heavy work is done by \ref sphere_KN, this file offers the inferace for passing parameters.

#### Syntax:
    - `sphere`
    - `sphere parafile`
    
`sphere` will read a parameter file that contains all the necessary parameters. If being called without any arguments,
`sphere` assumes a default parameter file name `params.txt`.

An example of the parameter file:
\code
[physical]
# blackbody temperature 
tbb = 0.026
# if thermal=false, i.e., powerlaw primary:
# gamma = 2.
# emin = 0.1
# emax = 20.
# electron temperature 
te = 100.
# Thompson optical depth
tau = 0.2

[gridsize]
# number of photons
nphoton = 1000000

[option]
# whether the primary radiation is thermal or powerlaw
thermal = true
# scatter type
stype = 2
# whether to include recoil effect; note that this option will be overrided by stype
recoil = false
# unit of photon and electron temperature; if true, both temperatures are dimensionless; otherwise in [keV]
me = false
# spatial distribution of primary photons.
pattern = 0
\endcode
For detailed description of *stype* and *pattern*, see \ref sphere_KN and \ref scatter_bias_photon_KN_pol.

#### Output files
The results are written into several binary files. One can use *calspec* provided by \ref calspec.cpp to calculate the spectrum.
    - `en0.dat`: dimensionless energies of the superphotons; double precision
    - `weight.dat`: statistical weights of the superphotons; double precision
    - `muinf.dat`: \f$k_z\f$, the \f$z-\f$component of wave vector. This is only useful when `pattern=2`, where
        the primary radiation are injected into the corona along \f$+z\f$ direction, such that in this case
        \f$k_z\f$ is cosine of the angle between the primary radiation and the escaped superphoton; double precision
    - `nsca.dat`: number of scattering; `int` type (note the actual size of `int` type depends on compiler and/or architect)

Besides the binary files, `sphere` will also write an ASCII log file `sphere.log`.
*/
#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <experimental/filesystem>
#include "tridgeo.h"
#include "photon_dist.h"
#include "utils.h"
#include "scatter.h"
#include "const.h"
#include "./external/INIReader.h"
#include "calhotcross.h"

namespace fs = std::experimental::filesystem;
    
const double maxdl = 0.01;
int main(int argc, char * argv[]) {
	init_hotx();
    std::vector<std::string> args(argv, argv + argc);
    std::string parafile;
    long nphoton;
    double tau, te, tbb, gamma, emin, emax, dl, te_emass, tbb_emass, emin_emass, emax_emass;
    bool me, recoil, thermal;
    long pattern, stype;
    
    std::string en_unit;
    
    std::ofstream logfile("sphere.log");
    
    if (argc == 1) {
        parafile = "params.txt";
    } else {
        parafile = args[1];
    }
    
    if (!fs::exists(parafile)) {
        std::cerr << parafile + " does not exist!" << std::endl;
        return 1;
    }
    
    INIReader parfile(parafile);
    
    if (!parfile.GetInteger("option", "stype", stype)) {
        std::cerr << "missing parameter stype" << std::endl;
        return 1;
    }
    
    logfile << "# stype: 0 - Thompson; 1 - Compton, Thompson cross-section; 2 - Klein-Nishina cross section" << std::endl;
    logfile << "stype = " << stype << std::endl;
    
    if (stype == 0) {
        recoil = false;
    } else if (stype == 2) {
        recoil = true;
    } else {
        if (!parfile.GetBoolean("option", "recoil", recoil)) {
            std::cerr << "missing parameter recoil" << std::endl;
            return 1;
        }
    }
    logfile << "recoil = " << recoil << std::endl;
    
    if (!parfile.GetInteger("option", "pattern", pattern)) {
        std::cerr << "missing parameter pattern" << std::endl;
        return 1;
    }
    
    logfile << "# pattern: spatial distribution of primary photon" << std::endl;
    logfile << "# 0 - center; 1 - uniformly distributed; 2 - bottom" << std::endl;
    logfile << "pattern = " << pattern << std::endl;
    
    if (!parfile.GetReal("physical", "tau", tau)) {
        std::cerr << "missing parameter tau" << std::endl;
        return 1;
    }
    logfile << "tau = " << tau << std::endl;
    dl = 0.01 / tau;
    if (dl > maxdl) dl = maxdl;
    
    if (!parfile.GetBoolean("option", "me", me)) {
        std::cerr << "missing parameter me" << std::endl;
        return 1;
    }
    logfile << "# me: Energy unit; 0 - keV; 1: electron rest mass" << std::endl;
    logfile << "me = " << me << std::endl;
    en_unit = me ? std::string(" me c^2") : std::string(" keV");
    
    if (!parfile.GetReal("physical", "te", te)) {
        std::cerr << "missing parameter te" << std::endl;
        return 1;
    }
    logfile << "Te = " << te << en_unit << std::endl;
    
    if (!parfile.GetBoolean("option", "thermal", thermal)) {
        std::cerr << "missing parameter thermal" << std::endl;
        return 1;
    }
    logfile << "# thermal_primary: wheter the priamry radiation is thermal or non-thermal (powerlaw)" << std::endl;
    logfile << "thermal_primary = " << thermal << std::endl;
    
    if (thermal) {
        if (!parfile.GetReal("physical", "tbb", tbb)) {
            std::cerr << "missing parameter tbb" << std::endl;
            return 1;
        }
        logfile << "Tbb = " << tbb << en_unit << std::endl;
    } else {
        if (!parfile.GetReal("physical", "gamma", gamma)) {
            std::cerr << "missing parameter gamma" << std::endl;
            return 1;
        }
        logfile << "Gamma = " << gamma << std::endl;
        if (!parfile.GetReal("physical", "emin", emin)) {
            std::cerr << "missing parameter emin" << std::endl;
            return 1;
        }
        logfile << "Emin = " << emin << en_unit << std::endl;
        if (!parfile.GetReal("physical", "emax", emax)) {
            std::cerr << "missing parameter emax" << std::endl;
            return 1;
        }
        logfile << "Emax = " << emax << en_unit << std::endl;
    }
    
    te_emass = me ? te : te / ME_KEV;
    
    if (thermal) {
        tbb_emass = me ? tbb : tbb / ME_KEV;
    } else {
        emin_emass = me ? emin : emin / ME_KEV;
        emax_emass = me ? emax : emax / ME_KEV;
    }
    
    logfile << "dl = " << dl << std::endl;
    if (!parfile.GetInteger("gridsize", "nphoton", nphoton)) {
        std::cerr << "missing parameter nphoton" << std::endl;
        return 1;
    }
    logfile << "Nphoton = " << nphoton << std::endl;
    logfile.close();
    
    std::unique_ptr<photon_dist> distptr(nullptr);
    if (thermal) {
        distptr.reset(new bb(tbb_emass));
    } else {
        distptr.reset(new pl(gamma, emin_emass, emax_emass));
    }
    
    std::vector<std::string> filenames = {"en0.dat", "weight.dat", "muinf.dat", "nsca.dat"};
        
    std::vector<std::ofstream> ofiles;
    for (unsigned long i = 0; i < filenames.size(); ++i) {
        std::ofstream fs(filenames[i], std::ios::out | std::ofstream::binary);
        ofiles.push_back(std::move(fs));
    }
    
    sphere_KN(true, recoil, stype, pattern, nphoton, tau, te_emass, dl, 1. / (double)(nphoton), 
        *distptr, ofiles);
    
    for (unsigned long i = 0; i < filenames.size(); ++i) {
        ofiles[i].close();
    }
    std::cout << std::endl;
    
}
