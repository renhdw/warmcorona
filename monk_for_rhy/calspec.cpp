/*!
\file calspec.cpp
Calculate energy spectrum at a given inclination angle bin. This program first looks for `en0.dat' and `weight.dat' under a
a given directory. If these two files are not found, the program will send an error and terminate.
If `qweight.dat' and `uweight.dat' exist under the given directory, the polarisation spectrum will be computed as well.
Its function depends on 
number of parameters.

#### Input files
- `en0.dat`: double-precision binary file containing array of dimensionless photon energy
- `weight.dat`: double-precision binary file of photon statistical weight
- `muinf.dat`: double-precision binary file of \f$\mu_\infty\f$. Only required if inclination option is on.
- `qeight.dat`: double-precision binary file of \f$w \delta {\rm cos}2\psi\f$, where \f$w\f$ is the statistical weight,
   \f$\delta\f$ is the polarisation degree, and \f$\psi\f$ is the polarisation angle. Only required if polarisation option is on.
- `ueight.dat`: double-precision binary file of \f$w \delta {\rm sin}2\psi\f$. Only required if polarisation option is on.
- `nsca.dat`: binary file of photon scattering numbers with `int` type. Only required if scattering number option is on.

@param dir the directory containing all the input files
@param ne |ne| is size of energy bin. If ne < 0: logarithm energy bin; otherwise linear bin.
@param emin lower boundary of energy bin
@param emax upper boundary of energy bin
@param imin lower boundary of inclination bin
@param imax upper boundary of inclination bin
@param nsca if `nsca >= 0`: calculate spectrum for photons having experienced `nsca' numbers of scattering. If `nsca` is negative, 
    calculate spectrum for photons having experienced at least `-nsca` number of scattering. For instance, if `nsca==-4`, `calspec`
    select photons scattered at least 4 times.

#### Syntax
- calspec dir ne emin emax 
- calspec dir ne emin emax nsca
- calspec dir ne emin emax imin imax
- calspec dir ne emin emax imin imax nsca

#### Options
note that all options are set explicitely. Their on/off state depends on number of arguments and/or existence of some input files
- polarisation calculation: on if both `qweight.dat` and `uweight.dat` exist; otherwise off. If on, `calspec` will calculate two Stokes parameters
    `Q` and `U`, and stores them in `qflux.dat` and `uflux.dat` or `qflux_$nsca$.dat` and `uflux_$nsca$.dat`, respectively. 
- inclination selection: on if `nargs >= 7`. If on, `calspec` will first check if `muinf.dat` exists. `calspec` will select photons
    with inclination \f$i \in (imin, imax]\f$, and normalise flux (and qflux, uflux if polarisation is on) 
    with a factor \f$4\pi/d\Omega = 2/[{\rm cos}(imin) - {\rm cos}(imax)]\f$.
- scattering number selection: on if `nargs == 6 || nargs == 8`. If on, `calspec` will first check if `nsca.dat` exists. 
    then select photons with: number of scattering equals to `nsca` (if `nsca >= 0`); or number of scattering no less then `nsca` (if `nsca < 0`). 

#### Output files: all output files will be put in the **current(!!)** directory.
- `en.dat`: double-precision binary file of energy bin centroid energy in [keV]
- `de.dat`: double-precision binary file of energy bin width [keV]
- `flux.dat`: double-precision binary file of flux of all scattering orders, 
       in \f$[\rm photons\ s^{-1}\ keV^{-1}]\f$; output if scattering number selection is off.
- `flux_$nsca$.dat`: double-precision binary file of flux of photons that scattered `nsca` times, in \f$[\rm photons\ s^{-1}\ keV^{-1}]\f$;
        output if scattering number selection is on. If `nsca=1`, then the file name would be `flux_1sca.dat`.
- `qflux.dat`: double-precision binary file of Stokes parameter \f$Q\f$. Output if polarisation option is on and scattering number option is off.
- `uflux.dat`: double-precision binary file of Stokes parameter \f$U\f$. Output if polarisation option is on and scattering number option is off.
- `qflux_$nsca$.dat`: double-precision binary file of Stokes parameter \f$Q\f$.
        Output if polarisation option is on and scattering number option is on.
- `uflux_$nsca$.dat`: double-precision binary file of Stokes parameter \f$U\f$
        Output if polarisation option is on and scattering number option is on.
*/

#include <algorithm>
#include <vector>
#include <string>
#include <cmath>
#include <experimental/filesystem>
#include "utils.h"
#include "const.h"
#include "detector.h"

namespace fs = std::experimental::filesystem;
int main(int argc, char * argv[]) {
    std::vector<std::string> arglist(argv, argv + argc);
    std::string dirname, enfile, wfile, qfile, ufile, nscafile, mufile;
    double emin, emax, mumin, mumax, imin, imax, de;
    long ne;
    std::vector<double> en0, weight, muinf, qweight, uweight, dearr;
    std::vector<int> nscaarr;
    long size_1part = 1e8;
    bool pol;
    std::vector<double> poldeg, polang;
    double polflux;
    int nsca;
    
    std::ofstream ofile("calspec.log");
    for (int i = 0; i < argc; ++i)
        ofile << argv[i] << " ";
    ofile << std::endl;
    ofile.close();
    
    dirname = arglist[1];
    enfile = dirname + "/en0.dat";
    wfile = dirname + "/weight.dat";
    qfile = dirname + "/qweight.dat";
    ufile = dirname + "/uweight.dat";
    nscafile = dirname + "/nsca.dat";
    mufile = dirname + "/muinf.dat";
    
    ne = std::stol(arglist[2]);
    emin = std::stod(arglist[3]);
    emax = std::stod(arglist[4]);
    
    if (!fs::exists(enfile) || !fs::exists(wfile)) {
        std::cerr << "en0.dat or weight.dat does not exist!" << std::endl;
        return 1;
    }
    pol = (fs::exists(qfile) && fs::exists(ufile));
    
    if (argc == 5 || argc == 6) {
        mumin = 0;
        mumax = 0.;
        if (argc == 6)
            nsca = std::stoi(arglist[5]);
    } else if (argc == 7) {
        imin = std::stod(arglist[5]);
        imax = std::stod(arglist[6]);
        mumin = std::cos(imax / 180. * M_PI);
        mumax = std::cos(imin / 180. * M_PI);
    } else {
        imin = std::stod(arglist[5]);
        imax = std::stod(arglist[6]);
        mumin = std::cos(imax / 180. * M_PI);
        mumax = std::cos(imin / 180. * M_PI);
        nsca = std::stoi(arglist[7]);
    }
    
    singledet_sp det(mumin, mumax, ne, emin, emax);
    singledet_sp qdet(mumin, mumax, ne, emin, emax);
    singledet_sp udet(mumin, mumax, ne, emin, emax);
    
    // get file size; if > 100M, we will split the files into segments of 100M
    auto filesize = get_filesize(enfile);
    auto wfilesize = get_filesize(wfile);
    if (filesize != wfilesize) {
        std::cerr << "Energy and weight file size not compatible!" << std::endl;
        return 1;
    }
    
    if (pol) {
        auto qfilesize = get_filesize(qfile);
        auto ufilesize = get_filesize(ufile);
        if (filesize != qfilesize || filesize != ufilesize) {
            std::cerr << "Energy and q/u weight file size not compatible!" << std::endl;
            return 1;
        }
    }
    
    if (filesize < size_1part) {
        rdoublevec(en0, enfile);
        rdoublevec(weight, wfile);
        if (pol) {
            rdoublevec(qweight, qfile);
            rdoublevec(uweight, ufile);
        }
        
        std::transform(en0.cbegin(), en0.cend(), en0.begin(), std::bind2nd(std::multiplies<double>(), ME_KEV));
        if (argc == 5) {
            det.calflux(en0, weight);
            if (pol) {
                qdet.calflux(en0, qweight);
                udet.calflux(en0, uweight);
            }
        } else if (argc == 6) {
            rdoublevec(nscaarr, nscafile);
            det.calflux(en0, weight, nscaarr, nsca);
            std::cout << "we are here" << std::endl;
        } else if (argc == 7) {
            rdoublevec(muinf, mufile);
            det.calflux(en0, weight, muinf);
            if (pol) {
                qdet.calflux(en0, qweight, muinf);
                udet.calflux(en0, uweight, muinf);
            }
        } else {
            std::cout << "here" << std::endl;
            rdoublevec(muinf, mufile);
            rdoublevec(nscaarr, nscafile);
            std::cout << "nscaarr.size() = " << nscaarr.size() << std::endl;
            std::cout << "muinf.size() = " << muinf.size() << std::endl;
            det.calflux(en0, weight, muinf, nscaarr, nsca);
            if (pol) {
                qdet.calflux(en0, qweight, muinf, nscaarr, nsca);
                udet.calflux(en0, uweight, muinf, nscaarr, nsca);
            }
        }
            
    } else {
        long npart = std::ceil((double)(filesize) / (double)(size_1part));
        std::cout << "Reading files part by part; size: " << size_1part / 1e6 << " MB" << std::endl;
        for (long i = 0; i < npart; ++i) {
            rdoublevec(en0, enfile, i * size_1part / 8, size_1part / 8);
            //std::cout << "en0 = " << en0 << std::endl;
            rdoublevec(weight, wfile, i * size_1part / 8, size_1part / 8);
            std::transform(en0.cbegin(), en0.cend(), en0.begin(), std::bind2nd(std::multiplies<double>(), ME_KEV));
            
            if (pol) {
                rdoublevec(qweight, qfile, i * size_1part / 8, size_1part / 8);
                rdoublevec(uweight, ufile, i * size_1part / 8, size_1part / 8);
            }
            
            if (argc == 5) {
                det.calflux(en0, weight);
                if (pol) {
                    qdet.calflux(en0, qweight);
                    udet.calflux(en0, uweight);
                }
            } else if (argc == 7) {
                rdoublevec(muinf, mufile, i * size_1part / 8, size_1part / 8);
                det.calflux(en0, weight, muinf);
                if (pol) {
                    qdet.calflux(en0, qweight, muinf);
                    udet.calflux(en0, uweight, muinf);
                }
                en0.clear();
                muinf.clear();
                weight.clear();
                if (pol) {
                    qweight.clear();
                    uweight.clear();
                }
            } else {
                rdoublevec(muinf, mufile, i * size_1part / 8, size_1part / 8);
                rdoublevec(nscaarr, nscafile, i * size_1part / 8, size_1part / 8);
                det.calflux(en0, weight, muinf, nscaarr, nsca);
                if (pol) {
                    qdet.calflux(en0, qweight, muinf, nscaarr, nsca);
                    udet.calflux(en0, uweight, muinf, nscaarr, nsca);
                }
                en0.clear();
                muinf.clear();
                weight.clear();
                if (pol) {
                    qweight.clear();
                    uweight.clear();
                }
            }
        }
    }
    if (argc != 5 && argc != 6) {
        det.norm();
        if (pol) {
            qdet.norm();
            udet.norm();
        }
    }
    
    det.norm_flux();
    if (argc == 6 || argc == 8) {
        std::cout << "nsca = " << nsca << std::endl;
        det.writename("flux_" + std::to_string(nsca) + "sca.dat");
    } else {
        det.write();
    }
        
    if (pol) {
        qdet.norm_flux();
        udet.norm_flux();
        if (argc != 8) {
            qdet.writename("qflux.dat");
            udet.writename("uflux.dat");
        } else {
            qdet.writename("qflux_" + std::to_string(nsca) + "sca.dat");
            udet.writename("uflux_" + std::to_string(nsca) + "sca.dat");
        }
        
        // we further calculate polarization degree and angle
        poldeg.resize(det.en.size());
        polang.resize(det.en.size());
        for (size_t ie = 0; ie < det.en.size(); ++ie) {
            polflux = std::sqrt(qdet.flux[ie] * qdet.flux[ie] + udet.flux[ie] * udet.flux[ie]);
            poldeg[ie] = polflux / det.flux[ie];
            polang[ie] = std::acos(qdet.flux[ie] / polflux) * 0.5;
            if (udet.flux[ie] < 0.) {
                polang[ie] = M_PI - polang[ie];
            }
        }
        if (argc != 8) {
            wdoublevec(poldeg, "poldeg.dat");
            wdoublevec(polang, "polang.dat");
        } else {
            wdoublevec(poldeg, "poldeg_" + std::to_string(nsca) + "sca.dat");
            wdoublevec(polang, "polang_" + std::to_string(nsca) + "sca.dat");
        }
    }
    
    if (ne > 0) {
        de = (emax - emin) / (double)(ne);
        dearr = std::vector<double>(ne, de);
    } else {
        double dloge = std::log10(det.en[1] / det.en[0]);
        double sqrt_Eratio = std::pow(10., 0.5 * dloge);
        sqrt_Eratio = sqrt_Eratio - 1. / sqrt_Eratio;
        
        dearr.resize(det.en.size());
        std::transform(det.en.cbegin(), det.en.cend(), dearr.begin(), std::bind2nd(std::multiplies<double>(), sqrt_Eratio));
    }
    
    wdoublevec(dearr, "de.dat");
}
