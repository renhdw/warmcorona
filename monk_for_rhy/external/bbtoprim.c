#define NE     500
#define NTHCOMP_DE 1e-2
#define NTHCOMP_EMIN0 1e-3
#define NTHCOMP_EMAX0 1e+3
#define E_MIN  1e-3
#define E_MAX  100.
//constants
#define H_KEVS 4.13566743e-18
#define K_KEVK 8.6174e-8
#define C_MS   2.99792458e8
#define MSOLAR 1.989e+30
#define G      6.6743e-11
#define SIGMA  5.6704e-8
#define SIGMAP 1.5204606e+15
#define YEAR   31557600.0
#define MPC_2  1.05e-49
//kpc in cm
#define KPC    3.0856776e+21
#define ERG 6.241509e8      //erg in keV
// Ledd is in keV/s (not W)
#define LEDD 7.8643017e+46
//kev in SI units (not in ergs!)
#define KEV 1.6022e-16
#define SIGMAT 0.66524574e-24//cm^(2)
#define PI 3.14159265358979
#define PI2 6.2831853071795865

// table file
#define LAMP "/home/wdzhang/data/bbtoprim/test/KBHlamp_q.fits"

static float  *radius;
static long   nrad;
static double *gfac, *cosin, *phiph, *transf_d;

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"
#include "../donthcomp.h"

void BBtoPrim(int locallx, int nthopt, int bolometric, const double gamma, const double a, const double h, double rout, const double m, 
    const double arate, const double theta_o, const double fcol, const double ecut, 
    const double lxemin, const double lxemax, const double lx, const long nr);

int main(int argc, char * argv[]) {
    int locallx = atoi(argv[1]);
    int nthopt = atoi(argv[2]);
    int bolometric = atoi(argv[3]);
    double h = atof(argv[4]);
    double gamma = atof(argv[5]);
    double lx = atof(argv[6]);
    double arate = 0.00693;
    //double arate = 6.83e-9;
    double a = 0.998, m = 1e7, theta_o = 30.;
    double fcol = 2.4, ecut = 100., lxemin = 2., lxemax = 10., rout = 1000.;
    long nr = 150;
    BBtoPrim(locallx, nthopt, bolometric, gamma, a, h, rout, m, arate, theta_o, fcol, ecut, lxemin, lxemax, lx, nr);
}

void BBtoPrim(int locallx, int nthopt, int bolometric, const double gamma, const double a, const double h, double rout, const double m, 
    const double arate, const double theta_o, const double fcol, const double ecut, 
    const double lxemin, const double lxemax, double lx, const long nr) {

/* Let's declare static variables whose values should be remembered for next
   run in XSPEC */
    //static char   pname[128]="KYDIR";
    static long   nrh, nh, nincl;
    static float  *r_horizon, *height, *incl, *dWadWo, *dWadSd, *q, *pr;
    double transf_o;

    FILE *fw;
    double rin, Fn_in, Fn_out, Fe_in, Fe_out;
    double Tbb = 0., Fmax;
    double Rg, D;
    double am2, pom, pom1, pom2, pom3, rms, r_plus, x, x0, x1, x2, x3, Ccal, Lcal, arcosa3, Tnorm;
    double h_rh;
    double ttmp, ttmp1, utmp, utmp1, vtmp, vtmp1, y1, y2, y3, y4, y5, y6, y7, y8;
    double pr_final, q_final;
    double r, r2, delta, ULt, tmp1, Ut, U_phi, U_r, Ur, UrU_r, Lms, Ems;
    double temp, tau, Ut_rms;
    double g_L, drad, radrat, rnow, radratio, factor;
    double lensing, gfactor;
    double lx0;
    int    ie, i, ig;
    int    imin, imax, irh0, ih0, ith0, ir, ir0;//, it0;
    //int    irh, ih;
    //nthcomp
    long nthcomp_ne;
    int option;
    float *nthcomp_ear, *nthcomp_far, *nthcomp_fer, nthcomp_param[6];
    double nthcomp_lum, tmp;
    const double gamt[10]={1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00, 3.25, 
                            3.50};
    const double taut[10]={2.500, 1.800, 1.200, 0.850, 0.600, 0.400, 0.260, 0.180, 
                        0.120, 0.075};

    // these are needed to work with a fits file...
    fitsfile *fptr;
    char     tables_file[255];
    int      hdutype = 2;
    int      colnum = 1;
    long     frow = 1, felem = 1, nelems, nrow;
    float    float_nulval = 0.;
    int      nelements1, nelements2;
    int      ihorizon, irow, anynul, status = 0;//, maxdim=1000, naxis;
    
    double ear[NE + 1], photar[NE], de[NE];//, photer[NE];
    for (ie = 0; ie <= NE; ie++) {
    //  ear[ie] = E_MIN + ie * (E_MAX-E_MIN) / NE;
        ear[ie] = E_MIN * pow(E_MAX/E_MIN, ((double) ie) / NE);
    }

    // Let's initialize parameters for subroutine ide()
    // a/M - black hole angular momentum
    am2 = a * a;
    r_plus = 1. + sqrt(1. - am2);
    //r_minus =1.-sqrt(1.-am2);
    pom1 = pow(1. + a, 1. / 3.);
    pom2 = pow(1. - a, 1. / 3.);
    pom3 = pow(1. - am2, 1. / 3.);
    pom = 1. + pom3 * (pom1 + pom2);
    pom1 = sqrt(3. * am2 + pom * pom);
    if (a >= 0.)
        rms= 3. + pom1 - sqrt((3. - pom) * (3. + pom + 2. * pom1));
    else 
        rms = 3. + pom1 + sqrt((3. - pom) * (3. + pom + 2. * pom1));
    //Ut_rms=(rms*rms-2.*rms+am*sqrt(rms))/rms/sqrt(rms*rms-3.*rms+2.*am*sqrt(rms));
    Ut_rms = (4 * (sqrt(rms) - a) + a) / sqrt(3.) / rms;
    x0 = sqrt(rms);
    arcosa3 = acos(a) / 3.;
    x1 = 2 * cos(arcosa3 - PI / 3.);
    x2 = 2 * cos(arcosa3 + PI / 3.);
    x3 = -2 * cos(arcosa3);
    // theta_o - observer inclination
    
    rin = rms - r_plus;
    rout = rout - r_plus;
    
    Rg = G * MSOLAR * m / (C_MS * C_MS);
    h_rh = h - r_plus;
    // PhoIndex - power-law energy index of the lamp emission
    // find out corresponding tau with interpolation
    ig=0;
    while (gamma > gamt[ig+1] && ig < 9) 
        ig++;
    tau = taut[ig]+(taut[ig+1]-taut[ig])/(gamt[ig+1]-gamt[ig])*(gamma-gamt[ig]);
    //fprintf(stdout,"tau = %4.3f\n",tau);
    
/******************************************************************************/
// Let's read the lamp post tables
// The status parameter must always be initialized.
    status = 0;
// Open the FITS file for readonly access
// - if set try KYDIR directory, otherwise look in the working directory
//   or in the xspec directory where tables are usually stored...
    //sprintf(kydir, "%s", FGMSTR(pname));
    //if (strlen(kydir) == 0)
    sprintf(tables_file, "%s", LAMP);
    
// Let's read the 'KBHlamp_q' fits file
// The status parameter must always be initialized.
    status = 0;
    ffopen(&fptr, tables_file, READONLY, &status);
    if (status) {
        //sprintf(tables_file, "%s%s", FGMODF(), LAMP);
        status = 0;
        ffopen(&fptr, tables_file, READONLY, &status);
    }
    if (status) {
        if (status) 
            ffrprt(stderr, status);
        ffclos(fptr, &status);
        //xs_write("\nkynBBtoPrim: set the KYDIR to the directory with the KY tables",5);
        return;
    }
// Let's read tables (binary tables => hdutype=2)
// Move to the extension 'r_horizon' and read its values
    ffmrhd(fptr, 1, &hdutype, &status);
    ffgnrw(fptr, &nrh, &status);
//******************************************************************************
//  fprintf(stdout,"nrh = %ld\n",nrh);
//******************************************************************************   
// Allocate memory for r_horizon...
    if ((r_horizon = (float *) malloc(nrh * sizeof(float))) == NULL) {
        //xs_write("kynBBtoPrim: Failed to allocate memory for tmp arrays.", 5);
        return;
    }
// Read the data in the 'r_horizon' table
    nelems = nrh;
// FTGCVE reads the VALUES from the first column.
    ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, r_horizon,
        &anynul, &status);
//******************************************************************************
//  for ( i=0; i<nrh; i++)fprintf(stdout,"%f\n",r_horizon[i]);
//******************************************************************************   
// Move to the extension 'height' and read its values
    ffmrhd(fptr, 1, &hdutype, &status);
    ffgnrw(fptr, &nh, &status);
//******************************************************************************
//  fprintf(stdout,"nh = %ld\n",nh);
//******************************************************************************   
// Allocate memory for height...
    if ((height = (float *) malloc(nh * sizeof(float))) == NULL) {
        //xs_write("kynBBtoPrim: Failed to allocate memory for tmp arrays.", 5);
        return;
    }
// Read the data in the 'height' table
    nelems = nh;
// FTGCVE reads the VALUES from the first column.
    ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, height, &anynul, &status);
//******************************************************************************
//  for ( i=0; i<nh; i++)fprintf(stdout,"%f\n",height[i]);
//******************************************************************************   
// Move to the extension 'inclination' and read its values
    ffmrhd(fptr, 1, &hdutype, &status);
    ffgnrw(fptr, &nincl, &status);
//******************************************************************************
   //fprintf(stdout,"nincl = %ld\n",nincl);
//******************************************************************************   
// Allocate memory for inclination...
    if ((incl = (float *) malloc(nincl * sizeof(float))) == NULL) {
        //xs_write("kynBBtoPrim: Failed to allocate memory for tmp arrays.", 5);
    return;
    }
// Read the data in the 'inclination' table
    nelems = nincl;
// FTGCVE reads the VALUES from the first column.
    ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, incl,
        &anynul, &status);
//******************************************************************************
//  for ( i=0; i<nincl; i++)fprintf(stdout,"%f\n",incl[i]);
//******************************************************************************   
// Move to the extension 'r_rh' and read its values
    ffmrhd(fptr, 1, &hdutype, &status);
    ffgnrw(fptr, &nrad, &status);
//******************************************************************************
    //fprintf(stdout,"nrad = %ld\n",nrad);
//******************************************************************************   
// Allocate memory for radius...
    if ((radius = (float *) malloc(nrad * sizeof(float))) == NULL) {
        //xs_write("kynBBtoPrim: Failed to allocate memory for tmp arrays.", 5);
        return;
    }
// Read the data in the 'r_rh' table
    nelems = nrad;
// FTGCVE reads the VALUES from the first column.
    ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, radius, &anynul, &status);
//******************************************************************************
//  for ( i=0; i<nrad; i++)fprintf(stdout,"%f\n",radius[i]);
//******************************************************************************   
    if ((gfac = (double *) malloc(nrad * sizeof(double))) == NULL) {
        //xs_write("kynBBtoPrim: Failed to allocate memory for tmp arrays.", 5);
        return;
    }
    if ((cosin = (double *) malloc(nrad * sizeof(double))) == NULL) {
        //xs_write("kynBBtoPrim: Failed to allocate memory for tmp arrays.", 5);
        return;    
    }
    if ((phiph = (double *) malloc(nrad * sizeof(double))) == NULL) {
        //xs_write("kynBBtoPrim: Failed to allocate memory for tmp arrays.", 5);
        return;
    }
    if ((transf_d = (double *) malloc(nrad * sizeof(double))) == NULL) {
        //xs_write("kynBBtoPrim: Failed to allocate memory for tmp arrays.", 5);
        return;    
    }
// Let's read the tables for dWadWo, q, p^r and dWadSd
// allocate memory for the arrays
    if ((dWadWo = (float *) malloc(nincl * nh * nrh * sizeof(float))) == NULL) {
        //xs_write("kynBBtoPrim: Failed to allocate memory for tmp arrays.", 5);
        return;
    }
    if ((q = (float *) malloc(nrad * nh * nrh * sizeof(float))) == NULL) {
        //xs_write("kynBBtoPrim: Failed to allocate memory for tmp arrays.", 5);
        return;
    }
    if ((pr = (float *) malloc(nrad * nh * nrh * sizeof(float))) == NULL) {
        //xs_write("kynBBtoPrim: Failed to allocate memory for tmp arrays.", 5);
        return;
    }
    if ((dWadSd = (float *) malloc(nrad * nh * nrh * sizeof(float))) == NULL) {
        //xs_write("kynBBtoPrim: Failed to allocate memory for tmp arrays.", 5);
        return;
    }
// read the tables
    for (ihorizon = 0; ihorizon < nrh; ihorizon++) {
        ffmrhd(fptr, 1, &hdutype, &status);
/*  to read the file only once we have to read in blocks (all columns
    from the extension are put to buffer together)
    let's find out how many rows are going to be read into the buffer */
        ffgrsz(fptr, &nrow, &status);
        nelements1 = nrow * nincl;
        nelements2 = nrow * nrad;
        for (irow = 0; irow < nh; irow += nrow) {
//    the last block to read may be smaller:
            if ((nh - irow) < nrow) {
                nelements1 = (nh - irow) * nincl;
                nelements2 = (nh - irow) * nrad;
            }
            ffgcv(fptr, TFLOAT, 1, irow+1, 1, nelements1, &float_nulval, 
                &dWadWo[irow * nincl + nh * nincl * ihorizon],&anynul, &status);
            ffgcv(fptr, TFLOAT, 2, irow+1, 1, nelements2, &float_nulval, 
                &q[irow * nrad + nh * nrad * ihorizon],&anynul, &status);
            ffgcv(fptr, TFLOAT, 3, irow+1, 1, nelements2, &float_nulval, 
                &pr[irow * nrad + nh * nrad * ihorizon],&anynul, &status);
            ffgcv(fptr, TFLOAT, 4, irow+1, 1, nelements2, &float_nulval, 
                &dWadSd[irow * nrad + nh * nrad * ihorizon],&anynul, &status);
        }
    }
// The FITS file must always be closed before exiting the program.
    ffclos(fptr, &status);
// Let's interpolate the tables to desired spin and height
    // given am->r_plus, find the corresponding index in r_horizon[]:
    imin = 0;
    imax = nrh;
    irh0 = nrh / 2;
    
    while ((imax - imin) > 1) {
        if (r_plus >= r_horizon[irh0-1])
            imin = irh0;
        else 
            imax = irh0;
        irh0 = (imin + imax) / 2;
    }
    if (irh0 == 0) 
        irh0 = 1;
//if ((imax == nrh) && (r_plus > r_horizon[nrh-1])) irh0 = nrh;
    ttmp = (r_plus - r_horizon[irh0 - 1])
        / (r_horizon[irh0] - r_horizon[irh0 - 1]);
    ttmp1 = 1. - ttmp;
// given h, find the corresponding index in height[]:
    imin = 0;
    imax = nh;
    ih0 = nh / 2;
    while ((imax - imin) > 1) {
        if (h_rh >= height[ih0 - 1])
            imin = ih0;
        else 
            imax = ih0;
        ih0 = (imin + imax) / 2;
    }
    if (ih0 == 0)
        ih0 = 1;
//if ((imax == nh) && (h_rh > height[nh-1])) ih0 = nh;
    utmp = (h_rh - height[ih0 - 1]) / (height[ih0] - height[ih0 - 1]);
    utmp1 = 1. - utmp;
    // given thetaO, find the corresponding index in incl[]:
    imin = 0;
    imax = nincl;
    ith0 = nincl / 2;
    while ((imax - imin) > 1) {
        if (theta_o >= incl[ith0 - 1])
            imin = ith0;
        else 
            imax = ith0;
        ith0 = (imin + imax) / 2;
    }
    
    if (ith0 == 0)
        ith0 = 1;
//if ((imax == nincl) && (thetaO > incl[nincl - 1])) ith0 = nincl;
    vtmp = (theta_o - incl[ith0 - 1]) / (incl[ith0] - incl[ith0 - 1]);
    vtmp1 = 1. - vtmp;
// transfer function from the axis to the observer
    y1 = dWadWo[ith0 - 1 + nincl * (ih0 - 1) + nincl * nh * (irh0 - 1)];
    y2 = dWadWo[ith0 - 1 + nincl * (ih0 - 1) + nincl * nh * irh0];
    y3 = dWadWo[ith0 - 1 + nincl * ih0 + nincl * nh * irh0];
    y4 = dWadWo[ith0 - 1 + nincl * ih0 + nincl * nh * (irh0 - 1)];
    y5 = dWadWo[ith0 + nincl * (ih0 - 1) + nincl * nh * (irh0 - 1)];
    y6 = dWadWo[ith0 + nincl * (ih0 - 1) + nincl * nh * irh0];
    y7 = dWadWo[ith0 + nincl * ih0 + nincl * nh * irh0];
    y8 = dWadWo[ith0 + nincl * ih0 + nincl * nh * (irh0 - 1)];
    transf_o = (vtmp1 * (utmp1 * (ttmp1 * y1 + ttmp * y2) + utmp *
        (ttmp * y3 + ttmp1 * y4)) + vtmp * (utmp1 *
        (ttmp1 * y5 + ttmp * y6) + utmp * (ttmp * y7 + ttmp1 * y8)));
    //fprintf(stdout, "transf_o = %f\n", transf_o);
    //fprintf(stderr, "%i %i\n", ih0, ith0);
        
    for (i = 0; i < nrad; i++) {
// q from the axis to the disc
        y1 = q[i + nrad * (ih0 - 1) + nrad * nh * (irh0 - 1)];
        y2 = q[i + nrad * (ih0 - 1) + nrad * nh * irh0];
        y3 = q[i + nrad * ih0 + nrad * nh * irh0];
        y4 = q[i + nrad * ih0 + nrad * nh * (irh0 - 1)];
        q_final = utmp1 * (ttmp1 * y1 + ttmp *y2) + 
            utmp * (ttmp * y3 + ttmp1 * y4);
// pr at the disc
        y1 = pr[i + nrad * (ih0 - 1) + nrad * nh * (irh0 - 1)];
        y2 = pr[i + nrad * (ih0 - 1) + nrad * nh * irh0];
        y3 = pr[i + nrad * ih0 + nrad * nh * irh0];
        y4 = pr[i + nrad * ih0 + nrad * nh * (irh0 - 1)];
        pr_final = utmp1 * (ttmp1 * y1 + ttmp * y2) + 
            utmp * (ttmp * y3 + ttmp1 * y4);
// Due to changing the direction of geodesic we have to change the sign in pr
        pr_final = -pr_final;
// temporary variables
        r = r_plus + radius[i];
        r2 = r * r;
        delta = r2 - 2. * r + am2;
        ULt = sqrt((h * h + am2) / (h * h - 2. * h + am2));
        if (r >= rms) {
            tmp1 = sqrt(r2 - 3. * r + 2. * a * sqrt(r));
            Ut = (r2 + a * sqrt(r)) / r / tmp1;
            U_phi = (r2 + am2 - 2. * a * sqrt(r)) / sqrt(r) / tmp1;
            U_r = 0.;
            Ur = 0.;
            UrU_r = 0.;
        }
        else {
            tmp1 = sqrt(rms * (rms - 3.) + 2. * a * sqrt(rms));
            Lms = (rms * rms + am2 - 2. * a * sqrt(rms)) / sqrt(rms) / tmp1;
            Ems = (rms * (rms - 2.) + a * sqrt(rms)) / rms / tmp1;
            Ut = (Ems * (r * (r2 + am2) + 2. * am2) - 2. * a * Lms) / r / delta;
            U_phi = Lms;
            UrU_r = -1. + ((r2 + am2 + 2. * am2 / r) * Ems * Ems - 4. * a /
                    r * Ems * Lms - (1. - 2. / r) * Lms * Lms) / delta;
            if (UrU_r < 0.)
                UrU_r = 0.;
            U_r = -sqrt(UrU_r / delta) * r;
            Ur = -sqrt(delta * UrU_r) / r;
        }
        tmp1 = Ut - pr_final * U_r;
// gfactor  from the disc to the axis
        gfac[i] = ULt / tmp1;
// cosin at the disc
        cosin[i] = q_final / r / tmp1;
// phip_i at the disc
        phiph[i] = atan2(-U_phi, r * (pr_final - Ur * tmp1));
// dWadSd from the axis to the disc
        y1 = dWadSd[i + nrad * (ih0 - 1) + nrad * nh * (irh0 - 1)];
        y2 = dWadSd[i + nrad * (ih0 - 1) + nrad * nh * irh0];
        y3 = dWadSd[i + nrad * ih0 + nrad * nh * irh0];
        y4 = dWadSd[i + nrad * ih0 + nrad * nh * (irh0 - 1)];
        transf_d[i] = utmp1 * (ttmp1 * y1 + ttmp * y2) + utmp
            * (ttmp * y3 + ttmp1 * y4);
    }
/******************************************************************************/
    //#ifdef OUTSIDE_XSPEC
    // Let's write input parameters to a text file
    fw = fopen("bbtoprim.log", "w");
    fprintf(fw, "a/M          %12.6f\n", a);
    fprintf(fw, "theta_o      %12.6f\n", theta_o);
    fprintf(fw, "rin          %12.6f\n", rin+r_plus);
    fprintf(fw, "rout         %12.6f\n", rout+r_plus);
    fprintf(fw, "BBmass       %12.6g\n", m);
    fprintf(fw, "arate        %12.6g\n", arate);
    fprintf(fw, "fcol        %12.6f\n", fcol);
    fprintf(fw, "height       %12.6f\n", h);
    fprintf(fw, "PhoIndex     %12.6f\n", gamma);
    fprintf(fw, "tau          %12.6f\n", tau);
    fprintf(fw, "sigma_e      %12e\n", tau/SIGMAT);
    fprintf(fw, "Ledd [erg/s] %12.6g\n", m*LEDD*KEV*1e+7);
    fprintf(fw, "eta*Mdot*c^2/Ledd %12.6g\n", arate*MSOLAR/YEAR*C_MS*C_MS*(1-Ut_rms)/m/LEDD/KEV);
    fprintf(fw, "Lx_obs/Ledd  %12.6g\n", lx / (LEDD*KEV*1e+7*m));
    fprintf(fw, "Ecut         %12.6f\n", ecut);
    fprintf(fw, "nrad         %12ld\n",  nrad);
    fprintf(fw, "r_horizon    %12.6f\n",r_plus);
    fprintf(fw, "r_ms         %12.6f\n",rms);
    //fprintf(fw, "emin         %12.6f\n",E_MIN);
    //fprintf(fw, "emax         %12.6f\n",E_MAX);
    fclose(fw);
    //#endif
/******************************************************************************/

// for (ie = 0; ie < ne; ie++) {
//   flx[ie] = 2 * pow((ear[ie] + ear[ie+1])/2. / H_KEVS / C_MS, 2.) / H_KEVS /
//             pow(f_col, 3.);
// }
    Tnorm = fcol * K_KEVK * pow(C_MS, 1.5) * pow(3. / (8. * PI * G * G * MSOLAR * YEAR * SIGMA), 0.25);

    g_L = sqrt(1. - 2. * h / (am2 + h * h));
//We work with local Lx in keV/s thus convert the units from Ledd or erg/s
    if (lx > 0.) {
        lx *= ERG;
    } else {
        lx *= (-1. * LEDD * m);
    }
//convert to intrinsic luminosity in case it is given at infinity
    //lx0 = 
    //fprintf(stdout, "Lamp-post frame: lx0 = %E keV/s\n", lx0);
    lx0 = (locallx != 0) ? lx : lx / g_L / g_L / transf_o;
    
    radratio=pow((rout+r_plus)/(rin+r_plus),1./nr);
//fw = fopen("kynbbtoprim_radius.dat", "w");
    Fn_in=0.;
    Fe_in=0.;
    for (ie = 0; ie < NE; ie++)
        photar[ie]=0.;
    for (ir = 0; ir < nr; ir++){
        radrat=pow(radratio,ir);
        rnow=0.5*(rin+r_plus)*(radratio+1.)*radrat;
        drad=(rin+r_plus)*(radratio-1.)*radrat;
        
        //factor0 = h / (rnow*rnow + h*h) / rnow;
        factor= 2. * PI * drad * h / (rnow*rnow + h*h);  
//given rnow, find corresponding indices in radius:
        imin = 0;
        imax = nrad;
        ir0 = nrad / 2;
        while ((imax - imin) > 1) {
            if (rnow >= (radius[ir0 - 1] + r_plus))
                imin = ir0;
            else 
                imax = ir0;
            ir0 = (imin + imax) / 2;
        }
        if ((ir0 == nrad) && (rnow == radius[nrad-1] + r_plus))
            ir0--;
        ttmp = (rnow - radius[ir0 - 1] - r_plus) / (radius[ir0] - radius[ir0 - 1]);
        ttmp1 = 1. - ttmp;
//Let's interpolate gfactor between two radii
        gfactor = ttmp * gfac[ir0] + ttmp1 * gfac[ir0-1];
//Let's interpolate lensing between two radii
        lensing = ttmp * transf_d[ir0] + ttmp1 * transf_d[ir0 - 1];
//Let's interpolate cosine between two radii
        
        x = sqrt(rnow);
  // Bcal = 1. + am / pow(x, 3.);
        Ccal = 1. - 3. / (x*x) + 2. * a / (x*x*x);
        Lcal = 1. / x * (x - x0 - 1.5 * a * log(x / x0) - 
            3. * pow(x1 - a,2.)/x1/(x1 - x2)/(x1 - x3) * log((x - x1)/(x0 - x1))-
            3. * pow(x2 - a,2.)/x2/(x2 - x1)/(x2 - x3) * log((x - x2)/(x0 - x2))-
            3. * pow(x3 - a,2.)/x3/(x3 - x1)/(x3 - x2) * log((x - x3)/(x0 - x3)));
        
        temp = Tnorm * pow(x, -1.5) * pow(arate, 0.25) * pow(m, -0.5) * sqrt(sqrt(Lcal / Ccal));
        temp *= gfactor;

        Fn_in += factor * lensing * SIGMAP / PI * pow(temp/K_KEVK/fcol,3.) / fcol;
        Fe_in += factor * lensing * SIGMA / PI * pow(temp/K_KEVK/fcol,4.);
        for (ie = 0; ie < NE; ie++)
            photar[ie] += factor * lensing / (expm1(ear[ie]/temp));
    }
        
    Fmax=0.;
    for (ie = 0; ie < NE; ie++) {
        photar[ie] *= pow(ear[ie],3);
        if (Fmax<=photar[ie]) {
            Fmax=photar[ie];
            Tbb=ear[ie];
        }
    }
    
    Tbb /= 2.82;
    fprintf(stdout, "Lamp-post sees blackbody of %f keV\n", Tbb);
    
//We find Lx and Fn_out from the temperature Tbb and nthcomp table values
//---------------------
//Let's compute normalisation of the nthcomp from given temperature, photon 
//index, Lx and Ec
    nthcomp_ne = (int) (log10(NTHCOMP_EMAX0/NTHCOMP_EMIN0)/NTHCOMP_DE)+1;
    if ((nthcomp_ear = (float *) malloc((nthcomp_ne+1) * sizeof(float))) == NULL) {
        return;
    }
    if ((nthcomp_far = (float *) malloc(nthcomp_ne * sizeof(float))) == NULL) {
        return;
    }
    if ((nthcomp_fer = (float *) malloc(nthcomp_ne * sizeof(float))) == NULL) {
        return;
    }
    
    for (ie = 0; ie <= nthcomp_ne; ie++)
        nthcomp_ear[ie] = (NTHCOMP_EMIN0 * 
            pow(NTHCOMP_EMAX0/NTHCOMP_EMIN0, ((float) ie) / nthcomp_ne));
    nthcomp_param[0] = gamma;    // Gamma
    nthcomp_param[1] = ecut; // kT_e
    nthcomp_param[2] = Tbb;    // kT_bb
    nthcomp_param[3] = 0.;   // int_type
    nthcomp_param[4] = 0.;   // Redshift
    nthcomp_param[5] = (float)(nthopt);   // output type
    option = 1;
    donthcomp1_(nthcomp_ear, &nthcomp_ne, nthcomp_param, &option, nthcomp_far, nthcomp_fer);
    //fprintf(stderr, "Nth_Te = %f\n", nthcomp_param[1]);
    
    Fn_out = Fe_out = 0.;
    for(ie=0; ie < nthcomp_ne; ie++){
        Fn_out += nthcomp_far[ie];
        Fe_out += nthcomp_far[ie]*(nthcomp_ear[ie]+nthcomp_ear[ie+1])/2.;
    }
    
  //photon and energy flux in the defined energy range
//  nthcomp_flux = 0.;
    nthcomp_lum = 0.;
    ie = (int) floor(nthcomp_ne * log10(lxemin/g_L/NTHCOMP_EMIN0)/log10(NTHCOMP_EMAX0/NTHCOMP_EMIN0));
    tmp = nthcomp_far[ie] * (nthcomp_ear[ie+1]-lxemin/g_L) / 
            (nthcomp_ear[ie+1]-nthcomp_ear[ie]);
//  nthcomp_flux = tmp;
        /* tmp in flux */
    nthcomp_lum = tmp * (lxemax/g_L+nthcomp_ear[ie+1])/2.;
    ie++;
    while (nthcomp_ear[ie+1] < lxemax/g_L && ie < nthcomp_ne) {
//    nthcomp_flux += nthcomp_far[ie];
        nthcomp_lum += nthcomp_far[ie]*(nthcomp_ear[ie]+nthcomp_ear[ie+1])/2.;
        ie++;
    }
    tmp = nthcomp_far[ie] * (lxemax/g_L-nthcomp_ear[ie]) / 
        (nthcomp_ear[ie+1]-nthcomp_ear[ie]);
//  nthcomp_flux += tmp;
    nthcomp_lum += tmp * (nthcomp_ear[ie]+lxemax/g_L)/2.; //i.e. in keV/s
    
// note that in the following we compute the total intrinsic X-ray flux
    /* convert from 2-10 keV to full band */
    if (bolometric == 0) 
        lx0 *= Fe_out / nthcomp_lum;
    //fprintf(stdout,"Lx(%.1f-%.1fkeV) = %e erg/s\n", lxemin, lxemax, lx0);
    //fprintf(stdout, "Bolometric to (%.1f-%.1f keV) flux ratio = %f\n", lxemin, lxemax, Fe_out / nthcomp_lum);
    /* Convert count flux from observer frame back to lamppost frame */
    /* Lx / Fe_out = area at lamppost frame */
    Fn_out *= lx0 / Fe_out;
    free((void *) nthcomp_ear); nthcomp_ear = NULL;
    free((void *) nthcomp_far); nthcomp_far = NULL;
//free((void *) photer1); photer1 = NULL;
//---------------------
    Fe_in *= Fn_out / Fn_in;
    //Ep = Lx / Fn_out;
//fprintf(stdout,"%E\t%E\n",Ebb,Ep);
    //photar_norm = 2. / pow(fcol*H_KEVS,3) / fcol / C_MS / C_MS * Fn_out / Fn_in / Lx;
    
    D = sqrt( 4. / PI  * g_L  * Fn_out / Fn_in / (1.-exp(-tau))) / Rg;
    //double dh = sqrt(g_L / PI  * Fn_out / Fn_in / (1.-exp(-tau))) / Rg;
    double dh = sqrt(g_L / PI  * Fn_out / Fn_in) / Rg;
    //D = sqrt( 4. / PI  * g_L  * Fn_out / Fn_in) / Rg;
    fprintf(stdout, "%f %E\n", h, dh);
    //fprintf(stderr, "g_L = %f", g_L); 
    /*
    fprintf(stdout,
        "h = %E, D = %E, Fe_in / Lx = %E, Tbb = %E keV\n\n",
        h, D, Fe_in/(KEV*Lx),Tbb);
        fprintf(fw_h,"%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\n",
        h, D, n_e, Fe_in/(KEV*Lx), Ebb, Tbb, Ep, Fn_out, 
        Lx/(BHmass*LEDD), g_L, transf_o, 
        Lx*KEV*g_L*g_L/arate/MSOLAR*YEAR/C_MS/C_MS/(1-Ut_rms),
        lnth_uv_sca*KEV*1e+7);
        */
    return;
}
