//! \file quadroots.cpp
//! Implementing `quadroots` objects
#include <iostream>
#include "quadroots.h"

/*!
Find out the roots. From `sim5`.
@param a spin
@param l angular momentum about z-axis
@param q Carter's constant
*/
quadroots::quadroots(const double & a, const double & l, const double & q) {
    double C, D, E, F, A, B, X, z, Z;
    std::complex<double> r1, r2, r3, r4, temp1, temp2;
    C = (a - l) * (a - l) + q;
    D = 2. * (q + l * l - a * a) / 3.;
    E = 9. * D * D / 4. - 12. * a * a * q;
    F = -27. * D * D * D / 4. - 108. * a * a * q * D + 108. * C * C;
    X = F * F - 4. * E * E * E;

    if (X >= 0) {
        A = (F>std::sqrt(X)?+1:-1)*1./3.*std::pow(std::fabs(F-std::sqrt(X))/2.,1./3.)+ 
            (F>-std::sqrt(X)?+1:-1)*1./3.*std::pow(std::fabs(F+std::sqrt(X))/2.,1./3.);
    } else {
        Z = std::sqrt(std::pow(F/54.,2) + std::pow(std::sqrt(-X)/54.,2));
        z = atan2(std::sqrt(-X)/54.,F/54.);
        A = std::pow(Z,1./3.)*2.*std::cos(z/3.);
    }

    B = std::sqrt(A + D);
    temp1 = -A + 2. * D - 4. * C / B;
    temp2 = -A + 2. * D + 4. * C / B;

    temp1 = std::sqrt(temp1);
    temp2 = std::sqrt(temp2);

    roots.push_back(B / 2. + .5 * temp1);
    roots.push_back(B / 2. - .5 * temp1);
    roots.push_back(-0.5 * B + .5 * temp2);
    roots.push_back(-0.5 * B - .5 * temp2);

    issorted = false;
}

void quadroots::sort() {
  //double complex rr1[N], rr2[N], tmp;
  nrr = 0;
  std::vector<std::complex<double>> rr2(4);
  std::complex<double> tmp;

  for (int i=0; i<4; ++i) {
      if (roots[i].imag() == 0.) {
          rr2[nrr] = roots[i];
          nrr += 1;
      }
  }

  int k = nrr;
  for (int i=0; i<4; i++) {
      if (roots[i].imag() != 0.) {
          rr2[k] = roots[i];
          k += 1;
      }
  }

  for (int i=0; i<nrr; i++) {
      for (int j=0; j<(nrr-i); j++) {
        if (rr2[i+j].real() > rr2[i].real()) {
            tmp=rr2[i+j];
            rr2[i+j]=rr2[i];
            rr2[i]=tmp;
	    }
      }
    }

    roots = rr2;
    issorted = true;
}

/*!
@param cout ostream object. 
@param roots \ref quadroots object.
 */
std::ostream & operator <<(std::ostream & cout, const quadroots & roots) {
  cout << "nrr = " << roots.nrr << std::endl;
  if (roots.nrr == 4) {
     cout << roots.roots[0].real() << " " << roots.roots[1].real() << " " << roots.roots[2].real() << " " << roots.roots[3].real() << std::endl;
  }
  if (roots.nrr == 2) {
     cout << roots.nrr << " " << roots.roots[0].real() << " " << roots.roots[1].real() << std::endl;
     cout << "One pair of conjugate roots: " << roots.roots[2].real() << " +/- " << roots.roots[2].imag() << std::endl;
  }
     
  if (roots.nrr == 0) {
     cout << "1st pair of conjugate roots: " << roots.roots[0].real() << " +/- " << roots.roots[0].imag() << std::endl;
     cout << "2nd pair of conjugate roots: " << roots.roots[2].real() << " +/- " << roots.roots[2].imag() << std::endl;
  }
  return cout;
}
