//! \file quadroots.h
//! Defines the `quadroots` object
#ifndef _QUADROOTS_H
#define _QUADROOTS_H
#include <complex>
#include <vector>
#include <iostream>

/*!
This object searches for roots of \f$R(r)=0\f$. The definition of \f$R(r)\f$
can be found in documents of \ref kerr_rintegrand().
*/
class quadroots {
public:
    //! number of real roots
    int nrr; 
    //! whether the roots are sorted
    bool issorted;
    //! roots
    std::vector<std::complex<double>> roots;
    //! constructor
    quadroots(const double & a, const double & l, const double & q);
    //! sorting function
    void sort();
};
//! Overloaded operator <<.
std::ostream & operator << (std::ostream &, const quadroots &);
#endif
