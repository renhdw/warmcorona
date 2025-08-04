//! \file rootsearch.h
//! modified version of rtbis() in sim5roots.c written by Michal. Now this function can handle
//! functions with additional parameter(s). To allow for different types and numbers of arguments,
//! I rewrote `rtbis()` into a function template. In this case additional parameters can
//! be passed within one class variable.
#ifndef _ROOTSEARCH_C
#define _ROOTSEARCH_C

const double MAX_STEPS = 1000;
template <typename T>
long rtbis(const double x1, const double x2, const double xacc, double (*fx)(double, const T&), const T & param, double & result)
//! Root finding.
//! Finds root of a function on an interval. Using bisection method, finds the root of a
//! function `fx` that is known to lie between `x1` and `x2`. The root, returned as `result`,
//! will be refined until its accuracy is +/- xacc.
//!
//! @param x1 left boundary of the interval where the root is searched for
//! @param x2 right boundary of the interval where the root is searched for
//! @param xacc accuracy
//! @param fx function
//! @param param parameter
//! @param result Returns 1 if OK and the root position in `result`, 0 if error.
{
	double dx, f, fmid, xmid, rtb;
	long j=0;

	fmid = (*fx)(x2, param);
	f    = (*fx)(x1, param);
	if ((f*fmid) >= 0.0) return(0);//error("rtbis: root is not bracketed");

	if (f < 0.0) {
		rtb = x1;
		dx  = x2-x1;
	}
	else {
		rtb = x2;
		dx  = x1-x2;
	}

	for (j=0; j<MAX_STEPS; j++) {
		dx = dx*0.5;
		xmid = rtb+dx;
		fmid = (*fx)(xmid, param);
		if (fmid <= 0.0) rtb = xmid;
		if ((fabs(dx) < xacc) || (fmid == 0.0)) break;
	}
	if (j >= MAX_STEPS) {
		std::cout << ("rtbis: too many steps") << std::endl;
		return(0);
	}

	result = rtb;
	return(1);
}

#endif
