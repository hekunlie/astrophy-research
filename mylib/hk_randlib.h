#ifndef __HK_RANDLIB__
#define __HK_RANDLIB__

#include <cmath>

template<class T>
inline const T SQR(const T a) {return a*a;}

struct Normaldev : Ran 
{
	double mu,sig;
	Normaldev(double mmu, double ssig, long seed)
	: Ran(seed), mu(mmu), sig(ssig){}
	double dev() {
		double u,v,x,y,q;
		do {
			u = double();
			v = 1.7156*(double()-0.5);
			x = u - 0.449871;
			y = fabs(v) + 0.386595;
			q = SQR(x) + y*(0.19600*y-0.25472*x);
		} while (q > 0.27597 && (q > 0.27846 || SQR(v) > -4.*log(u)*SQR(u)));
		return mu + sig*v/u;
	}
};

#endif