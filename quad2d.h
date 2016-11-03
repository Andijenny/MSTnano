#ifndef QUAD2D_H_
#define QUAD2D_H_

#include "nr3.h"
#include "qgaus.h"
//#include "krig.h"
#include "interp_2d.h"
//#include "romberg.h"


struct NRf2 {
	Doub xsav;
	VecDoub FBZ;
//	Bilin_interp myf1;
	VecDoub x1_coarse;
	VecDoub x2_coarse;
	MatDoub y_coarse;
	Doub (*func2d)(const VecDoub&, const VecDoub&, const MatDoub&, const Doub, const Doub);
//	Doub (*func2d)(Bilin_interp, const Doub, const Doub);
	Doub operator()(const Doub y)
	{
		return func2d(x1_coarse, x2_coarse, y_coarse, xsav,y);
//		return func2d(myf1, xsav,y);
	}
};


struct NRf1 {
//	Doub (*y1)(MatDoub, Doub);
//	Doub (*y2)(MatDoub, Doub);
	Doub (*y1)(VecDoub, Doub);
	Doub (*y2)(VecDoub, Doub);
	NRf2 f2;
//	NRf1(Doub yy1(MatDoub, Doub), Doub yy2(MatDoub, Doub)) : y1(yy1),y2(yy2) {}
	NRf1(Doub yy1(VecDoub, Doub), Doub yy2(VecDoub, Doub)) : y1(yy1),y2(yy2) {}
	Doub operator()(const Doub x)
	{
		f2.xsav=x;
		return qgaus(f2,y1(f2.FBZ,x),y2(f2.FBZ,x));
	//	return qromb(f2,y1(x),y2(x),1.0e-6);
//		return qtrap(f2,y1(x),y2(x),1.0e-6);
	}
};

template <class T>
Doub quad2d(T &func,VecDoub & x1_in,VecDoub & x2_in, MatDoub& y_in, VecDoub& BZ, 
//Doub quad2d(T &func,Bilin_interp myf1, VecDoub& BZ,
			const Doub x1, const Doub x2, 
		    Doub y1(VecDoub, Doub), Doub y2(VecDoub, Doub))
{
	NRf1 f1(y1,y2);
	f1.f2.FBZ = BZ;
	f1.f2.func2d = func;
	f1.f2.x1_coarse = x1_in;
	f1.f2.x2_coarse = x2_in;
	f1.f2.y_coarse = y_in;
	return qgaus(f1,x1,x2);
//	return qromb(f1,x1,x2,1.0e-5);
//	return qtrap(f1,x1,x2,1.0e-5);
}

#endif
