#ifndef QUAD3D_H_
#define QUAD3D_H_

#include "nr3.h"
#include "qgaus.h"
#include "krig.h"
//#include "romberg.h"


struct NRf3 {
	Doub xsav,ysav;
	MatDoub FBZ;
	MatDoub x_coarse;
	VecDoub y_coarse;
	Doub (*func3d)(const MatDoub&, const VecDoub&, const Doub, const Doub, const Doub);
	Doub operator()(const Doub z)
	{
		return func3d(x_coarse,y_coarse,xsav,ysav,z);
	}
};

struct NRf2 {
	NRf3 f3;
	Doub (*z1)(MatDoub, Doub, Doub);
	Doub (*z2)(MatDoub, Doub, Doub);
	NRf2(Doub zz1(MatDoub, Doub, Doub), Doub zz2(MatDoub, Doub, Doub)) : z1(zz1), z2(zz2) {}
	Doub operator()(const Doub y)
	{
		f3.ysav=y;
		return qgaus(f3,z1(f3.FBZ,f3.xsav,y),z2(f3.FBZ,f3.xsav,y));
	//	return qromb(f3,z1(f3.xsav,y),z2(f3.xsav,y),1.0e-6);
//		return qtrap(f3,z1(f3.xsav,y),z2(f3.xsav,y),1.0e-6);
	}
};

struct NRf1 {
	Doub (*y1)(MatDoub, Doub);
	Doub (*y2)(MatDoub, Doub);
	NRf2 f2;
	NRf1(Doub yy1(MatDoub, Doub), Doub yy2(MatDoub, Doub), Doub z1(MatDoub, Doub, Doub),
		Doub z2(MatDoub, Doub, Doub)) : y1(yy1),y2(yy2), f2(z1,z2) {}
	Doub operator()(const Doub x)
	{
		f2.f3.xsav=x;
		return qgaus(f2,y1(f2.f3.FBZ,x),y2(f2.f3.FBZ,x));
	//	return qromb(f2,y1(x),y2(x),1.0e-6);
//		return qtrap(f2,y1(x),y2(x),1.0e-6);
	}
};

template <class T>
Doub quad3d(T &func,MatDoub & x_in, VecDoub& y_in, MatDoub& BZ, 
			const Doub x1, const Doub x2, 
		    Doub y1(MatDoub, Doub), Doub y2(MatDoub, Doub),
			Doub z1(MatDoub, Doub, Doub), Doub z2(MatDoub, Doub, Doub))
{
	NRf1 f1(y1,y2,z1,z2);
	f1.f2.f3.FBZ = BZ;
	f1.f2.f3.func3d = func;
	f1.f2.f3.x_coarse = x_in;
	f1.f2.f3.y_coarse = y_in;
	return qgaus(f1,x1,x2);
//	return qromb(f1,x1,x2,1.0e-5);
//	return qtrap(f1,x1,x2,1.0e-5);
}

#endif
