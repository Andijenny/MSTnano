#ifndef _KKRCOEF_H_
#define _KKRCOEF_H_
#include "nr3.h"
//#define Int int
//#define Doub double
//#define Complex complex<double>

class KKRcoef{
//KKR conefficients
// from 'Green's functions for ordered and disordered systems' Pp375
// the coefficients give rise to the real-space structure constants of the lattice,
// and belong to the non-diagonal elements of the scatterring path operator matrix. 
public:
	KKRcoef(const Int ll1,const Int mm1,const Int ll2,const Int mm2, vector<int> Lv_, vector<int> Mv_):
		l1(ll1),m1(mm1),l2(ll2),m2(mm2), Lv(Lv_), Mv(Mv_) {};
	KKRcoef(){};
	Complex StuctConst(Doub & k, VecDoub & R);
	Complex StuctConst(Doub & k, vector<double> & R);
	Complex TranslationConst(Doub & k, vector<double> & R);
	~KKRcoef(){};
private:
	Int l1;
	Int m1;
	Int l2;
	Int m2;
	vector<int> Lv;
	vector<int> Mv;
};
#endif
