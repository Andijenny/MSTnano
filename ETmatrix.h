#ifndef _ETMATRIX_H_
#define _ETMATRIX_H_
#include "nr3.h"


class ulscf
{
private:
	Int l;
	Doub k;
	Doub dt;
	VecDoub r;
	VecDoub ur;//previous ur
	VecDoub qd1;//quadrature term 1
	VecDoub qd2;//quadrature term 2
	Doub qd3;//quadrature term 3
	VecDoub pot;
	VecDoub besjr;
	VecDoub besnr;
public:
	ulscf(Int ll,  VecDoub & rr, VecDoub urr,  VecDoub pott, 
		   VecDoub qdd1,  VecDoub qdd2,  Doub qdd3,
		   VecDoub jrr, VecDoub nrr, Doub kk, Doub dtt)
		: l(ll), r(rr), ur(urr), pot(pott),qd1(qdd1), qd2(qdd2), qd3(qdd3), 
		   besjr(jrr), besnr(nrr), k(kk), dt(dtt){};
	ulscf(){};
	Doub quadd();
	void next(const Doub & weight=1.0);
	VecDoub showu(){return ur;};
	Doub delta_new();
	~ulscf(){};
};


class Tmatrix
{
private:
	int LL;//angular momentum quantum number
	double kvalue;// k-vector value
	vector<double> radius;
	vector<double> potential; 
public:
	Tmatrix(int LL_,  double  kvalue_, vector<double> radius_, vector<double> potential_) 
		: LL(LL_), kvalue(kvalue_), radius(radius_), potential(potential_){};
	Tmatrix(){};
	complex<double> calc_tma(double MixWeight, double ConvergPrecision);
	~Tmatrix(){};
};


//#include "Tmatrix.cpp"// have to for template class
#endif


	
	
