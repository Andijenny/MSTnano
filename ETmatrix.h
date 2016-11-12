#ifndef _ETMATRIX_H_
#define _ETMATRIX_H_
#include "nr3.h"

class ulscf
{
private:
	int LL;
	double kvalue;
	double dt;
	VecDoub radius;
	VecDoub pot;
	VecDoub ul;
	VecDoub besjr;
	VecDoub besnr;
public:
	ulscf(int LL_, double kvalue_, double dt_, 
		  VecDoub radius_, VecDoub pot_, VecDoub ul_,
		  VecDoub besjr_, VecDoub besnr_): 
		  LL(LL_), kvalue(kvalue_), dt(dt_),
		  radius(radius_), pot(pot_), ul(ul_),
		  besjr(besjr_), besnr(besnr_){};
	ulscf(){};
	void next(const double weight=1.0);
	Doub delta_new(double dt0, double MixWeight);
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


	
	
