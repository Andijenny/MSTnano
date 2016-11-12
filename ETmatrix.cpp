//#ifndef _TMATRIX_CPP_
//#define _TMATRIX_CPP_//for template class not supporting the seperation of defination and implement. 
#include "ETmatrix.h"
#include "qgaus.h"
#include "interp_1d.h"
#include "besselfrac.h"
#include <iomanip>
#include <fstream>
#include <ctime>

class func
{
	protected:
		int LL;
		double kvalue;
		VecDoub radius;
		VecDoub pot;
		VecDoub ul;
	public:
		func(int LL_, double kvalue_, VecDoub radius_, 
				VecDoub pot_, VecDoub ul_):
			LL(LL_), kvalue(kvalue_), radius(radius_), pot(pot_), ul(ul_){};
		~func(){}
	virtual double operator()(double rr) = 0;
};

struct func1:func
{
	func1(int LL_, double kvalue_, VecDoub radius_, VecDoub pot_, VecDoub ul_):
			func(LL_, kvalue_, radius_, pot_, ul_){};
	double operator()(double rr)
	{
		Poly_interp fp1 = Poly_interp(radius, pot, 4);
		Poly_interp fp2 = Poly_interp(radius, ul, 4);

		Bessel besf;
		
		return besf.sphbesj(LL, rr*kvalue)
			   *fp1.interp(rr)*fp2.interp(rr)*rr;
	}
};

struct func2:func
{
	func2(int LL_, double kvalue_, VecDoub radius_, VecDoub pot_, VecDoub ul_):
			func(LL_, kvalue_, radius_, pot_, ul_){};
	double operator()(double rr)
	{
		Poly_interp fp1 = Poly_interp(radius, pot, 4);
		Poly_interp fp2 = Poly_interp(radius, ul, 4);

		Bessel besf;
	
		return besf.sphbesy(LL, rr*kvalue)
				*fp1.interp(rr)*fp2.interp(rr)*rr;
	}
};

void ulscf::next(const double weight)
{
	Bessel besf;
	func1 funcc1(LL, kvalue, radius, pot, ul);
	func2 funcc2(LL, kvalue, radius, pot, ul);

	for(int i = 0; i < ul.size(); i++)
		ul[i] = ul[i]
				+weight*(
					radius[i]*(besjr[i]*cos(dt)
								+kvalue*besnr[i]
								 *qgaus(funcc1, radius[0], radius[i])
								+kvalue*besjr[i]
								 *qgaus(funcc2, radius[i], radius[radius.size()-1]))
					-ul[i]);
}

double ulscf::delta_new(double dt0, double MixWeight)
{
	func1 funcc1 = func1(LL, kvalue, radius, pot, ul);
	
	double t1 = sin(dt0)+ MixWeight*(-sin(dt0) 
							- kvalue*qgaus(funcc1, radius[0], radius[radius.size()-1]));
	if(abs(t1) > 1.0) 
	{
		cout << "sin(dt): " << t1 << endl;
		cout << "Warning, phaseshift abnormal!" << endl;
		srand(time(0));
		t1=1.0-0.5*(rand()%11/10.0);
	}
	dt = asin(t1);
	return dt;
}


complex<double> Tmatrix::calc_tma(double MixWeight, double ConvergPrecision)
{
	
	VecDoub u0(radius.size(), 0.0);
	VecDoub j0(radius.size(), 0.0);
	VecDoub y0(radius.size(), 0.0);
	
	Bessel besj;
	Bessel besy;
	Doub tLL = 0.0;
	for(int i = 0; i < radius.size(); i++)
	{
		j0[i] = besj.sphbesj(LL, radius[i]*kvalue);
		y0[i] = besy.sphbesy(LL, radius[i]*kvalue);
		u0[i] = radius[i]*cos(0.0)*j0[i];
	}
	VecDoub CopyRadius(radius.size(), 0.0);
	VecDoub CopyPotential(potential.size(), 0.0);
	for(int i = 0; i < radius.size(); i++)
	{
		CopyRadius[i] = radius[i];
		CopyPotential[i] = potential[i];
	}
	
    ulscf ul = ulscf(LL, kvalue, 0.0, CopyRadius, CopyPotential, u0, j0, y0);	
	
	Doub t1 = 10.0;
	double t0 = ul.delta_new(0.0, MixWeight);
	
	int cc = 0;
	while(abs(t0 - t1) > abs(ConvergPrecision*t0))
	{
		cc++;
		t0 = t1;
		ul.next(MixWeight);
		t1 = ul.delta_new(t0, MixWeight);
	}
	
	tLL = sin(t1)/(-kvalue);
	cout<<"For angular quantum number l= "<< LL << endl;
	cout << "converged at the "<< cc <<"th step\n";
	cout << string('*',30) << endl;
	cout<<" kvector "<<"		tLL'		  "<<"	PhastShift		"<< endl; 
	cout << setw(10)<<setprecision(6)<<kvalue;
	cout <<"		"<< tLL <<"		"<< t1 <<endl;
	cout << string('*',30) << endl;
	return complex<double>(tLL, 0.0);
}


//#endif	
