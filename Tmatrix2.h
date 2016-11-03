#ifndef _TMATRIX2_H_
#define _TMATRIX2_H_

#include<iostream>
#include<vector>
#include<complex>
using namespace std;


class Tmatrix
{
private:
	int LL;
	double kvalue;
	vector<double> radius;
	vector<double> potential;
public:
	Tmatrix(int LL_,  double  kvalue_, vector<double> radius_, vector<double> potential_) 
		: LL(LL_), kvalue(kvalue_), radius(radius_), potential(potential_){};
	Tmatrix(){};
	double DlNext(vector<double> Rl);
	void RlNext(vector<double>& Rl, double Dl);
	complex<double> calc_tma(double MixWeight, double ConvergPrecision);
	~Tmatrix(){};
};

//void const TmatrixPrint(const vector<vector<vector<double> > > Tll);

//void TmatrixRead(vector<vector<vector<double> > >& Tll);

#endif
