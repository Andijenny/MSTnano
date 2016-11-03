#ifndef _GREENSFUNC_H_
#define _GREENSFUNC_H_
#include <iostream>
#include <vector>
#include <Dense>

//#include "nr3.h"
using namespace std;
using namespace Eigen;

int MinDis(vector<double> rr, vector<vector<double> > Xcor);

class ZFunc
{
	protected:
		int Lmax;
		vector<double> kvector;
		vector<int> Lv;
		vector<int> Mv;
		vector<vector<double> > Xcor;
		vector<vector<vector<vector<double> > > > tR;
	public:
		ZFunc(int Lmax_, vector<double> k_,vector<int> Lv_, vector<int> Mv_, 
				vector<vector<double> > Xcor_, 
				vector<vector<vector<vector<double> > > > tR_)
			: Lmax(Lmax_), kvector(k_), Lv(Lv_), Mv(Mv_), Xcor(Xcor_), tR(tR_){};
		ZFunc(){};
//		vector<double> operator()(vector<double> rr);
//		void operator()(vector<double> rr, vector<vector<complex> > & ZL);
		void operator()(int index, vector<double> rr, MatrixXcd & ZL);
		void JL(int index, vector<double> rr, double kv, VectorXcd & J0, VectorXcd & H0);
};


class GreenFunc : ZFunc
{
	private:
		
		vector<vector<vector<double> > > scatterR;
		vector<vector<vector<double> > > scatterI;
	
	public:

		GreenFunc(int Lmax_, vector<double> k_,vector<int> Lv_, vector<int> Mv_, 
				vector<vector<double> > Xcor_, 
				vector<vector<vector<vector<double> > > > tR_,
				vector<vector<vector<double> > > scatterR_,
				vector<vector<vector<double> > > scatterI_) :
				ZFunc(Lmax_, k_, Lv_, Mv_, Xcor_, tR_),
				scatterR(scatterR_), scatterI(scatterI_){};
		GreenFunc(){};
		~GreenFunc(){};
		void operator()(vector<double> rr1, vector<double> rr2, VectorXcd & GF);
	//	void ASALDOS(vector<double> IntegRadius, int row_index, vector<double> & RHO);//Atomic sphere approx.
};

#endif
