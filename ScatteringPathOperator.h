#ifndef SCATTERINGPATHOPERATOR_H_
#define SCATTERINGPATHOPERATOR_H_

#include <iostream>
#include <vector>

using namespace std;

class ScatteringPathOperator
{
	private:
		double ConvergePrec;
		vector<int> Lv;
		vector<int> Mv;
		vector<int> ParallelPar;
		vector<double> kvector;
		vector<vector<vector<double> > > Tll;
		vector<int> Kmesh;
		vector<vector<double> > RealSpaceBasis;
	
	public:
		ScatteringPathOperator(double ConvergePrec_,
							   vector<int> Lv_, vector<int> Mv_, 
							   vector<int> ParallelPar_, vector<double> kvector_,
							   vector<vector<vector<double> > > Tll_,
							   vector<vector<double> > RealSpaceBasis_,
							   vector<int> Kmesh_):
							   ConvergePrec(ConvergePrec_), 
							   Lv(Lv_), Mv(Mv_), 
							   ParallelPar(ParallelPar_), 
							   kvector(kvector_), Tll(Tll_),
							   RealSpaceBasis(RealSpaceBasis_),
							   Kmesh(Kmesh_){};
		ScatteringPathOperator(){};
		~ScatteringPathOperator(){};
		void operator()(vector<int> SL, vector<vector<double> > Xcor,
						vector<vector<vector<vector<double> > > >& tR,
						vector<vector<vector<double> > >& scatterR,
						vector<vector<vector<double> > >& scatterI,
						vector<vector<vector<vector<complex<double> > > > >& MKKRaux,
						vector<vector<vector<vector<complex<double> > > > >& MSPOaux);
		void ConstructTmatrix(vector<int> SL,
							  vector<vector<vector<vector<double> > > >& tR);
		void ImpuritiesSPO(vector<int> Lv_I, vector<int> Mv_I,
						   vector<int> SL_I, vector<vector<double> > Xcor_I,
						   vector<vector<vector<vector<double> > > >& tR_I,
						   vector<vector<vector<double> > >& scatterR_I,
						   vector<vector<vector<double> > >& scatterI_I);
		void ReNormInteractor(vector<int> Lv_I, vector<int> Mv_I,
							  vector<vector<double> >HostsXcor,
							  vector<vector<double> >ImpuritiesXcor,
							  vector<vector<vector<double> > >HostsScatterR,
							  vector<vector<vector<double> > >HostsScatterI,
							  vector<vector<vector<vector<complex<double> > > > >& MKKRaux,
							  vector<vector<vector<vector<complex<double> > > > >& MSPOaux,
							  vector<vector<vector<double> > >& ImpuritiesRIreal,
							  vector<vector<vector<double> > >& ImpuritiesRIimag);
};

#endif
