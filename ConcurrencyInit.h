#ifndef CONCURRENCYINIT_H_
#define CONCURRENCYINIT_H_

#include "nr3.h"

using namespace Eigen;

class TmatrixConcurrency
{
	private:
		vector<vector<double> >rfit;
		vector<vector<double> >pfit;
		vector<string> SA;
		int Lmax;
		Doub mw;
		Doub cp;
	public:
		TmatrixConcurrency(vector<vector<double> > rfit_, vector<vector<double> >pfit_,
						   vector<string> SA_, int Lmax_, Doub mw_, Doub cp_):
						   rfit(rfit_), pfit(pfit_), Lmax(Lmax_),
						   SA(SA_), mw(mw_), cp(cp_){};
		TmatrixConcurrency(){};
		~TmatrixConcurrency(){};
		void operator()(Doub k, vector<vector<double> >& Tll);
};
		

class ConcurrencyInit
{
	private:
		vector<double> kvector;
		int Nsc;
		vector<int> NSp;
		int Lmax;
		vector<vector<double> >Xcor ;
		vector<vector<double> >RBasis ;
		vector<vector<vector<vector<double> > > > tR;
		vector<int> Lv;
		vector<int> Mv;
		vector<int> Kmesh;
		double ConvergePrec;
	public:
		ConcurrencyInit(
				vector<double> kvector_, int Nsc_, vector<int> NSp_, int Lmax_,
				vector<vector<double> >Xcor_,
				vector<vector<vector<vector<double> > > > tR_,
				vector<int> Lv_, vector<int> Mv_, vector<int> Kmesh_,
				double ConvergePrec_, vector<vector<double> >RBasis_
				) : kvector(kvector_), Nsc(Nsc_), NSp(NSp_), Lmax(Lmax_),
				Xcor(Xcor_), tR(tR_), Lv(Lv_), Mv(Mv_), Kmesh(Kmesh_),
				ConvergePrec(ConvergePrec_), RBasis(RBasis_){};
		~ConcurrencyInit(){};
		int operator()( int k,
						 vector<vector<double> >& scatterR,
						 vector<vector<double> >& scatterI,
						 vector<vector<vector<complex<double> > > >& MKKRaux,
						 vector<vector<vector<complex<double> > > >& MSPOaux);
};

class ImpuritiesSPOConcurrency
{
	private:
		vector<double> kvector;
		vector<vector<double> >Xcor;
		vector<int> Lv;
		vector<int> Mv;
		vector<int> Kmesh;
		vector<vector<double> >RBasis ;
	public:
		ImpuritiesSPOConcurrency(vector<double> kvector_,
				vector<vector<double> >Xcor_,
				vector<int>Lv_, vector<int>Mv_,
				vector<int>Kmesh_,vector<vector<double> >RBasis_):
				kvector(kvector_),Xcor(Xcor_),Lv(Lv_),Mv(Mv_),
				Kmesh(Kmesh_),RBasis(RBasis_){};
		void operator()(
						vector<vector<vector<complex<double> > > >&MKKRaux_in,
						vector<vector<vector<complex<double> > > >&MSPOaux_in,
						vector<vector<double> >&MRI1R_out,
						vector<vector<double> >&MRI1I_out,
						vector<vector<double> >&MRI2R_out,
						vector<vector<double> >&MRI2I_out,
						vector<vector<double> >&MRI3R_out,
						vector<vector<double> >&MRI3I_out);

};

#endif
