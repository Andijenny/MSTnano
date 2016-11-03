#ifndef FFT3D_H_
#define FFT3D_H_

#include <iostream>
#include <Dense>
#include <vector>
#include <complex>
#include "nr3.h"

using namespace std;

using namespace Eigen;


class FFT3D
{
	private:
//		static const double pi;
		MatrixXd VR;
		VectorXi Nk;
		MatrixXd Xcor;
		vector<int> Lv;
		vector<int> Mv;
		vector<double> kvector;
		double volumeR;
		MatrixXd Kpoints;
		MatrixXd KBasis;
	public:
		FFT3D(MatrixXd BasisInRealSpace_, VectorXi KP_, MatrixXd Xcor_, 
			  vector<int> Lv_, vector<int> Mv_, vector<double> kvector_);
//		FFT3D(vector<vector<double> > BasisInRealSpace_, vector<int> KP_, 
//			  vector<vector<double> > Xcor_, vector<int> Lv_, vector<int> Mv_, 
//			  vector<double> kvector_);
		FFT3D(){};
		void GenBasisInBZ(MatrixXd & KBasis);
		void KpointSampling(MatrixXd & Kpoints);
//		void fftKKR(MatrixXcd  Func_in, VectorXcd & Func_fft, int kv);//{cout << "Vancy here\n";};
		void fftKKR(vector<vector<vector<complex<double> > > >& MKKRaux, int kv);
		void ifft(vector<vector<vector<complex<double> > > > MSPOaux,
				 vector<vector<double> >& scatterR,
				 vector<vector<double> >& scatterI);
//		void fft(vector<vector<complex<double> > > Func_in, 
//			 vector<complex<double> >& Func_fft, int kv);
//		void WeightM(MatrixXd & wt);
//		void ifft(VectorXcd & Func_in, MatrixXcd & Func_ifft, int Nlm);
//		void ifft(vector<complex<double> >& Func_in, 
//			vector<vector<complex<double> > >& Func_ifft, int Nlm);
//		void ifft2(VectorXcd & Func_in, MatrixXcd & Func_ifft, int Nlm);
		~FFT3D(){};
};


class fftKKR_divideAndConquer
{
	private:
		MatrixXd VR;
		MatrixXd Xcor;
		vector<int> Lv;
		vector<int> Mv;
		MatrixXd Kpoints;
		vector<double> kvector;
		int Naux;
	public:
		fftKKR_divideAndConquer(MatrixXd VR_, MatrixXd Xcor_, vector<int> Lv_,
								vector<int> Mv_, MatrixXd Kpoints_, 
								vector<double> kvector_, int Naux_):
								VR(VR_), Xcor(Xcor_), Lv(Lv_), Mv(Mv_), 
								Kpoints(Kpoints_), kvector(kvector_), Naux(Naux_){};
		void operator()(int kp, int i, int j, int k,
						vector<vector<complex<double> > >& MKKRaux, int kv,
						int MrowStart, int MrowEnd);
		~fftKKR_divideAndConquer(){};
};

class ifft_divideAndConquer
{
	private:
		vector<vector<vector<complex<double> > > > MSPOaux;
		MatrixXd Kpoints;
		MatrixXd KBasis;
		VectorXi Nk;
	public:
		ifft_divideAndConquer(vector<vector<vector<complex<double> > > > MSPOaux_,
							  MatrixXd Kpoints_, MatrixXd KBasis_, VectorXi Nk_):
							  MSPOaux(MSPOaux_), Kpoints(Kpoints_), 
							  KBasis(KBasis_), Nk(Nk_){};
		void operator()(vector<vector<double> >& scatterR,
						vector<vector<double> >& scatterI,
						int MrowStart, int MrowEnd);
};


#endif
