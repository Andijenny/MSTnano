//#include<mutex>
//#include <inverse.h>
#ifdef MULTI_THREAD
	#include<thread>
#endif
//#include "Matrix.h"
#include <Dense>
#include <iomanip>
#include "KKRcoef.h"
#include "Tmatrix2.h"
#include "Geometry.h"
#include "ConcurrencyInit.h"
#include <complex>
#include "FFT3D.h"


void TmatrixConcurrency::operator()(Doub k, vector<vector<double> >& Tll)
{
//	static const double pi=4.0*atan(1);
	
	Tmatrix t1;
	for(int i=0;i<SA.size();i++){
		for(int l=0;l<=Lmax;l++){
			t1 = Tmatrix(l, k, rfit[i], pfit[i]);
			Tll[i][l] = t1.calc_tma(mw,cp).real();
		}
	}
};




int ConcurrencyInit::operator()(int kv,
								 vector<vector<double> >& scatterR, 
								 vector<vector<double> >& scatterI,
								 vector<vector<vector<complex<double> > > >& MKKRaux,
								 vector<vector<vector<complex<double> > > >& MSPOaux)
{
	
//	vector<vector<vector<complex<double> > > > MKKRaux(NumK,
//				  vector<vector<complex<double> > >(Nsc,
//					       vector<complex<double> >(Nsc, complex<double>(0.0, 0.0))));
//	
//	vector<vector<vector<complex<double> > > > MSPOaux(NumK,
//				  vector<vector<complex<double> > >(Nsc,
//					       vector<complex<double> >(Nsc, complex<double>(0.0, 0.0))));

	cout << "Enter ConcurrencyInit.operator().\n";
	cout << string(60,'-') << endl;
	
	int Nlm = tR[0][0][0].size();
	int Nspecies = tR.size();
  	if(Nspecies == 1)
	{
		for(int i=0;i<Nlm;i++)
		{
			for(int j=0;j<Nlm;j++)
			{
				scatterR[i][j] = tR[0][kv][i][j];
				scatterI[i][j] = 0.0;
			}
		}
		return 0;
	}
	


//	Map<Matrix<double,Dynamic,Dynamic> > Xcor_(&(Xcor[0][0]),
//											   Xcor.size(),Xcor[0].size());
	
//	MatrixXd Xcor_in = Xcor_.transpose();
	
	MatrixXd Xcor_in(Xcor[0].size(),Xcor.size());
	for(int i=0;i<Xcor.size();i++)
		for(int j=0;j<Xcor[0].size();j++)
			Xcor_in(j,i) = Xcor[i][j];

//	Map<Matrix<double,Dynamic,Dynamic> > BR(&(RBasis[0][0]),3,3);
	
	MatrixXd BR(3,3);
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			BR(i,j) = RBasis[i][j];

//	BR << 2.942, 0.000, 0.00,
//		  0.000, 5.095, 0.00,
//		  0.000, 0.000, 20.0;
	
	VectorXi NKmesh(3);
	for(int i=0;i<3;i++)
		NKmesh(i) = Kmesh[i];

//	NKmesh << 5,5,1;

	int NumK = NKmesh.prod();


	FFT3D f1(BR, NKmesh, Xcor_in, Lv, Mv, kvector);
	

	int num_threads = (int)(std::thread::hardware_concurrency());
	clock_t start_t, finish_t;
//	cout << "Calculate KKR constants in k space.\n"<<endl;
	start_t = clock();
	f1.fftKKR(MKKRaux, kv);
	finish_t = clock();
	cout << "3D Fourier transform method FFT3D.fftKKR costs ";
	cout << (double)(finish_t - start_t)/CLOCKS_PER_SEC/num_threads << " s.\n";
	cout << endl << endl;
//	cout << "Return from FFT3D .fftKKR\n"<<endl;

	MatrixXcd tmpTotalKKR = MatrixXcd::Zero(Nsc, Nsc);
	MatrixXcd tmpTotalT = MatrixXcd::Zero(Nsc, Nsc);
	MatrixXcd tmpScatter = MatrixXcd::Zero(Nsc, Nsc);

	int Nsp = Xcor[0].size(); 
	
	
	start_t = clock();
	for(int i=0;i<Nsp;i++)
		for(int j=0;j<Nsp;j++)
			for(int lm1=0;lm1<Nlm;lm1++)
				for(int lm2=0;lm2<Nlm;lm2++)
					if(i==j)
						tmpTotalT(Nlm*i+lm1,Nlm*j+lm2) = tR[i][kv][lm1][lm2];
	
	for(int kp=0;kp<NumK;kp++)
	{
		for(int i=0;i<Nsc;i++)
			for(int j=0;j<Nsc;j++)
				tmpTotalKKR(i,j) = MKKRaux[kp][i][j];
		tmpScatter = tmpTotalT.inverse()-tmpTotalKKR;
		tmpScatter = tmpScatter.inverse();
		for(int i=0;i<Nsc;i++)
			for(int j=0;j<Nsc;j++)
				MSPOaux[kp][i][j] = tmpScatter(i,j);
	}
	
	finish_t = clock();
	cout << "Matrice inverse costs ";
	cout << (double)(finish_t - start_t)/CLOCKS_PER_SEC << " s.\n";
	cout << endl << endl;

	start_t = clock();
	f1.ifft(MSPOaux,scatterR,scatterI);
	finish_t = clock();
	cout << "Inverse 3D Fourier transform method FFT3D.ifft costs ";
	cout << (double)(finish_t - start_t)/CLOCKS_PER_SEC/num_threads << " s.\n";
	cout << string(60,'-') << endl;
	cout << "Return from ConcurrencyInit.operator().\n"<<endl<<endl;
//	cout << "Return from .ifft().\n";
}

void ImpuritiesSPOConcurrency::operator()(
						vector<vector<vector<complex<double> > > >&MKKRaux_in,
						vector<vector<vector<complex<double> > > >&MSPOaux_in,
						vector<vector<double> >&MRI1R_out,
						vector<vector<double> >&MRI1I_out,
						vector<vector<double> >&MRI2R_out,
						vector<vector<double> >&MRI2I_out,
						vector<vector<double> >&MRI3R_out,
						vector<vector<double> >&MRI3I_out)
{

	MatrixXd Xcor_in(Xcor[0].size(),Xcor.size());
	for(int i=0;i<Xcor.size();i++)
		for(int j=0;j<Xcor[0].size();j++)
			Xcor_in(j,i) = Xcor[i][j];

	MatrixXd BR(3,3);
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			BR(i,j) = RBasis[i][j];
	
	VectorXi NKmesh(3);
	for(int i=0;i<3;i++)
		NKmesh(i) = Kmesh[i];


	int NumK = MKKRaux_in.size();
	int Nsc = MKKRaux_in[0].size();

	MatrixXcd tmpMKKR(Nsc,Nsc);
	MatrixXcd tmpMSPO(Nsc,Nsc);

	MatrixXcd tmpM1(Nsc,Nsc);
	MatrixXcd tmpM2(Nsc,Nsc);
	MatrixXcd tmpM3(Nsc,Nsc);

	vector<vector<vector<complex<double> > > > tmpM1Full(NumK,
					   vector<vector<complex<double> > >(Nsc,
					            vector<complex<double> >(Nsc,complex<double>(0.0,0.0))));
	
	vector<vector<vector<complex<double> > > > tmpM2Full(NumK,
					   vector<vector<complex<double> > >(Nsc,
					            vector<complex<double> >(Nsc,complex<double>(0.0,0.0))));

	vector<vector<vector<complex<double> > > > tmpM3Full(NumK,
					   vector<vector<complex<double> > >(Nsc,
					            vector<complex<double> >(Nsc,complex<double>(0.0,0.0))));
	for(int k=0;k<NumK;k++)
	{
		for(int j=0;j<Nsc;j++)
		{
			for(int i=0;i<Nsc;i++)
			{
				tmpMKKR(j,i) = MKKRaux_in[k][j][i];
				tmpMSPO(j,i) = MSPOaux_in[k][j][i];
			}
		}
		tmpM1 = tmpMSPO*tmpMKKR;
		tmpM2 = tmpMKKR*tmpMSPO;
		tmpM3 = tmpMKKR*tmpMSPO*tmpMKKR;
		for(int j=0;j<Nsc;j++)
		{
			for(int i=0;i<Nsc;i++)
			{
				tmpM1Full[k][j][i] = tmpM1(j,i);
				tmpM2Full[k][j][i] = tmpM2(j,i);
				tmpM3Full[k][j][i] = tmpM3(j,i);
			}
		}
	}
	
	FFT3D ift(BR, NKmesh, Xcor_in, Lv, Mv, kvector);
	
	ift.ifft(tmpM1Full, MRI1R_out, MRI1I_out);
	ift.ifft(tmpM2Full, MRI2R_out, MRI2I_out);
	ift.ifft(tmpM3Full, MRI3R_out, MRI3I_out);
}
