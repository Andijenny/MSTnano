#include<algorithm>

#ifdef MULTI_THREAD

	#include<thread>

#endif

#include <iomanip>
#include <fstream>
#include <Dense>
//#include<mutex>
#include "GreensFunc.h"
#include "ConcurrencyInit.h"
#include "Geometry.h"
//#include "ParInitial.h"
//#include "ParaPrint.h"
#include "KKRcoef.h"
#include "ScatteringPathOperator.h"
#include "io.h"

//std::mutex some_mutex;
void ScatteringPathOperator::ConstructTmatrix(vector<int> SL, 
								vector<vector<vector<vector<double> > > >& tR)
{
	int Nspecies = tR.size();
	int Nk = tR[0].size();
	int Nlm = tR[0][0].size();
	int Lmax = (int)sqrt(Nlm)-1;
//	int Lmax = (int)Tll[0][0].size()-1;
	
	for(int i=0;i<Nspecies;i++)
	{
		for(int k=0;k<Nk;k++)
		{
			for(int lm1=0;lm1<Nlm;lm1++)
			{
				for(int lm2=0;lm2<Nlm;lm2++)
				{
					if((Lv[lm1]==Lv[lm2]) && (Mv[lm1]==Mv[lm2]))
					{
						tR[i][k][lm1][lm2] = Tll[k][SL[i]][Lv[lm2]];
					}
					else 
					{
						tR[i][k][lm1][lm2] = 0.0;
					}
				}
			}
		}
	}
}

void ScatteringPathOperator::operator()(
		vector<int> SL, vector<vector<double> > Xcor,
		vector<vector<vector<vector<double> > > >& tR,
		vector<vector<vector<double> > >& scatterR,
		vector<vector<vector<double> > >& scatterI,
		vector<vector<vector<vector<complex<double> > > > >& MKKRaux,
		vector<vector<vector<vector<complex<double> > > > >& MSPOaux)
{
	int Nspecies = tR.size();
	int Nk = tR[0].size();
	int Nlm = tR[0][0].size();
	int Lmax = (int)sqrt(Nlm)-1;
	int Nsc = Nspecies*Nlm;

	if((scatterR.size()!=Nk) || (scatterI.size()!=Nk)
	   ||(MKKRaux.size()!=Nk) || (MSPOaux.size()!=Nk))
		throw("Wrong input matrices size.");
    
	ConstructTmatrix(SL, tR);

	vector<int> NSp(Nsc);
	vector<int> NL(Nsc);
	vector<int> NM(Nsc);
	int c1=-1;
	for(int i=0;i<Nspecies;i++){
		for(int lm=0;lm<Nlm;lm++){
			c1++;
			NSp[c1] = i;
			NL[c1] = Lv[lm];
			NM[c1] = Mv[lm];
		}
	}
	
	ConcurrencyInit SPO_multi(kvector, Nsc, NSp, Lmax, Xcor,
							  tR, Lv, Mv, Kmesh, ConvergePrec, RealSpaceBasis);

	clock_t start_t, finish_t;
	double total_t;
	start_t = clock();

	#ifdef MULTI_THREAD

		int num_threads = ParallelPar[0];
		int Ncycle = ParallelPar[1];
		int Nexceed = ParallelPar[2];
	

		for(int j=0;j<Ncycle;j++)
		{
	//		std::vector<std::thread> threads_s(num_threads-1);
			std::vector<std::thread> threads_s(num_threads);
		//	for(int i=0;i<(num_threads-1);i++)
			for(int i=0;i<(num_threads);i++)
			{
	//			std::lock_guard<std::mutex>guard(some_mutex);
//				threads_s[i] = std::thread(SPO_multi,((i+1)+j*num_threads),
				threads_s[i] = std::thread(SPO_multi,((i)+j*num_threads),
									   std::ref(scatterR[(i)+j*num_threads]),
						               std::ref(scatterI[(i)+j*num_threads]),
									   std::ref(MKKRaux[(i)+j*num_threads]),
									   std::ref(MSPOaux[(i)+j*num_threads]));
			}
			
//			SPO_multi(j*num_threads, scatterR[j*num_threads], 
//					  scatterI[j*num_threads], MKKRaux[j*num_threads], 
//					  MSPOaux[j*num_threads]);
			
//			for(int i=0;i<(num_threads-1);i++) threads_s[i].join();
			for(int i=0;i<(num_threads);i++) threads_s[i].join();
			
		}
		
	
		std::vector<std::thread> threads_ds(Nexceed);
		for(int i=0;i<(Nexceed);i++)
			threads_ds[i]=std::thread(SPO_multi, (i+Ncycle*num_threads),
								   std::ref(scatterR[Ncycle*num_threads+i]),
						           std::ref(scatterI[Ncycle*num_threads+i]),
								   std::ref(MKKRaux[Ncycle*num_threads+i]),
								   std::ref(MSPOaux[Ncycle*num_threads+i]));
		for(int i=0;i<(Nexceed);i++) threads_ds[i].join();
	
	#else
		for(int k=0;k<Nk;k++)
			SPO_multi(k, scatterR[k], scatterI[k], MKKRaux[k], MSPOaux[k]);

	#endif

	finish_t = clock();
	total_t = (double)(finish_t - start_t)/CLOCKS_PER_SEC;

	cout <<"SPOMatrix at "<< pow(hbar*kvector[0],2.0)/2.0/me/ee <<" eV converged, ";
	cout <<"took "<< total_t << " s." << endl;
	cout << endl;
	
}

void ScatteringPathOperator::ImpuritiesSPO(vector<int> Lv_I, vector<int> Mv_I,
										   vector<int> SL_I, vector<vector<double> > Xcor_I,
										   vector<vector<vector<vector<double> > > >& tR_I,
										   vector<vector<vector<double> > >& scatterR_I,
										   vector<vector<vector<double> > >& scatterI_I)
{
	
	int Nspecies = tR_I.size();
	int Nk = tR_I[0].size(); 
	int Nlm = tR_I[0][0].size();
	
	
	ConstructTmatrix(SL_I, tR_I);
	
	if(Nspecies == 1)
	{
		for(int k=0;k<Nk;k++)
		{
			for(int i=0;i<Nlm;i++)
			{
				for(int j=0;j<Nlm;j++)
				{
					scatterR_I[k][i][j] = tR_I[0][k][i][j];
					scatterI_I[k][i][j] = 0.0;
				}
			}
		}
		return;
	}
	
	KKRcoef G1;
  	Complex ctmp;
  	
	vector<double>Xij(3,0.0);
  	vector<double>Rij(3,0.0);
	
	for(int k=0;k<Nk;k++)
	{
		for(int i=0;i<Nspecies;i++)
		{
			for(int j=0;j<Nspecies;j++)
			{
				if(i == j)
				{
					for(int lm1=0;lm1<Nlm;lm1++)
					{
						for(int lm2=0;lm2<Nlm;lm2++)
						{
							scatterR_I[k][Nlm*i+lm1][Nlm*j+lm2] = tR_I[i][k][lm1][lm2];
							scatterI_I[k][Nlm*i+lm1][Nlm*j+lm2] = 0.0;
						}
					}
				}
				else
				{
					for(int nr=0;nr<3;nr++)
						Xij[nr]=(Xcor_I[nr][j]-Xcor_I[nr][i]);
					Cart2Sph(Xij,Rij);	
					for(int lm1=0;lm1<Nlm;lm1++)
					{
						for(int lm2=0;lm2<Nlm;lm2++)
						{
							G1 = KKRcoef(Lv_I[lm1],Mv_I[lm1],Lv_I[lm2],Mv_I[lm2],Lv_I,Mv_I);
							ctmp = G1.StuctConst(kvector[k],Rij);	
							scatterR_I[k][i*Nlm+lm1][j*Nlm+lm2] = ctmp.real(); 
							scatterI_I[k][i*Nlm+lm1][j*Nlm+lm2] = ctmp.imag(); 
					//		if ((i == 257 && j == 374) || (i == 374 && j == 257))
					//		{
					//			cout << "Rij: "<< Rij[0] <<"	"<<Rij[1]<<"	"<<Rij[2]<< endl;
					//			cout << "ctmp: "<< ctmp << endl;
					//		}
						}
					}
				}
			}
		}
	}
//	ScatterMatrixPrint("ImscatterR", "ImscatterI", scatterR_I, scatterI_I);
}


void ScatteringPathOperator::ReNormInteractor(vector<int> Lv_I, vector<int> Mv_I,
											  vector<vector<double> >HostsXcor,
											  vector<vector<double> >ImpuritiesXcor,
											  vector<vector<vector<double> > >HostsScatterR,
											  vector<vector<vector<double> > >HostsScatterI,
											  vector<vector<vector<vector<complex<double> > > > >& MKKRaux,
											  vector<vector<vector<vector<complex<double> > > > >& MSPOaux,
											  vector<vector<vector<double> > >& ImpuritiesRIreal,
											  vector<vector<vector<double> > >& ImpuritiesRIimag)
{
	int Nk = kvector.size();
	int Nsc = HostsScatterR[0].size();
	int Nlm = Nsc/HostsXcor[0].size();
	int Nsc_I = ImpuritiesRIreal[0].size();
	int Nlm_I = Nsc_I/ImpuritiesXcor[0].size();

	vector<vector<vector<double> > > MRI1R(Nk,
	        	   vector<vector<double> >(Nsc, vector<double>(Nsc,0.0)));
	vector<vector<vector<double> > > MRI1I(Nk,
	        	   vector<vector<double> >(Nsc, vector<double>(Nsc,0.0)));
	vector<vector<vector<double> > > MRI2R(Nk,
	        	   vector<vector<double> >(Nsc, vector<double>(Nsc,0.0)));
	vector<vector<vector<double> > > MRI2I(Nk,
	        	   vector<vector<double> >(Nsc, vector<double>(Nsc,0.0)));
	vector<vector<vector<double> > > MRI3R(Nk,
	        	   vector<vector<double> >(Nsc, vector<double>(Nsc,0.0)));
	vector<vector<vector<double> > > MRI3I(Nk,
	        	   vector<vector<double> >(Nsc, vector<double>(Nsc,0.0)));

	
//	ScatterMatrixRead("MRI1R","MRI1I",1,MRI1R,MRI1I);	
//	ScatterMatrixRead("MRI2R","MRI2I",1,MRI2R,MRI2I);	
//	ScatterMatrixRead("MRI3R","MRI3I",1,MRI3R,MRI3I);	

	
	clock_t start_t, finish_t;
	double total_t;
	start_t = clock();

	ImpuritiesSPOConcurrency RI(kvector, HostsXcor, Lv, Mv, Kmesh, RealSpaceBasis);
	#ifdef MULTI_THREAD
	
		int num_threads = ParallelPar[0];
		int Ncycle = ParallelPar[1];
		int Nexceed = ParallelPar[2];
		
		for(int j=0;j<Ncycle;j++)
		{
	//		std::vector<std::thread> threads_s(num_threads-1);
			std::vector<std::thread> threads_s(num_threads);
//			for(int i=0;i<(num_threads-1);i++)
			for(int i=0;i<(num_threads);i++)
			{
				threads_s[i] = std::thread(RI,
						//			   std::ref(MKKRaux[i+1+j*num_threads]),
									   std::ref(MKKRaux[i+j*num_threads]),
						               std::ref(MSPOaux[i+j*num_threads]),
									   std::ref(MRI1R[i+j*num_threads]),
									   std::ref(MRI1I[i+j*num_threads]),
									   std::ref(MRI2R[i+j*num_threads]),
									   std::ref(MRI2I[i+j*num_threads]),
									   std::ref(MRI3R[i+j*num_threads]),
									   std::ref(MRI3I[i+j*num_threads]));
			}
			
			//master thread:	
//			RI(MKKRaux[j*num_threads], MSPOaux[j*num_threads], 
//			   MRI1R[j*num_threads], MRI1I[j*num_threads], 
//			   MRI2R[j*num_threads], MRI2I[j*num_threads], 
//			   MRI3R[j*num_threads], MRI3I[j*num_threads]);
			
//			for(int i=0;i<(num_threads-1);i++) threads_s[i].join();
			for(int i=0;i<(num_threads);i++) threads_s[i].join();
		}
	
	//		ScatterMatrixPrint("MRI1R", "MRI1I", MRI1R, MRI1I);
	//		ScatterMatrixPrint("MRI2R", "MRI2I", MRI2R, MRI2I);
	//		ScatterMatrixPrint("MRI3R", "MRI3I", MRI3R, MRI3I);
			
		std::vector<std::thread> threads_ds(Nexceed);
		for(int i=0;i<(Nexceed);i++)
			threads_ds[i] = std::thread(RI,
								   std::ref(MKKRaux[Ncycle*num_threads+i]),
					               std::ref(MSPOaux[Ncycle*num_threads+i]),
								   std::ref(MRI1R[Ncycle*num_threads+i]),
								   std::ref(MRI1I[Ncycle*num_threads+i]),
								   std::ref(MRI2R[Ncycle*num_threads+i]),
								   std::ref(MRI2I[Ncycle*num_threads+i]),
								   std::ref(MRI3R[Ncycle*num_threads+i]),
								   std::ref(MRI3I[Ncycle*num_threads+i]));
		for(int i=0;i<(Nexceed);i++) threads_ds[i].join();
	
	#else
		
		for(int k=0;k<Nk;k++)
			RI(MKKRaux[k], MSPOaux[k], MRI1R[k], MRI1I[k], 
			   MRI2R[k], MRI2I[k], MRI3R[k], MRI3I[k]);
	
	#endif
	
		
	finish_t = clock();
	total_t = (double)(finish_t - start_t)/CLOCKS_PER_SEC;
	cout <<"Auxiliary Matrice at "<< pow(hbar*kvector[0],2.0)/2.0/me/ee <<" eV converged, ";
	cout <<"took "<< total_t << " s.\n" << endl;
	cout << endl;


	KKRcoef G1;
	vector<double> Xij(3), Rij(3);

	MatrixXcd tmpM0(Nsc,Nsc);
	MatrixXcd tmpM1(Nsc,Nsc);
	MatrixXcd tmpM2(Nsc,Nsc);
	MatrixXcd tmpM3(Nsc,Nsc);
	MatrixXcd tmpGstuct(Nsc,Nsc_I);
	MatrixXcd tmpGtrans(Nsc,Nsc_I);
	MatrixXcd tmpMtotal(Nsc_I,Nsc_I);
	
	
	for(int k=0;k<Nk;k++)
	{
		for(int m=0;m<HostsXcor[0].size();m++)
		{
			for(int n=0;n<ImpuritiesXcor[0].size();n++)
			{
				for(int nr=0;nr<3;nr++)
					Xij[nr] = HostsXcor[nr][m] - ImpuritiesXcor[nr][n];
		
				Cart2Sph(Xij, Rij);
				for(int lm1=0;lm1<Nlm;lm1++)
				{
					for(int lm2=0;lm2<Nlm_I;lm2++)
					{
						G1 = KKRcoef(Lv[lm1],Mv[lm1],Lv_I[lm2],Mv_I[lm2],Lv_I,Mv_I);
						tmpGstuct(m*Nlm+lm1,n*Nlm_I+lm2) = 
												G1.StuctConst(kvector[k],Rij);
						tmpGtrans(m*Nlm+lm1,n*Nlm_I+lm2) = 
												G1.TranslationConst(kvector[k],Rij);
	//					cout << tmpGtrans(m*Nlm+lm1, n*Nlm+lm2) << endl;
					}
				}
			}
		}
		
		for(int i=0;i<Nsc;i++)
		{
			for(int j=0;j<Nsc;j++)
			{
				tmpM0(i,j) = complex<double>(HostsScatterR[k][i][j],HostsScatterI[k][i][j]);
				tmpM1(i,j) = complex<double>(MRI1R[k][i][j],MRI1I[k][i][j]);
				tmpM2(i,j) = complex<double>(MRI2R[k][i][j],MRI2I[k][i][j]);
				tmpM3(i,j) = complex<double>(MRI3R[k][i][j],MRI3I[k][i][j]);
			}
		}
		tmpMtotal =	tmpGstuct.transpose()*tmpM0*tmpGstuct +
					tmpGstuct.transpose()*tmpM1*tmpGtrans +
					tmpGtrans.transpose()*tmpM2*tmpGstuct +
					tmpGtrans.transpose()*tmpM3*tmpGtrans;

		for(int i=0;i<Nsc_I;i++)
		{
			for(int j=0;j<Nsc_I;j++)
			{
				ImpuritiesRIreal[k][i][j] = tmpMtotal(i,j).real();
				ImpuritiesRIimag[k][i][j] = tmpMtotal(i,j).imag();
			}
		}
	}
}

