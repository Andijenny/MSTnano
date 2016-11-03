#ifdef MULTI_THREAD
	#include<thread>
#endif
#include "ParInitial.h"
#include "io.h"
#include "ParallelTmatrix.h"

//#include "Tmatrix2.h"
//#include <algorithm>

#include <iomanip>
//#include <fstream>
#include <Dense>
//#include<mutex>
//#include "Geometry.h"
#include "GreensFunc.h"
#include "ScatteringPathOperator.h"
//#include "ReNormInteractor.h"
//#include "cdefs.h"

//std::mutex some_mutex;

extern "C"{
	extern char _binary_VersionInfo_start;
	extern char _binary_VersionInfo_end;
}


int main()
{	
//	__COPYRIGHT("@(#)Copyright(c) 2016, 2166\n\
//			The Electronic Department of the Peking University.\n\
//			All rights reserved.\n");
	
	 char*  p = &_binary_VersionInfo_start;
	 while ( p != &_binary_VersionInfo_end ) cout << *p++;
	
	clock_t start_t;
	start_t = clock();
    cout << "*******************Main Program*****************\n";
	cout << endl;
	cout << endl;
	cout << "-------------------Starting in -----------------\n";
    cout << endl;
	time_t now_time;
	char tmp[64];
	now_time = time(0);
	strftime(tmp,sizeof(tmp),"\t\t%Y-%m-%d %A \n\t\t %H:%M:%S ",localtime(&now_time));
	puts(tmp);
//	cout << now_time;
//	CTime tt = Ctime::GetCurrentTime();
//	system("time");
//	cout << t.GetMonth() << "." << t.GetDay() << ", " << t.GetYear();
//	cout << " " << GetHour() <<":"<< GetMinute()<<":"<<GetSecond()<<endl; 
    cout << endl;
    cout << endl;
	ParInitial pinit;
	char *parameterFile = "par.in";
	bool Lpar=pinit.ParIn(parameterFile);
	if(Lpar) cout << "ParIn found parameters from par.in.\n";
	cout << endl;
	cout << endl;
	pinit.PrintParIn();	
	
	int NN = pinit.NPotentialPoints;
	
	vector<vector<double> >Xcor(3,vector<double>(10000,0.0));  //coordinates
	vector<string> SA(10000);  //ChemicalSpecies
	vector<int> SL(10000, 0);     //SpeciesLabel
	vector<vector<double> >rfit(10000,vector<double>(NN,0.0)); //radius r
	vector<vector<double> >pfit(10000,vector<double>(NN,0.0)); //potential(r)
	vector<double> IntegRadius(10000,0.0);
	char *geometryFile = "g.in";
	pinit.GeomIn(geometryFile,Xcor, SA, SL, rfit, pfit, IntegRadius);

	vector<vector<double> >HostsXcor(3,vector<double>(Xcor[0].size(),0.0));
	vector<vector<double> >ImpuritiesXcor(3,vector<double>(Xcor[0].size(),0.0));
	vector<int> HostsSL(SL.size());     
	vector<int> ImpuritiesSL(SL.size());     

	int ImpL[2];
	for(int i=0;i<2;i++)
		ImpL[i] = pinit.ImpuritiesLabel[i];
	int cH = 0, cI = 0;
	for(int i=0;i<Xcor[0].size();i++)
	{
		if(((i+1) > ImpL[1]) || ((i+1) < ImpL[0]))
		{
			for(int j=0;j<3;j++)
				HostsXcor[j][cH] = Xcor[j][i];
			HostsSL[cH] = SL[i];     
			cH++;
		}
		else
		{
			for(int j=0;j<3;j++)
				ImpuritiesXcor[j][cI] = Xcor[j][i];
	        ImpuritiesSL[cI] = SL[i];     
			cI++;
		}
	}
	for(int i=0;i<3;i++)
	{
		HostsXcor[i].resize(cH);
		ImpuritiesXcor[i].resize(cI);
	}
	HostsSL.resize(cH);
	ImpuritiesSL.resize(cI);     


	cout << "Finish reading geometric information from g.in. \n";
	cout << endl;
	cout << endl;
	cout << "Have loaded "<< Xcor[0].size() << " atoms ";
	cout << "that pertain to "<<SA.size() <<" species,\n" << endl;
	cout << "including "<< HostsXcor[0].size() << " host atoms ";
	cout << "and "<< ImpuritiesXcor[0].size() << " impurities atoms.\n ";
	cout << endl;
	cout << endl;
	cout << "Finish Interpolating "<< NN <<" points to potential(r).\n";
	cout << endl;
	cout << endl;

	
	
	int Nk = pinit.Nkpoints;

	vector<int> parallelPar(3,1);
	
	#ifdef MULTI_THREAD

//		int num_threads=min((int)(std::thread::hardware_concurrency()),Nk);
		int num_threads=min((int)(std::thread::hardware_concurrency())-1,Nk);
	
		if(num_threads > 1)
			cout << "Multithread traggered, running on " << num_threads <<" threads.\n";
//		cout << "num_threads " << num_threads << endl;
		cout << endl;
		cout << endl;
		
		if(Nk > num_threads)
			Nk = Nk-Nk%num_threads;// concurrent more effeciently
//		int Ncycle = Nk/num_threads;
//		int Nexceed = Nk-Ncycle*num_threads;
		parallelPar[0] = num_threads;
		parallelPar[1] = Nk/num_threads;
		parallelPar[2] = Nk-parallelPar[1]*num_threads;

	#endif

	
	double Ema = pinit.Emax;
	double Emi = pinit.Emin; // in unit of eV
	VecDoub kvector(Nk);
	double Etmp[Nk];
	vector<double> CopyKvector(Nk);
	for(int i=0;i<Nk;i++)
	{
		if(Nk > 1)
		{
			Etmp[i] = Emi+i*(Ema-Emi)/(Nk-1);
			kvector[i] = sqrt(2*me*ee*Etmp[i])/hbar;
		}
		else
		{
			Etmp[i] = Emi;
			kvector[i] = sqrt(2*me*ee*Emi)/hbar;
		}
		CopyKvector[i] = kvector[i];
	}
	cout << "Initialized the range of energy E as "<< Emi<<" to\n";
	cout << endl;
	cout << Ema << " eV, altogther "<<Nk<<" points.\n";
	cout << endl;
	cout << endl;
	int Lmax=pinit.Lmax;
	int Lmax_I=pinit.ImpuritiesLmax;
//	int Nlm=(Lmax+1);
	
	cout << "The Angular quantum numbers l's truncate to "<<Lmax<<" for host atoms,\n";
	cout << "while trancate to "<<Lmax_I<<" for impurities.\n";
	cout << endl;
	cout << endl;

	vector<vector<vector<double> > > Tll(Nk,vector<vector<double> >
			(SA.size(),vector<double>(max(Lmax, Lmax_I)+1, 0.0)));

	if(!pinit.LReadTMatrix)
	{
		cout << "Begin self-consistent calculation for t-matrix\n";	
		cout << endl;
		cout << endl;
		double mw = pinit.WeightMixing;
		double cp = pinit.TMConvergePrecision;
//	cout << Tll.size() <<" * "<<Tll[0].size()<<Tll[0][0].size()<<endl;
		cout << "-----------------------------" <<endl;
		cout << "|                           |" <<endl;
		for(int i=0;i<SA.size();i++)
		{
			cout << "|       load " << setw(4)<<i<<": "<<SA[i]<<"    |" <<endl;
//		cout << endl;
		}
		cout << "|                           |" <<endl;
		cout << "-----------------------------" <<endl;
		cout << endl;
		cout << endl;
	
		ParallelTmatrix(parallelPar, rfit, pfit, SA, Lmax, mw, cp, kvector, Tll);
	
		cout << "End t-matrix calculation\n";	
		cout << endl;
		cout << endl;
		TmatrixPrint(Tll);  
	}
	else
	{
		TmatrixRead(Tll);
		cout << "Have read t-matrix from Tll_matrix\n";
		cout << endl;
		cout << endl;
	}
	
	
	double SPOcp = pinit.SPOConvergePrecision;
	vector<vector<double> > RBasis(3, vector<double>(3));
	vector<int> Kmesh(3);
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			RBasis[i][j]= pinit.RealSpaceBasis[i][j];
	for(int i=0;i<3;i++)
		    Kmesh[i] = pinit.Kmesh[i];

	
	int Nlm=(Lmax+1)*(Lmax+1);
	vector<int> Lv(Nlm);
	vector<int> Mv(Nlm);
	for(int l = 0, c1 = -1;l<=Lmax;l++){
		for(int m=-l;m<=l;m++){
			c1++;
			Lv[c1] = l;
			Mv[c1] = m;
//			Mv[c1] = 0;
		}
	}

	
	
	int Nlm_I=(Lmax_I+1)*(Lmax_I+1);
	vector<int> Lv_I(Nlm_I);
	vector<int> Mv_I(Nlm_I);
	for(int l = 0, c1 = -1;l <= Lmax_I; l++)
		for(int m = -l; m <= l; m++)
		{
			c1++;
			Lv_I[c1] = l;
			Mv_I[c1] = m;
		}
	

	int Nspecies=HostsXcor[0].size();
	int Nsc=Nspecies*Nlm;


    ScatteringPathOperator Tao(SPOcp, Lv, Mv, parallelPar, CopyKvector,
								Tll, RBasis, Kmesh); 	
	
	vector<vector<vector<vector<double> > > > HoststR(Nspecies,
					 vector<vector<vector<double> > >(Nk,
							  vector<vector<double> >(Nlm,
									   vector<double>(Nlm,0.0))));
    
	vector<vector<vector<double> > > HostsScatterR(Nk,
						   vector<vector<double> >(Nsc,
									vector<double>(Nsc,0.0)));
	vector<vector<vector<double> > > HostsScatterI(Nk,
						   vector<vector<double> >(Nsc,
									vector<double>(Nsc,0.0)));

	int NKinBZ = Kmesh[0]*Kmesh[1]*Kmesh[2];

	vector<vector<vector<vector<complex<double> > > > > MKKRaux(Nk,
		   vector<vector<vector<complex<double> > > >(NKinBZ,
				  vector<vector<complex<double> > >(Nsc,
					       vector<complex<double> >(Nsc, complex<double>(0.0, 0.0)))));
	
	vector<vector<vector<vector<complex<double> > > > > MSPOaux(Nk,
		   vector<vector<vector<complex<double> > > >(NKinBZ,
				  vector<vector<complex<double> > >(Nsc,
					       vector<complex<double> >(Nsc, complex<double>(0.0, 0.0)))));
    
	
	std::thread::id master_thread = std::this_thread::get_id();
	
	
	
	char FileHostsSR[100];
	char FileHostsSI[100];
	sprintf(FileHostsSR,"HostsScatterMR%.3f",Etmp[0]);
	sprintf(FileHostsSI,"HostsScatterMI%.3f",Etmp[0]);
	
	
	if(!pinit.LReadHostsSMatrix)
	{
		cout << "Start iteration of the host scattering path operator.\n";
		cout << endl;
		cout << endl;
		Tao(HostsSL, HostsXcor, HoststR, HostsScatterR, HostsScatterI, MKKRaux, MSPOaux);
		
		cout << "The expansion of the host scattering path operator converged.\n";
		cout << endl;
		cout << endl;
		
		if(pinit.LWriteHostsSMatrix)
			ScatterMatrixPrint(FileHostsSR, FileHostsSI, HostsScatterR, HostsScatterI);
	}
	else
	{
		ScatterMatrixRead(FileHostsSR, FileHostsSI, pinit.nKstart , HostsScatterR, HostsScatterI);
		Tao.ConstructTmatrix(HostsSL, HoststR);
		cout << "Have read Scattering Path Matrix for Host atoms.\n"<< endl << endl;
	}
	
	int Nspecies_I = ImpuritiesXcor[0].size();
	int Nsc_I = Nspecies_I*Nlm_I;

	vector<vector<vector<vector<double> > > > ImpuritiestR(Nspecies_I,
						  vector<vector<vector<double> > >(Nk,
								   vector<vector<double> >(Nlm_I,
											vector<double>(Nlm_I, 0.0))));

	vector<vector<vector<double> > > ImpuritiesScatterR(Nk,
								vector<vector<double> >(Nsc_I,
										 vector<double>(Nsc_I,0.0)));
	vector<vector<vector<double> > > ImpuritiesScatterI(Nk,
								vector<vector<double> >(Nsc_I,
										 vector<double>(Nsc_I,0.0)));



	vector<vector<vector<double> > > ImpuritiesRIreal(Nk,
							  vector<vector<double> >(Nsc_I,
									   vector<double>(Nsc_I,0.0)));
	vector<vector<vector<double> > > ImpuritiesRIimag(Nk,
							  vector<vector<double> >(Nsc_I,
									   vector<double>(Nsc_I,0.0)));
	char FileImpuritiesSR[100];
	char FileImpuritiesSI[100];
	sprintf(FileImpuritiesSR,"ImpuritiesScatterMR%.3f",Etmp[0]);
	sprintf(FileImpuritiesSI,"ImpuritiesScatterMI%.3f",Etmp[0]);
	
	if (Nspecies_I > 0)
	{
		if(!pinit.LReadHostsSMatrix)
		{
	
			cout << "Iteration of the impurities' scattering path operator.\n";
			cout << endl;
			cout << endl;
			Tao.ImpuritiesSPO(Lv_I, Mv_I, ImpuritiesSL, ImpuritiesXcor, ImpuritiestR, 
								  ImpuritiesScatterR, ImpuritiesScatterI);
		//	Tao.ImpuritiesSPO(ImpuritiesSL, ImpuritiesXcor, ImpuritiestR, ImpuritiesScatterR, ImpuritiesScatterI);
			
			cout << "Determination of renormonization interactor for impurities."<< endl;
			cout << endl;
			cout << endl;
			Tao.ReNormInteractor(Lv_I, Mv_I, HostsXcor, ImpuritiesXcor, HostsScatterR, HostsScatterI,
								 MKKRaux, MSPOaux, ImpuritiesRIreal, ImpuritiesRIimag);
	//		Tao.ReNormInteractor(HostsXcor, ImpuritiesXcor, HostsScatterR, HostsScatterI,
	//							 MKKRaux, MSPOaux, ImpuritiesRIreal, ImpuritiesRIimag);
			
			
		
	//		cout << "Impurities RI stores over k points."<<endl;	
	//		cout << endl;
	//		cout << endl;
			MatrixXcd scattertmp_I(Nsc_I, Nsc_I);
			MatrixXcd RItmp_I(Nsc_I, Nsc_I);
			MatrixXcd Mtmp_I(Nsc_I, Nsc_I);
	
			for(int k=0;k<CopyKvector.size();k++)
			{
				for(int i=0;i<Nsc_I;i++)
				{
					for(int j=0;j<Nsc_I;j++)
					{
						Mtmp_I(i,j) = complex<double>(ImpuritiesScatterR[k][i][j], ImpuritiesScatterI[k][i][j]);
						RItmp_I(i,j) = complex<double>(ImpuritiesRIreal[k][i][j], ImpuritiesRIimag[k][i][j]);
					}
				}
				scattertmp_I = Mtmp_I.inverse()-RItmp_I;
				scattertmp_I = scattertmp_I.inverse();
				for(int i=0;i<Nsc_I;i++)
				{
					for(int j=0;j<Nsc_I;j++)
					{
						ImpuritiesScatterR[k][i][j] = scattertmp_I(i,j).real();
						ImpuritiesScatterI[k][i][j] = scattertmp_I(i,j).imag();
					}
				}
			}
			
			cout << "The expansion of the impurities' scattering path operator converged.\n"; 
			cout << endl;
			cout << endl;
			if(pinit.LWriteImpuritiesSMatrix)
				ScatterMatrixPrint(FileImpuritiesSR, FileImpuritiesSI, 
								   ImpuritiesScatterR, ImpuritiesScatterI);
		}
		else
		{
			ScatterMatrixRead(FileImpuritiesSR, FileImpuritiesSR, pinit.nKstart,
							  ImpuritiesScatterR, ImpuritiesScatterI);
			Tao.ConstructTmatrix(ImpuritiesSL, ImpuritiestR);
			
			cout << "Have read Scattering Path Matrix for Impurities.\n"<<endl<<endl;
		}
	}

//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------

	cout << "Generate Green's function in real space." << endl;
	cout << endl;
	cout << endl;
	
	GreenFunc GFobj;
	if(Nspecies_I > 0)
	{
		GFobj = GreenFunc(Lmax_I, CopyKvector, Lv_I, Mv_I, ImpuritiesXcor,
						  ImpuritiestR, ImpuritiesScatterR, ImpuritiesScatterI);

//		GFobj = GreenFunc(Lmax, CopyKvector, Lv, Mv, Xcor_t,
//						  tR_t, ImpuritiesScatterR, ImpuritiesScatterI);
		cout << "export Green's function around impurities."<<endl;
		cout << endl;
		cout << endl;
	}
	else
	{
		GFobj = GreenFunc(Lmax, CopyKvector, Lv, Mv, HostsXcor,
						  HoststR, HostsScatterR, HostsScatterI);
		cout << "export Green's function over undisturbed host atoms."<<endl;
		cout << endl;
		cout << endl;
	}
	
	cout <<"Write RHO ...\n"<< endl;
	cout << endl;
    
	double rtip1[3];
    double rtip2[3];
	for(int i=0;i<3;i++)
	{
		rtip1[i] = pinit.TipPosition1[i];
		rtip2[i] = pinit.TipPosition2[i];
	}
	int Np = pinit.NTipPoints;

	RHOPrint(Nk, Np, rtip1, rtip2, Etmp, GFobj);	
	
	

	clock_t finish_t=clock();
    double total_t=(double)(finish_t-start_t)/CLOCKS_PER_SEC;
	cout <<(int)total_t/3600 <<" hours "<<((int)total_t-(int)(total_t/3600)*3600)/60;
	cout <<" minutes "<<(int)total_t%60 <<" s elapsed.\n"<<endl;
	cout << endl;
	cout << "End  Program in "<<endl;
	cout << endl;
	now_time = time(0);
	strftime(tmp,sizeof(tmp),"\t\t\t%Y-%m-%d %A \n\t\t\t %H:%M:%S ",localtime(&now_time));
	puts(tmp);
	cout << endl;
	cout << endl;
	cout << "End of file. "<<endl;
	return 0;
}
