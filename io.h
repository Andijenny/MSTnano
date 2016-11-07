#ifndef IOSCATTERM_H_
#define IOSCATTERM_H_

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <Dense>
#include "GreensFunc.h"
#include "constant.h"

inline void ScatterMatrixPrint(
				const char* FileSR, const char * FileSI,
				std::vector<std::vector<std::vector<double> > > scatterR, 
				std::vector<std::vector<std::vector<double> > > scatterI)
{
	int rs[3];
	rs[0] = scatterR.size();
	rs[1] = scatterR[0].size();
	rs[2] = scatterR[0][0].size();
//	int is[3] = {scatterI.size(), scatterI[0].size(), scatterI[0][0].size()};
	int prs = rs[0]*rs[1]*rs[2];
	double *tmpsR = new double[prs];
	double *tmpsI = new double[prs];
//	double tmpsI[prs];
//	double tmpsI[prs];
	int cc = 0;
	for(int i=0;i<rs[0];i++)
	{
		for(int j=0;j<rs[1];j++)
		{
			for(int k=0;k<rs[2];k++)
			{
				tmpsR[cc] = scatterR[i][j][k];
				tmpsI[cc] = scatterI[i][j][k];
				cc++;
			}
		}
	}
//	for(int i=0;i<3;i++)
	std::ofstream foutR(FileSR,std::ios::out|std::ios::binary);
	std::ofstream foutI(FileSI,std::ios::out|std::ios::binary);
	
	foutR.write((char*)(rs),sizeof(int)*3);
	foutR.write((char*)tmpsR,sizeof(double)*prs);
	
	foutI.write((char*)(rs),sizeof(int)*3);
	foutI.write((char*)tmpsI,sizeof(double)*prs);
	
	foutR.close();
	foutI.close();
	delete [] tmpsR;
	delete [] tmpsI;
}

inline void ScatterMatrixRead(
					char* FileSR, char * FileSI, int nKstart, 
					std::vector<std::vector<std::vector<double> > >& scatterR, 
					std::vector<std::vector<std::vector<double> > >& scatterI)
{
	/*FileSR, FileSI: binary file which we want to read from;
	 * nKstart      : the value of the first dimension of stored matrix
	 *				  we determine to start to read when recalculate,
	 *				  e.g., assuming we have finished calculating and storing 
	 *				  4 energy points in the first run, then if we want to 
	 *				  read from the 2ed to the 4th, just set nKstart as 2.
	 *				  */
	
	std::ifstream finR(FileSR,std::ios::in|std::ios::binary);
	std::ifstream finI(FileSI,std::ios::in|std::ios::binary);
	if(!finR)
		throw("Failed to open file ScatterMR.");
	if(!finI)
		throw("Failed to open file ScatterMI.");
	
	int rs[3], is_i[3], rs_i[3];
	rs[0] = scatterR.size();
	rs[1] = scatterR[0].size();
	rs[2] = scatterR[0][0].size();
	
	finR.read((char*)rs_i,sizeof(int)*3);
	finI.read((char*)is_i,sizeof(int)*3);
	if((rs_i[0] < rs[0])|| 
	   (rs_i[1] < rs[1])||
	   (rs_i[2] < rs[2])  ) throw("wrong matrix size when read from ScatterMR.");


	int prs = rs[0]*rs[1]*rs[2];
	double tmpsR[prs];
	double tmpsI[prs];
	
	finR.read((char*)tmpsR, sizeof(double)*prs);
	finI.read((char*)tmpsI, sizeof(double)*prs);

	int Ncyc = nKstart;
	while(Ncyc > 1)
	{
		finR.read((char*)tmpsR, sizeof(double)*prs);
		finI.read((char*)tmpsI, sizeof(double)*prs);
		--Ncyc;
	}
	
	finR.close();
	finI.close();
	
	
	int cc = 0;
	for(int i=0;i<rs[0];i++)
	{
		for(int j=0;j<rs[1];j++)
		{
			for(int k=0;k<rs[2];k++)
			{
				scatterR[i][j][k] = tmpsR[cc];
				scatterI[i][j][k] = tmpsI[cc];
				cc++;
			}
		}
	}

}


inline void const RHOPrint(int Nk, int Np, 
						   double* rtip1, double* rtip2,
						   double* Etmp, GreenFunc GFobj)
{
//	static const double ee = 1.6e-19;
//	static const double pi = 4.0*atan(1.0);
//	const extern double ee;
//	const extern double pi;
	
	VectorXcd GF(Nk);
	VectorXd rho(Nk);
    
	vector<double> rr1(3);
	int Npos = Np*Np*3;
	int Nrho = Np*Np*Nk;
	double Spos[Npos];
	double Srho[Nrho];

	char FileRHO[100];
	sprintf(FileRHO,"RHO%.3f",Etmp[0]);
	ofstream rho_o(FileRHO,ios::out|ios::binary);

	int cc1 = 0,cc2 = 0;
	for(int j1=0;j1<Np;j1++)
	{
		rr1[2] = rtip1[2];
		if(Np > 1)
		{ rr1[0] = rtip1[0]+(j1)*(1.0)/(Np-1)*(rtip2[0]-rtip1[0]);}
		else
		{ rr1[0] = rtip1[0];}
		for(int j2=0;j2<Np;j2++)
		{	
			rr1[1] = rtip1[1]+(j2+1)*(1.0)/Np*(rtip2[1]-rtip1[1]);
			for(int l=0;l<3;l++)
				Spos[cc1++] = rr1[l];
			
			GFobj(rr1, rr1, GF);
	
			for(int i=0;i<Nk;i++)
			{
				rho(i) = -1.0/pi*GF(i).imag()*ee;//Unit states/eV
				Srho[cc2++] = rho(i);
			}

		}
	}
				

	rho_o.write((char*)&Npos,sizeof(int));
	rho_o.write((char*)&Nk,sizeof(int));
	rho_o.write((char*)&Nrho,sizeof(int));
	rho_o.write((char*)Spos,sizeof(double)*Npos);
	rho_o.write((char*)Etmp,sizeof(double)*Nk);
	rho_o.write((char*)Srho,sizeof(double)*Nrho);

	rho_o.close();
}

//inline void const LDOSPrint(int Nk, int Nsp, double * Etmp,
//							vector<double> IntegRadius, GreenFunc GFobj)
//{
//	vector<vector<double> > tmpLDOS(Nsp, vector<double> (Nk,0.0));
//	double Sld[Nsp*Nk];
//	int c1 = 0;
//	for(int i=0;i<Nsp;i++)
//	{
//		GFobj.ASALDOS(IntegRadius,i,tmpLDOS[i]);
//		for(int k=0;k<Nk;k++)
//		{
//			Sld[c1++] = tmpLDOS[i][k];
//		}
//	}
//	
//	ofstream ldos_o("LDOS",ios::out|ios::binary);
//	ldos_o.write((char*)&Nsp,sizeof(int));
//	ldos_o.write((char*)&Nk,sizeof(int));
//	ldos_o.write((char*)Etmp,sizeof(double)*Nk);
//	ldos_o.write((char*)Sld,sizeof(double)*Nsp*Nk);
//	ldos_o.close();		
//
//}

#endif
