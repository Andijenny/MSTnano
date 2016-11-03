#ifndef PARINIT_H_GG
#define PARINIT_H_
//#include "nr3.h"
#include "Geometry.h"
#include<iostream>
#include <iomanip>
#include<string>
#include<string.h>
#include<fstream>
#include<sstream>
#include<vector>
#include "potential.h"
#include "interp_1d.h"


using namespace std;

const Doub me=9.10938188e-31;
const Doub hbar=1.055571596e-34;
const Doub ee=1.602176462e-19;
const Doub Br=5.291772083e-11;
const Doub pi=4.0*atan(1.0);

struct ParInitial
{
	bool LReadTMatrix;
	bool LReadHostsSMatrix;
	bool LWriteHostsSMatrix;
	bool LReadImpuritiesSMatrix;
	bool LWriteImpuritiesSMatrix;
	int nKstart;
	double WeightMixing;
	double TMConvergePrecision;
	double SPOConvergePrecision;
	int	   Nkpoints;
	double Emin;
	double Emax;
	int Lmax;
	int ImpuritiesLmax;
	int NPotentialPoints;
	int NTipPoints;
	double TipPosition1[3];
	double TipPosition2[3];
	int ImpuritiesLabel[2];
	double RealSpaceBasis[3][3];
	int Kmesh[3];
	ParInitial(){};
	inline void PrintParIn();
	inline bool ParIn(char* filename);
	inline void GeomIn(char* fn, vector<vector<double> >& Xcor,
				vector<string>& ChemicalSpecies,
				vector<int>& SpeciesLabel,
				vector<vector<double> >& rfit,
				vector<vector<double> >& pfit,
				vector<double>& IntegRadius);
	inline string printLogical(bool x);
};

bool ParInitial::ParIn(char* filename)
{
	ifstream fin;
	fin.open(filename);
	string buf;
	string sTm="LReadTMatrix";
	string sSrm="LReadHostsSMatrix";
	string sSwm="LWriteHostsSMatrix";
	string sSrmi="LReadImpuritiesSMatrix";
	string sSwmi="LWriteImpuritiesSMatrix";
	string sNst="nKstart";
	string sWm="WeightMixing";
	string sCp1="TMConvergePrecision";
	string sCp2="SPOConvergePrecision";
	string sNk="Nkpoints";
	string sNtip="NTipPoints";
	string sLm="Lmax";
	string sLmi="ImpuritiesLmax";
	string sEmi="Emin";
	string sEma="Emax";
	string sNPp="NPotentialPoints";
	string sTp1="TipPosition1";
	string sTp2="TipPosition2";
	string sRsb="RealSpaceBasis";
	string sKm="Kmesh";
	string sImp="ImpuritiesLabel";

	WeightMixing = 0.2;
	TMConvergePrecision = 1.0e-5;
	SPOConvergePrecision = 1.0e-15;
	Nkpoints = 100;
	Lmax = 1;
	ImpuritiesLmax = 2;
	NPotentialPoints = 10000;
	NTipPoints = 1000;
	Emin = 0.0;
	Emax = 1.0;
	LReadTMatrix = false;
	LReadHostsSMatrix = false;
	LWriteHostsSMatrix = false;
	LReadImpuritiesSMatrix = false;
	LWriteImpuritiesSMatrix = false;
	nKstart = 1;
	for(int i=0;i<3;i++)
	{
		TipPosition1[i] = 0.0;
		TipPosition2[i] = 1.0e-10;
		Kmesh[i] = 2;
	}
	Kmesh[2] = 1;
	ImpuritiesLabel[0] = 10000;
	ImpuritiesLabel[1] = 20000;
	
	bool fflag = false;
	string str;
	double ptemp[3];
	while(getline(fin,buf))
	{
		if(buf == "") continue;
		istringstream instr(buf);
		instr >> str;
//		cout << "success"<< endl;
//		cout << str[i]<< par[i]<<ptemp[0]<<ptemp[1]<< endl;
		if (str == sTm)
		{
			instr >> str;
			if((str == "true") || (str == "True") || (str == "TRUE"))
				LReadTMatrix = true;
			else
				LReadTMatrix = false;
			fflag = true;
		}
		else if (str == sSrm)
		{
			instr >> str;
			if((str == "true") || (str == "True") || (str == "TRUE"))
				LReadHostsSMatrix = true;
			else
				LReadHostsSMatrix = false;
			fflag = true;
		}
		else if (str == sSrmi)
		{
			instr >> str;
			if((str == "true") || (str == "True") || (str == "TRUE"))
				LReadImpuritiesSMatrix = true;
			else
				LReadImpuritiesSMatrix = false;
			fflag = true;
		}
		else if (str == sSwm)
		{
			instr >> str;
			if((str == "true") || (str == "True") || (str == "TRUE"))
				LWriteHostsSMatrix = true;
			else
				LWriteHostsSMatrix = false;
			fflag = true;
		}
		else if (str == sSwmi)
		{
			instr >> str;
			if((str == "true") || (str == "True") || (str == "TRUE"))
				LWriteImpuritiesSMatrix = true;
			else
				LWriteImpuritiesSMatrix = false;
			fflag = true;
		}
		else if (str == sNst)
		{
			instr >> ptemp[0];
			nKstart = (int)ptemp[0];
			fflag = true;
		}
		else if (str == sWm)
		{
			instr >> ptemp[0];
			WeightMixing = ptemp[0];
			fflag = true;
		}
		else if(str == sCp1)
		{
			instr >> ptemp[0];
			TMConvergePrecision = ptemp[0];
			fflag = true;
		}
		else if(str == sCp2)
		{
			instr >> ptemp[0];
			SPOConvergePrecision = ptemp[0];
			fflag = true;
		}
		else if(str == sNtip)
		{
			instr >> ptemp[0];
			NTipPoints = (int)ptemp[0];
			fflag = true;
		}
		else if(str == sNk)
		{
			instr >> ptemp[0];
			Nkpoints = (int)ptemp[0];
			fflag = true;
		}
		else if(str == sLm)
		{
			instr >> ptemp[0];
			Lmax = (int)ptemp[0];
			fflag = true;
		}
		else if(str == sLmi)
		{
			instr >> ptemp[0];
			ImpuritiesLmax = (int)ptemp[0];
			fflag = true;
		}
		else if(str == sNPp)
		{
			instr >> ptemp[0];
			NPotentialPoints = (int)ptemp[0];
			fflag = true;
		}
		else if(str == sEmi)
		{
			instr >> ptemp[0];
			Emin = ptemp[0];
			fflag = true;
		}
		else if(str == sEma)
		{
			instr >> ptemp[0];
			Emax = ptemp[0];
			fflag = true;
		}
		else if(str == sTp1)
		{
			instr >> ptemp[0] >> ptemp[1] >> ptemp[2];
			for(int i=0;i<3;i++)
				TipPosition1[i] = ptemp[i]*1.0e-10;
			fflag = true;
		}
		else if(str == sTp2)
		{
			instr >> ptemp[0] >> ptemp[1] >> ptemp[2];
			for(int i=0;i<3;i++)
				TipPosition2[i] = ptemp[i]*1.0e-10;
			fflag = true;
		}
		else if(str == sImp)
		{
			instr >> ptemp[0] >> ptemp[1];
			ImpuritiesLabel[0]=(int)ptemp[0];
			ImpuritiesLabel[1]=(int)ptemp[1];
			fflag = true;
		}
		else if(str == sRsb)
		{	
			for(int i=0;i<3;i++)
			{
				getline(fin,buf);
				istringstream instr(buf);
				instr >> RealSpaceBasis[i][0] 
					  >> RealSpaceBasis[i][1]
					  >> RealSpaceBasis[i][2];
				for(int j=0;j<3;j++)
				{
					 RealSpaceBasis[i][j] *= 1.0e-10;
				
				}
			}
			fflag = true;
		}
		else if(str == sKm)
		{
			instr >> Kmesh[0] >> Kmesh[1] >> Kmesh[2];
			fflag = true;
		}
		else
		{
			instr >> str;
		}
	}
	return fflag;
}

void ParInitial::GeomIn(char* fn,
						vector<vector<double> >& Xcor,
						vector<string>& ChemicalSpecies,
						vector<int>& SpeciesLabel,
						vector<vector<double> >& rfit,
						vector<vector<double> >& pfit,
						vector<double>& IntegRadius)
{
//	char * fn="g.in";
	Geom gg(fn);
	gg.Rcoor(Xcor, ChemicalSpecies, SpeciesLabel, IntegRadius);
	int Nspecies = ChemicalSpecies.size();
//	cout << Nspecies << endl;
	int NN = NPotentialPoints;
	vector<vector<double> >rin(Nspecies,vector<double>(NN,0.0));
	vector<vector<double> >pin(Nspecies,vector<double>(NN,0.0));
	rfit.resize(Nspecies);
	pfit.resize(Nspecies);
//	cout << "rfit size: "<< rfit.size() << "*" << rfit[0].size() << endl;
	pot p1;
	char* potfile;
	string sufx=".pot";
	Poly_interp fp;
	VecDoub rtemp;
	VecDoub ptemp;
	vector<double> tmpR(3,0.0);
	for(int i=0;i<ChemicalSpecies.size();i++){
		ChemicalSpecies[i] = ChemicalSpecies[i]+sufx;
		potfile = (char*)ChemicalSpecies[i].c_str();
		p1 = pot(potfile);
		p1.rv(rin[i]);
		p1.pv(pin[i]);
		rtemp = VecDoub(rin[i].size(),0.0);
		ptemp = VecDoub(pin[i].size(),0.0);
//		for(int m=0;m<3;m++) tmpR[m] = Xcor[m][i];
//		cout <<"tmpR: "<<tmpR[0] << tmpR[1] << tmpR[2]<< endl;
		for(int j=0;j<rin[i].size();j++){
			rin[i][j]*=1.0e-10;
			pin[i][j]*=2*me/hbar/hbar*ee;
			rtemp[j]=rin[i][j];
			ptemp[j]=pin[i][j];
		}
		fp=Poly_interp(rtemp,ptemp,4);
    	for(int j=0;j<NN;j++){
			rfit[i][j]=rin[i][0]+(double)(j*1.0/(NN-1))*(rin[i][rin[i].size()-1]-rin[i][0]);
	//		rfit[i][j]=rin[i][0]+(double)(j*1.0/(NN-1))*(IntegRadius[i]-rin[i][0]);
			pfit[i][j]=fp.interp(rfit[i][j]);
//			cout <<" aa: "<<rfit[i][j]<<"		"<<pfit[i][j] << endl;
		}

	}
}

string ParInitial::printLogical(bool x)
{
	if(x)
	{
		return "true";
	}
	else
	{
		return "false";
	}
}


void ParInitial::PrintParIn()
{

	cout << setw(20)<<"NpotentialPoints:	"<< NPotentialPoints << endl;
	cout << setw(20)<<"RealSpaceBasis:\n";
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
			cout << setw(20)<< RealSpaceBasis[i][j];
		cout << endl;
	}
	cout << setw(20)<<"NKmeshPoints:	"<< Kmesh[0]<<"	";
	cout << Kmesh[1]<<"	" << Kmesh[2]<< endl;
	cout << setw(20)<<"Lmax:	"<< Lmax << endl;
	cout << setw(20)<<"WeightMixing:	"<< WeightMixing << endl;
	cout << setw(20)<<"TMConvPrecision:	"<< TMConvergePrecision << endl;
	cout << setw(20)<<"SPOConvPrecision:	"<< SPOConvergePrecision << endl;
	cout << setw(20)<<"LReadTMatrix:	" << printLogical(LReadTMatrix) << endl;
	cout << setw(20)<<"LReadHostsSMatrix:	" << printLogical(LReadHostsSMatrix) << endl;
	cout << setw(20)<<"nKstart:	" << nKstart << endl;
	cout << setw(20)<<"LWriteHostsSMatrix:	" << printLogical(LWriteHostsSMatrix) << endl;
	cout << setw(20)<<"LReadImpuritiesSMatrix:	" << printLogical(LReadImpuritiesSMatrix) << endl;
	cout << setw(20)<<"LWriteImpuritieSMatrix:	" << printLogical(LWriteImpuritiesSMatrix) << endl;
	cout << setw(20)<<"Emin(eV):	"<< Emin << endl;
	cout << setw(20)<<"Emax(eV):	"<< Emax << endl;
	cout << setw(20)<<"NEnergyPoints:	"<< Nkpoints << endl;
	cout << setw(20)<<"ImpuritiesLabel:	";
	cout << ImpuritiesLabel[0]<<"	" << ImpuritiesLabel[1]<< endl;
	cout << setw(20)<<"NTipPoints:	"<< NTipPoints << endl;
	cout << setw(20)<<"PositionTip1(m):	"<< TipPosition1[0]<<"	";
	cout << TipPosition1[1]<<"	" << TipPosition1[2]<< endl;
	cout << setw(20)<<"PositionTip2(m):	"<< TipPosition2[0]<<"	";
	cout << TipPosition2[1]<<"	" << TipPosition2[2]<< endl;
}

#endif
