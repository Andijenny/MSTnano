#include "Geometry.h"
#include <fstream>
#include <sstream>
#include <string>
#include <Dense>
using namespace Eigen;
void Geom::Rcoor(vector<vector<double> > & rr, vector<string> & SA,
				 vector<int> & SL, vector<double> &  IntegRadius){
	ifstream fin;
	fin.open(fn);
	string buf;
//	char *p;
//	vector<vector<double> > cor(3,vector<double>(1000,0));
//	char species[1000];
	vector<string> ss(100000);
	int i=0,fl=0;
	int cc=-1;
	while(getline(fin,buf))
	{
		if(buf == "")
			continue;
		istringstream instr(buf);
		instr >> rr[0][i] >> rr[1][i] >> rr[2][i] >> ss[i];
		if(i==0){
			SL[0]=0;
			cc++;
			SA[cc]=ss[i];
		}
		for(int j=0;j<i;j++){
			if(ss[i]!=ss[j])
			{ 
				fl = 1; 
				SL[i] = (SL[i]>SL[j]?SL[i]:SL[j]);
			}
			else
			{
				fl = 0;
				SL[i] = SL[j];
				break;
			}
		}
		if(fl==1){
			SL[i]++;
			cc++;
			SA[cc]=ss[i];
		}
		i++;
	}
	
	ss.resize(i);
	SL.resize(i);
	SA.resize(cc+1);
	IntegRadius.resize(i);
	
	for(int m=0;m<3;m++) rr[m].resize(i);
	

	for(int k=0;k<rr.size();k++)
		for(int j=0;j<rr[0].size();j++)
			rr[k][j]=rr[k][j]*1.0e-10;
	
	for(int k=0;k<rr[0].size();k++)
		IntegRadius[k] = MinDis2(k, rr)/2.0;
	
}

double Geom::MinDis2(int index, vector<vector<double> > Xcor)
{
	int Ncols = Xcor[0].size();

	VectorXd dis(Ncols);
	
	for(int i=0; i< Ncols; i++)
	{
		if(i != index)
		{
			dis(i) = sqrt(pow(Xcor[0][index]-Xcor[0][i],2)+pow(Xcor[1][index]-Xcor[1][i],2)
						  +pow(Xcor[2][index]-Xcor[2][i],2));
		}
		else
		{
			dis(i) = 100.0;
		}
	}
	int in_d;
	double ms = dis.minCoeff(& in_d);
	return ms;
}


void Cart2Sph(vector<double> & x, vector<double> & r)
{
	r[0]=norm_hw(x);
	if(abs(r[0]) < 1.0e-20) throw("The input radius is zero in Cart2Sph().");
	r[1]=acos(x[2]/r[0]);
//	for(int i=0;i<3;i++)
//		tmp[i]=x[i]*sin(r[1]);
//	r[2]=acos(tmp[0]/norm_hw(tmp));
	double dis = sqrt(x[0]*x[0]+x[1]*x[1]);
	if(abs(dis) > 1.0e-20){ r[2] = acos(x[0]/dis);}
	else				  { r[2] = 0.0;}
}

void Cart2Sph(VectorXd & x, VectorXd & r)
{
	r(0) = sqrt(x.squaredNorm());
	r(1) = acos(x(2)/r(0));
	r(2) = acos(x(0)/r(0)/sin(r(1)));
}
