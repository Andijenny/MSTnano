//#include <cmath>
#include "GreensFunc.h"
#include "Geometry.h"
#include "besselfrac.h"
#include "plegendre.h"
#include "quadd.h"

const double pi = 4.0*atan(1.0);

int MinDis(vector<double> rr, vector<vector<double> > Xcor)
{
	int Ncols = Xcor[0].size();
	VectorXd dis(Ncols);
	for(int i=0; i< Ncols; i++)
	{
		dis(i) = sqrt(pow(rr[0]-Xcor[0][i],2)+pow(rr[1]-Xcor[1][i],2)
				 +pow(rr[2]-Xcor[2][i],2));
	}
	int index;
//	double ms = dis.minCoeff(& index);
	dis.minCoeff(& index);
	return index;
}

void ZFunc::JL(int index, vector<double> rr, double kv, VectorXcd & J0, VectorXcd & H0)
{
	int Nlm=Lv.size();
	
//	int index = MinDis(rr, Xcor);
	vector<double>Rrel(3,0.0);
	for(int i=0;i<3;i++)
//		Rrel[i] = rr[i];
		Rrel[i]=rr[i]-Xcor[i][index];

	vector<double>Rij(3,0.0);
  	Cart2Sph(Rrel,Rij);	
	
	Bessel besj;
	Bessel besy;
	for(int lm=0; lm<Nlm; lm++)
	{
		J0(lm)=besj.sphbesj(Lv[lm],kv*Rij[0])*SpHamnx(Lv[lm],Mv[lm],Rij[1],Rij[2]);
		H0(lm)=complex<double>(besj.sphbesj(Lv[lm],kv*Rij[0]),besy.sphbesy(Lv[lm],kv*Rij[0]))
				*SpHamnx(Lv[lm],Mv[lm],Rij[1],Rij[2]);
	}
}


//vector<double> ZFunc:: operator()(vector<double> rr)
//void ZFunc:: operator()(vector<double> rr, vector<vector<complex> > & ZL)
void ZFunc::operator()(int index, vector<double> rr, MatrixXcd & ZL)
{
	int Ncols = Xcor[0].size();
//	int index = MinDis(rr, Xcor);
	if(Ncols != tR.size())throw("Wrong Nspecies.");
	if(kvector.size() != tR[0].size())throw("Wrong k-point numbers.");
	if (index > tR.size())throw("index exceed tR dimensions.");
	int Nlm=Lv.size();
	VectorXcd J0(Nlm);
	VectorXcd H0(Nlm);

	MatrixXcd tmp(Nlm,Nlm);
	VectorXcd tmpz0(Nlm);

	for(int k=0; k < kvector.size(); k++)
	{
		for(int lm=0; lm<Nlm; lm++)
		{
			for(int lm2=0; lm2<Nlm; lm2++)
			{
				tmp(lm,lm2)=tR[index][k][lm][lm2];
			}
		}
		JL(index, rr,kvector[k], J0, H0);
		tmpz0=tmp*J0+complex<double>(0.0,(-1.0)*kvector[k])*H0;
		for(int lm=0; lm<Nlm; lm++)
		{
			ZL(k,lm)=tmpz0(lm);
		}
	}
}



void GreenFunc::operator()(vector<double> rr1, vector<double> rr2, VectorXcd & GF)
{
	
	GF = VectorXcd::Zero(GF.size());
	
	int Nlm = Lv.size();
	int Nk =kvector.size();
	MatrixXcd tmpZ(Nlm,Nlm);
	ZFunc Zf(Lmax, kvector, Lv, Mv, Xcor, tR);
	MatrixXcd ZL1(Nk,Nlm);
	MatrixXcd ZL2(Nk,Nlm);
	MatrixXcd ZL3(Nk,Nlm);
	
//	int ind1 = MinDis(rr1,Xcor);
//	Zf(ind1, rr1, ZL1);
	
	vector<double> mrr2(3);
	for(int i=0;i<3;i++) mrr2[i] = (-1.0)*rr2[i];
//	int indm2 = MinDis(mrr2,Xcor);
//	Zf(indm2, mrr2, ZL2);
	
//	int ind2 = MinDis(rr2,Xcor);
//	Zf(ind2, rr2, ZL3);
	
	VectorXcd tmpZL1(Nlm);
	VectorXcd tmpZL2(Nlm);
	VectorXcd tmpZL3(Nlm);
	
	VectorXcd J0(Nlm);
	VectorXcd H0(Nlm);
	complex<double> ctmp;
	double tha;
	int ind2;
	for(int k=0;k < Nk;k++)
	{
		for(int ind1=0;ind1<Xcor[0].size();ind1++)
		{
			Zf(ind1, rr1, ZL1);
		//	for(int ind2=0;ind2<Xcor[0].size();ind2++)
		//	{
				ind2 = ind1;
				Zf(ind2, mrr2, ZL2);
				Zf(ind2, rr2, ZL3);
				for(int lm=0; lm<Nlm; lm++)
				{
					for(int lm2=0; lm2<Nlm; lm2++)
					{
						tmpZ(lm,lm2)=complex<double>(scatterR[k][ind1*Nlm+lm][ind2*Nlm+lm2],
								scatterI[k][ind1*Nlm+lm][ind2*Nlm+lm2]);
					}
				
		//		cout << "In GreesFunc.cpp succeed "<<endl;
			
					tmpZL1(lm) = ZL1(k,lm);
					tmpZL2(lm) = ZL2(k,lm);
					tmpZL3(lm) = ZL3(k,lm);
				}
		//	cout << "ZL" << ZL1(k,lm)<< endl;
		//	cout << "ZL2" << ZL2(k,lm)<< endl;
		//	cout << "ZL3" << ZL3(k,lm)<< endl;
				ctmp =	tmpZL1.adjoint()*tmpZ*tmpZL2;
				GF(k) = ctmp + GF(k);
//		cout << "GF(k)" << GF(k)<< endl;
		
				tha=norm_hw(rr1)-norm_hw(rr2);
				if(tha > 0.0)
				{
					JL(ind1, rr1, kvector[k], J0, H0);
					ctmp = tmpZL3.adjoint()*J0;
					GF(k) = ctmp + GF(k);
				}
				else if(tha < 0.0)
				{
					JL(ind2, rr2, kvector[k], J0, H0);
					ctmp = tmpZL1.adjoint()*J0;
					GF(k) = GF(k)+ctmp;
				}
		}
	}
}

//void GreenFunc::ASALDOS(vector<double> IntegRadius, 
//						int row_index, vector<double> & RHO)//Atomic sphere approx.
//{
//	int Nk = kvector.size();
//	int Nlm = Lv.size();
//	vector<vector<double> > FL(Nk, vector<double>(Nlm,0.0));
//	
//	ZFunc Zf(Lmax, kvector, Lv, Mv, Xcor, tR);
//	
//	MatrixXcd ZL(Nk,Nlm);
//	int NN = 10000;
//	vector<vector<double> > tmpR(NN, vector<double>(3,0.0));
//	double dr = IntegRadius[row_index]/(NN*1.0);
//	
//	cout <<"IntegRadius: "<< IntegRadius[row_index] << endl;
//	
//	vector<vector<vector<complex<double> > > > 
//						tmpZL(Nk,vector<vector<complex<double> > >
//							 (Nlm, vector<complex<double> >(NN,0.0)));
//	
//	for(int i=0;i< NN;i++)
//	{
//		tmpR[i][0] = 1.0*i*dr;
//		Zf(tmpR[i], ZL);
//		for(int k=0;k<Nk;k++)
//		{
//			for(int j=0;j<Nlm;j++)
//			{
//				tmpZL[k][j][i] = ZL(k,j);
//			}
//		}
//	}
//	
//	vector<double> tmpy(NN);
//	for(int k=0;k<Nk;k++)
//	{
//		for(int j=0;j<Nlm;j++)
//		{
//			for(int i=0;i<NN;i++)
//			{
//				tmpy[i] = pow(tmpR[i][0],2.0)*
//						  ( pow(tmpZL[k][j][i].real(),2.0)
//						   +pow(tmpZL[k][j][i].imag(),2.0)
//					      );
//				cout << "tmpy: " << tmpy[i] << endl;
//			}
//			FL[k][j] = simpInt(dr, tmpy);
//			cout << "FL: " << FL[k][j] << endl;
//		}
//	}
//	
//	int lm = -1;
//	for(int k=0;k<Nk;k++)
//	{
//		lm = -1;
//		for(int l=0;l<Lmax;l++)
//		{
//			for(int m=-l;m<=l;m++)
//			{
//				lm++;
//				RHO[k] = RHO[k]-(1.0)/pi*FL[k][lm]
//					*scatterI[k][row_index*Nlm+lm][row_index*Nlm+lm];
//			}
//		}
//	}
//}
