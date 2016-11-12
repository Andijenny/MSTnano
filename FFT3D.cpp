#ifdef MULTI_THREAD
	
	#include <thread>

#endif

#include "FFT3D.h"
//#include <Dense>
#include <cmath>
#include "quad2d.h"
#include "interp_2d.h"
#include "KKRcoef.h"
#include "Geometry.h"
#include "nr3.h"
#include "constant.h"
//#include <fstream>
//using namespace Eigen;
//using namespace std;


//const double FFT3D::pi = 4.0*atan(1);
//const extern double pi;

FFT3D::FFT3D(MatrixXd BasisInRealSpace_, VectorXi KP_, MatrixXd Xcor_,
			 vector<int> Lv_, vector<int> Mv_, vector<double> kvector_):
				VR(BasisInRealSpace_), Nk(KP_), Xcor(Xcor_), 
				Lv(Lv_), Mv(Mv_), kvector(kvector_)
{
	MatrixXd Kpoints(Nk.prod(),3);
	KpointSampling(Kpoints);

}
void FFT3D::GenBasisInBZ(MatrixXd & KBasis)
{
	Vector3d vx = VR.row(0);
	Vector3d vy = VR.row(1);
	Vector3d vz = VR.row(2);
	Vector3d vkx;
	Vector3d vky;
	Vector3d vkz;
	
	volumeR = vx.dot(vy.cross(vz));
	
	vkx=2.0*pi*vy.cross(vz)/volumeR;
	vky=2.0*pi*vz.cross(vx)/volumeR;
	vkz=2.0*pi*vx.cross(vy)/volumeR;

	KBasis << vkx,
		  	  vky,
		      vkz;
}

void fk(int tNi, VectorXd & coi)
{
	double tmp1;
	int j = 0;
	for(int i=1;i<=tNi;i++)
	{
		tmp1 = (2.0*i-tNi-1)*1.0/(2.0*tNi);
		if( abs(tmp1) <= 0.5 )
		{
			coi(j) = tmp1;
			j++;
		}
	}
	coi.resize(j);
}

//Monkhorst and Pack(see PRB, 13, 5188)
void FFT3D::KpointSampling(MatrixXd & Kpoints)
{
	if( Kpoints.cols() != 3) throw("wrong Kpoints size!");
	
	MatrixXd KBasis(3,3);
	GenBasisInBZ(KBasis);
	
	VectorXd u1(Nk(0));
	VectorXd u2(Nk(1));
	VectorXd u3(Nk(2));
	fk(u1.size(), u1);
	fk(u2.size(), u2);
	fk(u3.size(), u3);
	
	int cc = 0;
	for(int i=0;i<Nk(0);i++)
		for(int j=0;j<Nk(1);j++)
			for(int k=0;k<Nk(2);k++)
				Kpoints.row(cc++) = u1(i)*KBasis.row(0)+u2(j)*KBasis.row(1)
									+u3(k)*KBasis.row(2);				
		//		Kpoints.row(cc++) = -0.5*KBasis.row(0)+1.0*(i)/(Nk(0)*1.0)*KBasis.row(0)
		//						    -0.5*KBasis.row(1)+1.0*(j)/(Nk(1)*1.0)*KBasis.row(1)
		//						    -0.5*KBasis.row(2)+1.0*(k)/(Nk(2)*1.0)*KBasis.row(2);
}

/*	Fourier transform of a matrix Func_in of size N*N.Assuming the original 
 *	matrix func_in is M000,then we need to consider M00-1 M000 M001 ... M111(27 
 *	matrices in whole) under the first neighboring approximation, the subscript
 *	means the displacement along VR.
*/

void fftKKR_divideAndConquer::operator()(int kp, int i, int j, int k,
				  vector<vector<complex<double> > >& MKKRaux, int kv,
										 int MrowStart, int MrowEnd)
{
	VectorXd tmpRmn(3);
	vector<double> Xij(3), Rij(3);
	KKRcoef G1;
	double tmpPhase;
	int Nlm = Lv.size();	
	
	int Nsp = Xcor.rows();
	MatrixXd Xcor_aux(Nsp, 3);
	
	for(int m = 0; m < 3; m++)
		tmpRmn(m) = i*VR(0,m)+j*VR(1,m)+k*VR(2,m);
//	tmpRmn << i, j, k;
//	tmpRmn =  VR*tmpRmn;

	for(int l = 0;l < Xcor_aux.rows();l++)
		for(int m = 0 ; m < 3; m++)
			Xcor_aux(l,m) = Xcor(l,m) + tmpRmn(m);
	
	tmpPhase = Kpoints(kp,0)*tmpRmn(0)+
			   Kpoints(kp,1)*tmpRmn(1)+
			   Kpoints(kp,2)*tmpRmn(2);
	for(int m = MrowStart; m <= MrowEnd; m++)// cycle for Xcor
	{
		for(int n = 0; n < Nsp; n++)//cycle for Xcor_aux
		{
			if((i==0)&&(j==0)&&(k==0)&&(m != n))
			{
				for(int l=0;l<3;l++)
				{
					Xij[l] = Xcor(n,l) - Xcor(m,l);
				}
			}
			else if(!((i==0)&&(j==0)&&(k==0)))
			{
				for(int l=0;l<3;l++)
					Xij[l] = Xcor_aux(n,l) - Xcor(m,l);
			}
			if(!((n==m)&&(i==0)&&(j==0)&&(k==0)))
			{
				
				Cart2Sph(Xij, Rij);
				
				for(int lm1=0;lm1<Nlm;lm1++)
				{
					for(int lm2=0;lm2<Nlm;lm2++)
					{
						G1 = KKRcoef(Lv[lm1],Mv[lm1],Lv[lm2],Mv[lm2],Lv,Mv);
						MKKRaux[Nlm*m+lm1][Nlm*n+lm2] = 
								MKKRaux[Nlm*m+lm1][Nlm*n+lm2]+		
								complex<double>(1.0/Naux,0.0)*
								exp(complex<double>(0.0,1.0)*tmpPhase)*
								G1.StuctConst(kvector[kv],Rij);
					}
				}
			}
		}
	}
}

void FFT3D::fftKKR(vector<vector<vector<complex<double> > > >& MKKRaux, int kv)
{

	
//	vector<vector<vector<complex<double> > > > MKKRaux(NumK,
//				  vector<vector<complex<double> > >(Nsc,
//					       vector<complex<double> >(Nsc, complex<double>(0.0, 0.0)));

	
	if (Xcor.cols() != 3) throw("wrong atom coordinators");


	MatrixXd Kpoints(Nk.prod(),3);
	KpointSampling(Kpoints);
	

	int NumK = Kpoints.rows();
	int Nne = 2;// neighbouring order
	int Naux = pow(2*Nne+1,2);
	int Nsp = Xcor.rows();

	
	int num_threads=min((int)(std::thread::hardware_concurrency()), Nsp);
	int Ninv = Nsp/num_threads;
	int Nexd = Nsp%num_threads;
	fftKKR_divideAndConquer	fkParall(VR, Xcor, Lv, Mv, Kpoints, kvector, Naux);
	
	for(int kp=0;kp<NumK;kp++)
	{
		for(int i=-Nne;i<=Nne;i++)
		{
			for(int j=-Nne;j<=Nne;j++)
			{
				for(int k = 0; k < 1; k++)
				{
				#ifdef MULTI_THREAD
					std::vector<std::thread> threads(num_threads-1);
					if(Ninv > 0)
					{
						for(int n1 = 0; n1 < (num_threads-1); n1++)
							threads[n1] = std::thread(fkParall, kp, i, j, k,
													  std::ref(MKKRaux[kp]), kv, 
													  ((n1+1)*Ninv), 
													  ((n1+2)*Ninv-1));
					
						fkParall(kp, i, j, k, MKKRaux[kp], kv, 0, Ninv-1);
					
						for(int i = 0;i < (num_threads-1); i++) threads[i].join();
//						std::vector<std::thread> threads ;
//							for(int n1 = 0;n1 < num_threads;j++)
//							        threads.push_back(std::thread(fkParall, kp, i, j, k,
//													  std::ref(MKKRaux[kp]), kv, 
//													  ((n1)*Ninv), 
//													  ((n1+1)*Ninv-1)));
//
//						std::for_each(threads.begin(),threads.end(),
//							                  std::mem_fn(&std::thread::join));

					}
					if(Nexd > 0)
					{
						for(int n1 =0; n1 < Nexd; n1++)
							threads[n1] = std::thread(fkParall, kp, i, j, k,
													  std::ref(MKKRaux[kp]), kv, 
													  (Ninv*num_threads+n1), 
													  (Ninv*num_threads+n1));
					
						for(int i = 0;i < Nexd; i++) threads[i].join();
					}
				#else
					
					fkParall(kp, i, j, k, MKKRaux[kp], kv, 0, Nsp-1);
				
				#endif
				
				}
			}
		}
	}
	
}

Doub func1(const VecDoub& x1_coarse, const VecDoub& x2_coarse, const MatDoub& y_coarse, const Doub x, const Doub y)
{
	Bilin_interp myf1(x1_coarse,x2_coarse,y_coarse);
//	Spline2D_interp myf1(x1_coarse,x2_coarse,y_coarse);
	return myf1.interp(x,y);	
}

Doub y1(VecDoub BZ, const Doub x)
{
	return BZ[0];
};

Doub y2(VecDoub BZ, const Doub x)
{
	return BZ[1];
};

void ifft_divideAndConquer::operator()(vector<vector<double> >& scatterR,
									   vector<vector<double> >& scatterI,
									   int MrowStart, int MrowEnd)
{
	
	VecDoub Kpx(Nk(0)), Kpy(Nk(1));
	for(int i=0;i<Kpx.size();++i)
		Kpx[i] = Kpoints(i*Nk(1)*Nk(2),0);
	for(int i=0;i<Kpy.size();++i)
		Kpy[i] = Kpoints(i*Nk(2),1);
	
	int Nsc = MSPOaux[0].size();


	VecDoub Dky(2);
	Doub Dkx_l = KBasis(0,0)*(-0.5);
	Doub Dkx_u = KBasis(0,0)*(0.5);
	Dky[0] = KBasis(1,1)*(-0.5);
	Dky[1] = KBasis(1,1)*(0.5);
	Doub OBZ = KBasis(1,1)*KBasis(0,0);
	int ck;
	
	MatDoub tmpTijR(Kpx.size(),Kpy.size()), tmpTijI(Kpx.size(),Kpy.size());
	
	for(int i = MrowStart; i <= MrowEnd; i++)
	{
		for(int j=0;j<Nsc;j++)
		{
			ck = 0;
			for(int k1=0;k1<Kpx.size();k1++)
			{
				for(int k2=0;k2<Kpy.size();k2++)
				{
					tmpTijR[k1][k2] = MSPOaux[ck][i][j].real();
					tmpTijI[k1][k2] = MSPOaux[ck][i][j].imag();
					ck++;
				}
			}
			scatterR[i][j] = 1.0/OBZ*quad2d(func1,Kpx,Kpy,tmpTijR,Dky,Dkx_l,Dkx_u,y1,y2);
			scatterI[i][j] = 1.0/OBZ*quad2d(func1,Kpx,Kpy,tmpTijI,Dky,Dkx_l,Dkx_u,y1,y2);
		}
	}
}


void FFT3D::ifft(vector<vector<vector<complex<double> > > > MSPOaux,
				 vector<vector<double> >& scatterR,
				 vector<vector<double> >& scatterI)
{
	
	MatrixXd Kpoints(Nk.prod(),3);
	KpointSampling(Kpoints);
	
	MatrixXd KBasis(3,3);
	GenBasisInBZ(KBasis);

	int Nsc = MSPOaux[0].size();
	
	ifft_divideAndConquer iftParall(MSPOaux, Kpoints, KBasis, Nk);
	
	#ifdef MULTI_THREAD
		
		int num_threads=min((int)(std::thread::hardware_concurrency()), Nsc);
		int Ninv = Nsc/num_threads;
		int Nexd = Nsc%num_threads;
		std::vector<std::thread> threads(num_threads-1);
		
		if(Ninv > 0)
		{
			for(int n1 =0; n1 < (num_threads-1); n1++)
				threads[n1] = std::thread(iftParall, std::ref(scatterR), 
										  std::ref(scatterI),
										  ((n1+1)*Ninv), 
										  ((n1+2)*Ninv-1));
					
			iftParall(scatterR, scatterI, 0, Ninv-1);
					
			for(int i = 0;i < (num_threads-1); i++) threads[i].join();
		}
		if(Nexd > 0)
		{
			for(int n1 =0; n1 < Nexd; n1++)
				threads[n1] = std::thread(iftParall, std::ref(scatterR), 
										  std::ref(scatterI),
										  (Ninv*num_threads+n1), 
										  (Ninv*num_threads+n1));
					
			for(int i = 0;i < Nexd; i++) threads[i].join();
		}
		
	#else
			
		iftParall(scatterR, scatterI, 0, Nsp-1);
					
	#endif

}


