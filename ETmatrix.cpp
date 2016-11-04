//#ifndef _TMATRIX_CPP_
//#define _TMATRIX_CPP_//for template class not supporting the seperation of defination and implement. 
#include "ETmatrix.h"
#include "quadd.h"
#include "besselfrac.h"
#include <iomanip>
#include <fstream>
#include <ctime>
Doub ulscf::quadd(){
	if (qd1.size() != qd2.size())throw("different range of quadrature!");
	VecDoub y1(r.size()),y2(r.size());
	for(Int i=0;i<r.size();i++){
		y1[i]=r[i]*besjr[i]*pot[i]*ur[i];
		y2[i]=r[i]*besnr[i]*pot[i]*ur[i];
	}
	Doub h=r[1]-r[0];
	VecDoub vv1;
	VecDoub vv2;
	quadsimp qs1;
	quadsimp qs2;
	quadsimp qs3(h,y1);
	for(Int i=0;i<r.size();i++){
		vv1=y1;
		vv2=y2;
		vv1.resize(0,i);
		vv2.resize(i,r.size()-1);
		qs1= quadsimp(h,vv1);
		qs2= quadsimp(h,vv2);
		qd1[i]=qs1.intd();
		qd2[i]=qs2.intd();
	}
	qd3=qs3.intd();
	
	return qd3;
}
		

void ulscf::next(const Doub & weight){
	for(Int i=0;i<r.size();i++){
		ur[i]=ur[i]+weight*(r[i]*(besjr[i]*cos(dt)+k*besnr[i]*qd1[i]
			+k*besjr[i]*qd2[i])-ur[i]);
	}
}

Doub ulscf::delta_new(){
	Doub t1=-k*qd3;
	if(abs(t1)>1.0) 
	{
		cout << "Warning !" << endl;
		cout << "dt "<< dt << endl; 
		cout << "t1 "<< t1 << endl;
		cout << "k "<< k <<endl;
		cout << "qd3 " << qd3 << endl;
		srand(time(0));
		t1=1.0-0.05*(rand()%11/10.0);
	}
	dt=asin(t1);
	 return dt;
}


complex<double> Tmatrix::calc_tma(double MixWeight, double ConvergPrecision)
{
	
	VecDoub u0(radius.size());
	VecDoub j0(radius.size());
	VecDoub y0(radius.size());
	VecDoub qd10(radius.size());
	VecDoub qd20(radius.size());
	Bessel besj;
	Bessel besy;
	ulscf ul;
	Doub t0;
	Doub t1 = 10.0;
	Int cc;
	Doub tLL = 0.0;
	for(int i = 0; i < radius.size(); i++)
	{
		j0[i] = besj.sphbesj(LL,radius[i]*kvalue);
		y0[i] = besy.sphbesy(LL,radius[i]*kvalue);
		u0[i] = radius[i]*cos(t1)*j0[i];
	}
	VecDoub CopyRadius(radius.size(), 0.0);
	VecDoub CopyPotential(potential.size(), 0.0);
	for(int i = 0; i < radius.size(); i++)
	{
		CopyRadius[i] = radius[i];
		CopyPotential[i] = potential[i];
	}
	
    ul = ulscf(LL, CopyRadius, u0, CopyPotential, qd10,
			  qd20, tLL, j0, y0, kvalue, t1);	
	ul.quadd();
	t0 = ul.delta_new();
	Doub dtt = abs(t0-t1);
	cc = 0;
	while((dtt > ConvergPrecision))
	{
		cc++;
		t0 = t1;
		ul.next(MixWeight);
		tLL = ul.quadd();
		t1 = ul.delta_new();
		dtt = abs(t0-t1);
	}
	
	cout<<"For angular quantum number l= "<< LL << endl;
	cout << "converged at: \n";
	cout<<" kvector "<<"		tLL'		  "<<"	PhastShift		"<< endl; 
	cout << setw(10)<<setprecision(6)<<kvalue;
	cout <<"		"<< tLL <<"		"<< t1 <<endl;
	return complex<double>(tLL, 0.0);
}


//#endif	
