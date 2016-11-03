#include "Tmatrix2.h"
#include "quadd.h"
#include "besselfrac.h"
#include <iomanip>
#include <fstream>
#include<sstream>


double Tmatrix::DlNext(vector<double> Rl)
{
	Bessel besj;
	int Nr = radius.size();
	vector<double> Intmp(Nr,0.0);
	
	if(Nr != Rl.size())
		throw("Rl has wrong size.");
	double j0;
	for(int i=0; i < Nr; i++)
	{
		j0 = besj.sphbesj(LL, kvalue*radius[i]);
		Intmp[i] = j0*potential[i]*Rl[i]*radius[i]*radius[i];
	//	cout << "intmp[i]" << Intmp[i]<<endl;
	}
	double h = radius[1]-radius[0];
	double ss = simpInt(h, Intmp);
	ss = (-1.0)*kvalue*ss;
	return ss;
}


void Tmatrix::RlNext(vector<double>& Rl, double Dl)
{
	Bessel besj;
	Bessel besy;
	double j0; 
	double y0; 
	for(int i=0;i<Rl.size();i++)
	{
		j0 = besj.sphbesj(LL, kvalue*radius[i]);
		y0 = besy.sphbesy(LL, kvalue*radius[i]);
		Rl[i] = j0-y0*Dl;
	}
};



complex<double> Tmatrix::calc_tma(double MixWeight, double ConvergPrecision)
{
	vector<double> Rl(radius.size(),0.0);
	double t0 = 10.0;
	double t1 = 0.0;
	double dtt = abs(t0-t1);
	int cc=0;
	while( dtt > ConvergPrecision)
	{
		t0 = t1;
		RlNext(Rl, t1);
		t1 = DlNext(Rl);
		dtt=abs(t0-t1);
		cout << "step: "<< ++cc << "	t1" << t1 <<"		dt"<< dtt << endl;
	}
	complex<double> ss;
	double atm = atan(t1);
//	atm = -1.39;
	ss = (-1.0)*complex<double>(cos(atm), sin(atm))*sin(atm)/(kvalue);
	
	cout<<"For angular quantum number l= "<< LL << endl;
	cout << "converged at: \n";
	cout<<" kvector "<<"		tLL'		  "<<"	PhastShift		"<< endl; 
	cout << setw(10)<<setprecision(6)<<kvalue<<"		"<<ss<<"		"<<atm<<endl;
	return ss;
}

//#endif	
