#include "KKRcoef.h"
#include "besselfrac.h"
#include "plegendre.h"
#include "constant.h"
//const Doub pi=4.0*atan(1.0);
//const extern double pi;
//static const Int Lmax=3;
// structure constants
// to get the result, we have to first calculate the Gaunt numbers.
// The process also invokes  Hankel functions(j(kr)+i*n(kr)) multiplied by spherical 
// harmonics.
// The inputing parameter R contains three elements, (r,theta,phi), for calculating spherical
// harminics.
Complex KKRcoef::StuctConst(Doub & k, VecDoub & R){
	Doub GauntNumber(Int l1,Int l2,Int l3,Int m1,Int m2,Int m3);
	Bessel besj;
	Bessel besy;
	Doub gtNum;
	Complex cc(0.0,-4*pi*k);
	Complex besh;
	Complex sum(0,0);
	Complex ci;
	int Nlm = Lv.size();
	int l3, m3;
    for(Int lm3=0;lm3<Nlm;lm3++)
	{
//  for(Int l3=0;l3<=Lmax;l3++){
//		for(Int m3=-l3;m3<=l3;m3++){
			l3 = Lv[lm3];
			m3 = Mv[lm3];
			gtNum = GauntNumber(l1,l2,l3,m1,m2,m3);
//		cout <<"GauntNumber "<< gtNum << endl;
			besh = Complex(besj.sphbesj(l3,k*R[0]),besy.sphbesy(l3,k*R[0]));
//			cout << "besh: "<< besh.real() << " +i" << besh.imag() << endl;
			ci = Complex(0.0,1.0);
			ci = pow(ci,l2+l3-l1);
			sum = sum + ci*gtNum*besh*SpHamnx(l3, m3, R[1], R[2]); //Pp 375 (B.11)
//		}
	}
	return (sum*cc);//p494
}

Complex KKRcoef::TranslationConst(Doub & k, vector<double> & R){
	Doub GauntNumber(Int l1,Int l2,Int l3,Int m1,Int m2,Int m3);
	Bessel besj;
	Doub gtNum;
	Complex cc(4*pi,0.0);
	Complex besh;
	Complex sum(0,0);
	Complex ci; 
	int Nlm = Lv.size();
	int l3, m3;
    for(Int lm3=0;lm3<Nlm;lm3++)
	{
			l3 = Lv[lm3];
			m3 = Mv[lm3];
 //   for(Int l3=0;l3<=Lmax;l3++){
//		for(Int m3=-l3;m3<=l3;m3++){
			gtNum=GauntNumber(l1,l2,l3,m1,m2,m3);
//		cout <<"GauntNumber "<< gtNum << endl;
			besh=Complex(besj.sphbesj(l3,k*R[0]),0.0);
//			cout << "besh: "<< besh.real() << " +i" << besh.imag() << endl;
			ci=Complex(0.0,1.0);
			ci=pow(ci,l2+l3-l1);
			sum=sum+ci*gtNum*besh*SpHamnx(l3,m3,R[1],R[2]); //Pp 374 (B.5)
//		}
	}
	return (sum*cc);
}


Complex KKRcoef::StuctConst(Doub & k, vector<double> & R){
	VecDoub rr(3,0.0);
	for(int i=0;i<3;i++){
		rr[i]=R[i];
	}
	Complex ss(0.0,0.0);
	ss=StuctConst(k,rr);
	return ss;
}
