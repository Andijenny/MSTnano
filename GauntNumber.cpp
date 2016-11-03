//#include<iostream>
//#include<cmath>
#include "nr3.h"
#include "constant.h"
// Perform the Integral of three spherical harmonics
// cf. Computer Physics Communications 25 (1982) 81-85
// edited by H.Wang
// May 17, 2016
Doub GauntNumber(Int l1,Int l2,Int l3,Int m1,Int m2,Int m3){
	Bool istriangle(Int l1, Int l2, Int l3);
	Doub SumOverGamma(Int l1,Int l2,Int l3, Int m1, Int m2);
	Doub prefactor1(Int l1, Int l2, Int l3);
	Doub prefactor2(Int l1, Int l2, Int l3, Int m1, Int m2, Int m3);
	if (m1+m2+m3!=0) return 0;
	if (!(istriangle(l1,l2,l3))){
//		std::cout<< "1 running wrong"<< std::endl;
		return 0;}
	if (((l1+l2+l3)%2==1)){
//		std::cout<< "2 running wrong"<< std::endl;
		return 0;}
	// l1, l2, l3 form a triangle with even perimeter
	Doub F3Y;
//	std::cout<<"part:"<<prefactor1(l1,l2,l3)<<" "<<prefactor2(l1,l2,l3,m1,m2,m3);
//	std::cout<<"  "<<SumOverGamma(l1,l2,l3,m1,m2)<< endl;
	F3Y=prefactor1(l1,l2,l3)*prefactor2(l1,l2,l3,m1,m2,m3)
		*SumOverGamma(l1,l2,l3,m1,m2);
	return F3Y;
}
//-----------------

Bool istriangle(Int l1, Int l2, Int l3){
	// The criterion is not valid for formming a triangle,
	// which also contains
	// l2+l3=l1 || l1+l2=l3 || l1+l3=l2
//	if (l2+l3<=l1 || l1+l2<=l3 || l1+l3<=l2) return false;
	if (l2+l3<l1 || l1+l2<l3 || l1+l3<l2) return false;
	else return true;
}	
//-----------------
//factorial of n
Int fac(Int n){return n<=1?1:fac(n-1)*n;}
//-----------------

Doub SumOverGamma(Int l1,Int l2,Int l3, Int m1, Int m2){
	Int fac(Int n);
	Int g, M, M1, M2;
	Int t[5]={l3-l1-m2,l2+m2,l1-m1,l3-l2+m1,l1+l2-l3};
	M=std::min(t[1],t[2]);
	M=std::min(M,t[4]);
	M=std::max(M,0);
	// we should also set the lowest limit of gamma so that
	// no factorial of negative numbers appear.
	// this is missed in the reference.
	if (t[0] < 0) M1=0-t[0];
	else M1=0;
	if (t[3] < 0) M2=0-t[3];
	else M2=0;
	M2=std::max(M1,M2);
//	std::cout<< "M= " << M << std::endl;
	for(Int i=0;i<4;i++){
//	std::cout<< "t["<< i <<"]= " << t[i] << std::endl;
	}
	Doub sum1 = 0.0;
	for(g=M2;g<=M;g++){
		sum1=sum1+pow(-1,g)/fac(g)/fac(t[0]+g)/fac(t[1]-g)
            /fac(t[2]-g)/fac(t[3]+g)/fac(t[4]-g);
	}
	return sum1;
}
//-----------------

Doub prefactor1(Int l1, Int l2, Int l3){
	Int fac(Int );
//	const Doub pi=4*atan(1);
//	const extern double pi;
	Doub pf;
	pf=(2*l1+1)*(2*l2+1)*(2*l3+1)/4/pi;
	pf=pow(pf,0.5);
	Int la=(l1+l2+l3)/2;//the sum is even in GauntNumber
	pf=pf*1.0*fac(la)/fac(la-l1)/fac(la-l2)/fac(la-l3)*
		fac(l2-l1+l3)*fac(l1-l2+l3)*fac(l1+l2-l3)/fac(2*la+1);
	return pf;
}
//-----------------

Doub prefactor2(Int l1, Int l2, Int l3, Int m1, Int m2, Int m3){
	Doub pf2;
	Int la=(l1+l2+l3)/2;//the sum is even in GauntNumber
//	cout<< "l1-m1:"<< fac(l1-m1)<<endl;
//	cout<< "l1+m1:"<< fac(l1+m1)<<endl;
//	cout<< "l2-m2:"<< fac(l2-m2)<<endl;
//	cout<< "l2+m2:"<< fac(l2+m2)<<endl;
//	cout<< "l3-m3:"<< fac(l3-m3)<<endl;
//	cout<< "l3+m3:"<< fac(l3+m3)<<endl;
//	cout<< "la+l3+m1-m2:"<<pow(-1,la+l3+m1-m2)<<endl;
//	Doub pf3=sqrt(fac(l2-m2)*fac(l3+m3))*sqrt(fac(l1-m1)*fac(l1+m1))*sqrt(fac(l2+m2)*fac(l3-m3));
//	cout<< "pf3 "<< pf3<<endl;
	pf2=pow(-1.0,la+l3+m1-m2)*sqrt(fac(l1-m1)*fac(l1+m1))
		*sqrt(fac(l2-m2)*fac(l2+m2))*sqrt(fac(l3-m3)*fac(l3+m3));
//	cout<< "pf2 "<< pf2<<endl;
	return pf2;
}
