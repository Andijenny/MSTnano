#include "potential.h"
#include "interp_1d.h"
#include<sstream>


//constructor of class pot
//read potential(r) from outer readable file
//and return actual size of data.
pot::pot(const char * fname)
{
	ifstream fin;
	fin.open(fname);
	VecDoub rr(10000,0.0);
	VecDoub U(10000,0.0);
	string buf;
	Int i=0;
	while(getline(fin,buf)){
		istringstream instr(buf);
		instr >> rr[i] >> U[i] ;
//		cout << rr[i] <<" "<<U[i] <<endl;
		i++;
	}
	fin.close();
	rr.resize(0,i-1);
	U.resize(0,i-1);
	radius=rr;
	potential=U;
	n=i;
}

//constructor for interpolation
//pot & p : already known potential(r)
//VecDoub & r: interpolating point
//store interpolationg result
pot::pot(pot & p, const VecDoub & r)
{   
	VecDoub pfit(r.size(),0.0);
	Poly_interp fpot(radius,potential,4);
	for(Int i=0;i<r.size();i++){
		if (radius.max()> r[i]) pfit[i]=fpot.interp(r[i]);
	}
	radius=r;
	potential=pfit;
	n=r.size();
}

//derived class of pot class
//overloaded of ()
//so-called functor
//for the case that function as a input paremeter
Doub potfit::operator()(const Doub r)
{
//	for(int i=0;i<length();i++)
//		cout << rv()[i] << "  "<< pv()[i]<< endl;
	Poly_interp fpot(radius,potential,4);
	if (rmax()> r){ return fpot.interp(r);}
	else {return 0.0;}
}

void pot::rv(vector<double> & rr){
	int i;
	for(i=0;i<radius.size();i++)
		rr[i]=radius[i];
	rr.resize(i);
}

void pot::pv(vector<double> & pp){
	int i;
	for(i=0;i<potential.size();i++)
		pp[i]=potential[i];
	pp.resize(i);
}
