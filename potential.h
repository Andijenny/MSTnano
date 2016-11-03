#ifndef _POTENTIAL_H_
#define _POTENTIAL_H_
#include "nr3.h"

//extern char * potfile;

class pot
{
protected:
	VecDoub radius;
	VecDoub potential;
	Int n;
public:
	pot(){};
	pot(const char * fname);
	pot(pot & p, const VecDoub & r);
	int length(){return n;};
	Doub rmax(){return radius.max();};
	void rv(vector<double> & rr);
	void pv(vector<double> & pp);
	~pot(){};
};

class potfit : public pot
{
public:
	potfit(pot & pp): pot(pp){};//constructor
	Doub operator()(const Doub r);
};

#endif
