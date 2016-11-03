#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_
//#include "nr3.h"
#include<iostream>
#include<vector>
#include<cmath>
#include <Eigen>

using namespace std;
using namespace Eigen;

class Geom{

public:
	Geom(char * fname):fn(fname){};
	~Geom(){};
	void Rcoor(vector<vector<double> > & rr, vector<string> & SA,
			   vector<int> & SL, vector<double> & IntegRadius);
	double MinDis2(int index, vector<vector<double> > Xcor);
private:
	char * fn;
};

inline double norm_hw(vector<double> & r){
	return (sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]));
};

void Cart2Sph(vector<double> & x, vector<double> & r);

void Cart2Sph(VectorXd & x, VectorXd & r);

#endif
