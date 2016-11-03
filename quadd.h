#ifndef _QUADSIMP_H_
#define _QUADSIMP_H_
#include "nr3.h"
// extended Simpson's rule for vector integral over equal interval
class quadsimp
{
	public:
		quadsimp(Doub hh, VecDoub yy):h(hh), y(yy){};
		quadsimp(){};
		Doub intd();
		~quadsimp(){};
	private:
		Doub h;//x interval
		VecDoub y;
};

inline double simpInt(double h_, vector<double> yy_)
{
	VecDoub ytmp(yy_.size());
	for(int i=0;i<yy_.size();i++) ytmp[i] = yy_[i];
	quadsimp qt(h_, ytmp);
	double ss = qt.intd();
	return ss;
};




#endif
