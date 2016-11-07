#include "quadd.h"

Doub quadsimp::intd(){
	Doub s;
	Int n=y.size();
//	cout <<"n: "<< n<<endl;
	VecDoub co(n);
	if(n==0){throw("Empty vector!");}
	else if(n==1) {s=0.0;}
	else if(n==2) {s=h*(y[0]+y[1])/2.0;}
	else if(n==3) {s=h*(1.0/3.0*y[0]+4.0/3.0*y[1]+1.0/3.0*y[2]);}
	else{
		co[0]=co[n-1]=1.0/3.0;
		for(Int i=1;i<Int((n-2)/2)+1;i++){
			co[2*i-1]=4.0/3.0;
			co[2*i]=2.0/3.0;
		}
		if (n & 1) co[n-2]=4.0/3.0; //n is odd
		//for(Int i=0;i<n;i++)cout << co[i] << endl;
	    co=co*y;// vector mulitiplication, defined in nr3.h
		//for(Int i=0;i<n;i++)cout << co[i] << endl;
		s=h*co.sum();
	}
	return s;
}
