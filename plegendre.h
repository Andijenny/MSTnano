inline Doub plegendre(const Int l, const Int m, const Doub x) {
	static const Doub PI=3.141592653589793;
	Int i,ll;
	Doub fact,oldfact,pll,pmm,pmmp1,omx2;
	if (m < 0 || m > l || abs(x) > 1.0)
		throw("Bad arguments in routine plgndr");
	pmm=1.0;
	if (m > 0) {
		omx2=(1.0-x)*(1.0+x);
		fact=1.0;
		for (i=1;i<=m;i++) {
			pmm *= omx2*fact/(fact+1.0);
			fact += 2.0;
		}
	}
	pmm=sqrt((2*m+1)*pmm/(4.0*PI));
	if (m & 1)
		pmm=-pmm;
	if (l == m)
		return pmm;
	else {
		pmmp1=x*sqrt(2.0*m+3.0)*pmm;
		if (l == (m+1))
			return pmmp1;
		else {
			oldfact=sqrt(2.0*m+3.0);
			for (ll=m+2;ll<=l;ll++) {
				fact=sqrt((4.0*ll*ll-1.0)/(ll*ll-m*m));
				pll=(x*pmmp1-pmm/oldfact)*fact;
				oldfact=fact;
				pmm=pmmp1;
				pmmp1=pll;
			}
			return pll;
		}
	}
}
inline Doub plgndr(const Int l, const Int m, const Doub x)
{
	const Doub PI=3.141592653589793238;
	if (m < 0 || m > l || abs(x) > 1.0){
		cout << m <<"  "<< l <<"  "<< abs(x) <<endl;
		throw("Bad arguments in routine plgndr");
	}
	Doub prod=1.0;
	for (Int j=l-m+1;j<=l+m;j++)
		prod *= j;
	return sqrt(4.0*PI*prod/(2*l+1))*plegendre(l,m,x);
}

inline Int Factorial(Int n){ return (n>1?Factorial(n-1)*n:1);}


inline Complex SpHamnx(const Int l, const Int m, Doub theta, Doub phi){
	const Doub PI=3.141592653589793238;
	if(l<m) throw("bad arguments in routine SpHamnx");
	Complex epi(cos(m*phi),sin(m*phi));
//	cout << "epi: "<<epi <<endl;
	Complex re;
	re=1.0*sqrt((2.0*l+1.0)/4.0/PI*Factorial(l-abs(m))/Factorial(l+abs(m)))*plgndr(l,abs(m),cos(theta))*epi;
	if(m<0) re=pow(-1,abs(m))*Complex(re.real(),-re.imag());
	return re;
}


