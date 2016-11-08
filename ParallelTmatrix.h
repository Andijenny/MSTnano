#ifdef MULTI_THREAD
	#include <thread>
#endif
#include "ConcurrencyInit.h"
#include "nr3.h"

void ParallelTmatrix(vector<int> ParallelPar, vector<vector<double> >rfit, vector<vector<double> >pfit,
					 vector<string>SA, int Lmax, double mw, double cp, VecDoub kvector, 
					 vector<vector<vector<double> > >& Tll)
{

	TmatrixConcurrency t_multi(rfit, pfit, SA, 
							   Lmax, mw, cp);
	#ifdef MULTI_THREAD
		int num_threads = ParallelPar[0];
		int Ncycle = ParallelPar[1];
		int Nexceed = ParallelPar[2];
	
		
		for(int j=0;j<Ncycle;j++)
		{
			std::vector<std::thread> threads(num_threads-1);
			for(int i=0;i<(num_threads-1);i++)
			{
				threads[i]=std::thread(t_multi, kvector[(i+1)+j*num_threads],
						               std::ref(Tll[(i+1)+j*num_threads]));
			}
	//master_thread:		
			t_multi(kvector[j*num_threads], Tll[j*num_threads]);
			
			for(int i=0;i<(num_threads-1);i++) threads[i].join();
		}
		
		if(Nexceed > 0)
		{
			std::vector<std::thread> threads_d(Nexceed);
		  	for(int i=0;i<(Nexceed);i++)
		  	{
		  		threads_d[i]=std::thread(t_multi, kvector[i+Ncycle*num_threads],
		  							   std::ref(Tll[Ncycle*num_threads+i]));
			}
		  	for(int i=0;i<(Nexceed);i++)
		  	{
				threads_d[i].join();
			}
		}
	#else
		for(int k=0;k<Tll.size();k++)
			t_multi(kvector[k],Tll[k]);

	#endif
}


void const TmatrixPrint(char* FileTll, const vector<vector<vector<double> > > Tll)
{
	int cc = 1;
	ofstream ftll;
//	ftll.open("Tll_matrix");
	ftll.open(FileTll);
    for(int i=0;i<Tll.size();i++)
	{
		for(int j=0;j<Tll[0].size();j++)
		{
			for(int k=0;k<Tll[0][0].size();k++)
			{
				if(cc > 8) 
				{
					ftll << endl;
					cc = 1;
				}
				ftll << Tll[i][j][k]<<"   ";
				cc++;
			}
		}
	}
    ftll.close();
}

void TmatrixRead(char* FileTll, vector<vector<vector<double> > >& Tll)
{
	ifstream ftll;
//	ftll.open("Tll_matrix");
	ftll.open(FileTll);
	string buf;
	int ts = Tll.size()*Tll[0].size()*Tll[0][0].size();
	vector<double> tmp(ts);
	
	int cc=0;
	while(getline(ftll,buf))
	{
		istringstream instr(buf);
		for(int i=0;i<8;i++)
		{
			if(cc >= ts) break;
			instr >> tmp[cc];
			cc++;
		}
	}
	if(cc != ts)
	{   
		cout << "The total number read from Tll_matrix is : " << cc << endl;
		cout << "( != " << ts << " )" << endl;
		throw("error: Tll has wrong size.\n");
	}
	
	cc = 0;
    for(int i=0;i<Tll.size();i++)
	{
		for(int j=0;j<Tll[0].size();j++)
		{
			for(int k=0;k<Tll[0][0].size();k++)
			{
				Tll[i][j][k] = tmp[cc];
				cc++;
			}
		}
	}
    ftll.close();
}

