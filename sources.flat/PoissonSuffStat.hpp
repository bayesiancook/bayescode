
#ifndef POISUFFSTAT_H
#define POISUFFSTAT_H

#include <map>
using namespace std;

class PoissonSuffStat	{

	public:

	PoissonSuffStat() {}
	~PoissonSuffStat() {}

	void Clear()	{
		rootcount.clear();
		paircount.clear();
		waitingtime.clear();
	}

	map<int,int> rootcount;
	map<pair<int,int>,int> paircount;
	map<int,double> waitingtime;
};


#endif

