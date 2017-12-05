
#ifndef PATHSUFFSTAT_H
#define PATHSUFFSTAT_H

#include "SuffStat.hpp"
#include "SubMatrix.hpp"
#include "Array.hpp"
#include "BidimArray.hpp"
#include <map>


class PathSuffStat : public SuffStat	{

	public:

	PathSuffStat() {}
	~PathSuffStat() {}

	void Clear()	{
		rootcount.clear();
		paircount.clear();
		waitingtime.clear();
	}

	void IncrementRootCount(int state)	{
		rootcount[state]++;
	}

	void IncrementPairCount(int state1, int state2)	{
		paircount[pair<int,int>(state1,state2)]++;
	}

	void AddRootCount(int state, int in)	{
		rootcount[state] += in;
	}

	void AddPairCount(int state1, int state2, int in)	{
		paircount[pair<int,int>(state1,state2)] += in;
	}
	
	void AddWaitingTime(int state, double in)	{
		waitingtime[state] += in;
	}

	void AddTo(PathSuffStat& suffstat) const {
		for (std::map<int,int>::const_iterator i = rootcount.begin(); i!= rootcount.end(); i++)	{
			suffstat.AddRootCount(i->first,i->second);
		}
		for (std::map<pair<int,int>, int>::const_iterator i = paircount.begin(); i!= paircount.end(); i++)	{
            suffstat.AddPairCount(i->first.first,i->first.second,i->second);
		}
		for (std::map<int,double>::const_iterator i = waitingtime.begin(); i!= waitingtime.end(); i++)	{
            suffstat.AddWaitingTime(i->first,i->second);
        }
	}

	int GetRootCount(int state) const {
        std::map<int,int>::const_iterator i = rootcount.find(state);
        if (i == rootcount.end())   {
            return 0;
        }
        return i->second;
	}

	int GetPairCount(int state1, int state2) const  {
        std::map<pair<int,int>,int>::const_iterator i = paircount.find(pair<int,int>(state1,state2));
        if (i == paircount.end())   {
            return 0;
        }
        return i->second;
	}

	double GetWaitingTime(int state) const	{
        std::map<int,double>::const_iterator i = waitingtime.find(state);
        if (i == waitingtime.end()) {
            return 0;
        }
        return i->second;
	}
	
	double GetLogProb(const SubMatrix& mat) const {
		double total = 0;
		auto stat = mat.GetStationary();
		for (std::map<int,int>::const_iterator i = rootcount.begin(); i!= rootcount.end(); i++)	{
			total += i->second * log(stat[i->first]);
		}
		for (std::map<int,double>::const_iterator i = waitingtime.begin(); i!= waitingtime.end(); i++)	{
			total += i->second * mat(i->first,i->first);
		}
		for (std::map<pair<int,int>, int>::const_iterator i = paircount.begin(); i!= paircount.end(); i++)	{
			total += i->second * log(mat(i->first.first, i->first.second));
		}
		return total;
	}

    const std::map<int,int>& GetRootCountMap() const {return rootcount;}
    const std::map<pair<int,int>,int>& GetPairCountMap() const {return paircount;}
    const std::map<int,double>& GetWaitingTimeMap() const {return waitingtime;}

	private:

	std::map<int,int> rootcount;
	std::map<pair<int,int>,int> paircount;
	std::map<int,double> waitingtime;
};

class PathSuffStatArray : public SimpleArray<PathSuffStat>	{

	public:

	PathSuffStatArray(int insize) : SimpleArray<PathSuffStat>(insize) {}
	~PathSuffStatArray() {}

	void Clear()	{
		for (int i=0; i<GetSize(); i++)	{
			(*this)[i].Clear();
		}
	}

	double GetLogProb(const Selector<SubMatrix>& matrixarray) const	{

		double total = 0;
		for (int i=0; i<GetSize(); i++)	{
			total += GetVal(i).GetLogProb(matrixarray.GetVal(i));
		}
		return total;
	}

	void AddToComponents(Array<PathSuffStat>& suffstatarray, const Array<int>& alloc)	const {
		for (int i=0; i<GetSize(); i++)	{
			GetVal(i).AddTo(suffstatarray[alloc.GetVal(i)]);
		}
	}
};

class PathSuffStatBidimArray : public SimpleBidimArray<PathSuffStat>	{

	public:

	PathSuffStatBidimArray(int inncol, int innrow) : SimpleBidimArray<PathSuffStat>(inncol,innrow,PathSuffStat()) {}
	~PathSuffStatBidimArray() {}

	void Clear()	{
        for (int i=0; i<this->GetNrow(); i++)  {
            for (int j=0; j<this->GetNcol(); j++)   {
                (*this)(i,j).Clear();
            }
        }
	}

    double GetLogProb(const ConstBidimArray<SubMatrix>& matrixarray) const  {
        double total = 0;
        for (int j=0; j<this->GetNcol(); j++)   {
            total += GetLogProb(j,matrixarray);
        }
        return total;
    }

    double GetLogProb(int j, const ConstBidimArray<SubMatrix>& matrixarray) const   {
        double total = 0;
        for (int i=0; i<this->GetNrow(); i++)  {
            total += GetVal(i,j).GetLogProb(matrixarray.GetVal(i,j));
        }
        return total;
    }

    double GetLogProb(int j, const vector<int>& flag, const ConstBidimArray<SubMatrix>& matrixarray) const   {
        double total = 0;
        for (int i=0; i<this->GetNrow(); i++)  {
            if (flag[i])    {
                total += GetVal(i,j).GetLogProb(matrixarray.GetVal(i,j));
            }
        }
        return total;
    }
};

#endif
