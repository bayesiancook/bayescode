#pragma once
#include "PathSuffStat.hpp"


class MeanPathSuffStat {
  public:
    MeanPathSuffStat(int inNstate) : 
	    Nstate(inNstate),
	    rootcount(Nstate,0),
	    paircount(Nstate*Nstate,0),
	    waitingtime(Nstate,0) {}

    ~MeanPathSuffStat() {}

    int GetNstate() const {
        return Nstate;
    }

    void Add(const PathSuffStat &suffstat) {
        for (std::map<int, double>::const_iterator i = suffstat.GetRootCountMap().begin();
             i != suffstat.GetRootCountMap().end(); i++) {
            rootcount[i->first] += i->second;
        }
        for (std::map<pair<int, int>, double>::const_iterator i = suffstat.GetPairCountMap().begin();
             i != suffstat.GetPairCountMap().end(); i++) {
            paircount[i->first.first * Nstate + i->first.second] += i->second;
        }
        for (std::map<int, double>::const_iterator i = suffstat.GetWaitingTimeMap().begin();
             i != suffstat.GetWaitingTimeMap().end(); i++) {
            waitingtime[i->first] += i->second;
        }
    }

    void Normalize(double f)    {
        for (auto x : rootcount) {
            x *= f;
        }
        for (auto x : paircount)	{
            x *= f;
        }
        for (auto x : waitingtime)	{
            x *= f;
        }
    }

    const vector<float>& GetPairCounts() const	{
	    return paircount;
    }

    const vector<float>& GetWaitingTimes() const {
	    return waitingtime;
    }

    const vector<float>& GetRootCounts() const	{
	    return rootcount;
    }

    vector<float> GetPostProbs() const	{
	    float tot = 0;
	    for (int i=0; i<Nstate; i++)	{
		    tot += waitingtime[i];
	    }
	    vector<float> v(Nstate,0);
        if (tot > 0)    {
            for (int i=0; i<Nstate; i++)	{
                v[i] = waitingtime[i] / tot;
                if (isnan(v[i]))    {
                    cerr << "error nan in post probs\n";
                    exit(1);
                }
                if (isinf(v[i]))    {
                    cerr << "error inf in post probs\n";
                    exit(1);
                }
            }
        }
	    return v;
    }

    vector<float> GetWeightedWaitingTimes(const double* weight) const	{
	    vector<float> v(Nstate*Nstate,0);
	    for (int i=0; i<Nstate; i++)	{
		    for (int j=0; j<Nstate; j++)	{
			    if (i != j)	{
				    v[i*Nstate + j] = waitingtime[i] * weight[j];
                    if (isnan(v[i*Nstate + j]))    {
                        cerr << "error nan in post probs\n";
                        exit(1);
                    }
                    if (isinf(v[i*Nstate + j]))    {
                        cerr << "error inf in post probs\n";
                        exit(1);
                    }
			    }
		    }
	    }
	    return v;
    }

  private:
    int Nstate;
    std::vector<float> rootcount;
    std::vector<float> paircount;
    std::vector<float> waitingtime;
};

class MeanPathSuffStatBidimArray : public SimpleBidimArray<MeanPathSuffStat> {
  public:
    MeanPathSuffStatBidimArray(int inncol, int innrow, int inNstate)
        : SimpleBidimArray<MeanPathSuffStat>(inncol, innrow, MeanPathSuffStat(inNstate)), Nstate(inNstate) {}
    ~MeanPathSuffStatBidimArray() {}

    //! reduce path suffstats based on bidim allocations
    void Add(const PathSuffStatBidimArray& suffstatarray)	{
        for (int i=0; i<this->GetNrow(); i++) {
            for (int j=0; j<this->GetNcol(); j++) {
                (*this)(i,j).Add(suffstatarray.GetVal(i,j));
            }
        }
    }

    void Normalize(double f) {
        for (int i=0; i<this->GetNrow(); i++) {
            for (int j=0; j<this->GetNcol(); j++) {
                (*this)(i, j).Normalize(f);
            }
        }
    }

    vector<float> GetAllPostProbs() const	{
	    vector<float> v;
	    v.reserve(this->GetNrow() * this->GetNcol() * Nstate);
	    for (int i=0; i<this->GetNrow(); i++) {
		    for (int j=0; j<this->GetNcol(); j++) {
			    auto w = this->GetVal(i,j).GetPostProbs();
                for (auto x : w)    {
                    if (isnan(x))   {
                        cerr << "nan in all post probs\n";
                        exit(1);
                    }
                    if (isinf(x))   {
                        cerr << "inf in all post probs\n";
                        exit(1);
                    }
                    if (x<0)   {
                        cerr << "negative post prob\n";
                        exit(1);
                    }
                }
			    v.insert(v.end(), w.begin(), w.end());
		    }
	    }
	    return v;
    }

    vector<float> GetAllWeightedWaitingTimes(const double* weight)	const {
	    vector<float> v;
	    v.reserve(this->GetNrow() * this->GetNcol() * Nstate * Nstate);
	    for (int i=0; i<this->GetNrow(); i++) {
		    for (int j=0; j<this->GetNcol(); j++) {
			    auto w = this->GetVal(i,j).GetWeightedWaitingTimes(weight);
                for (auto x : w)    {
                    if (isnan(x))   {
                        cerr << "nan in all weighted waiting times\n";
                        exit(1);
                    }
                    if (isinf(x))   {
                        cerr << "inf in all weighted waiting times\n";
                        exit(1);
                    }
                    if (x<0)   {
                        cerr << "negative weighted waiting time\n";
                        exit(1);
                    }
                }
			    v.insert(v.end(), w.begin(), w.end());
		    }
	    }
	    return v;
    }

    vector<float> GetAllPairCounts()	const {
	    vector<float> v;
	    v.reserve(this->GetNrow() * this->GetNcol() * Nstate * Nstate);
	    for (int i=0; i<this->GetNrow(); i++) {
		    for (int j=0; j<this->GetNcol(); j++) {
			    auto w = this->GetVal(i,j).GetPairCounts();
                for (auto x : w)    {
                    if (isnan(x))   {
                        cerr << "nan in all pair counts\n";
                        exit(1);
                    }
                    if (isinf(x))   {
                        cerr << "inf in all pair counts\n";
                        exit(1);
                    }
                    if (x<0)   {
                        cerr << "negative pair count\n";
                        exit(1);
                    }
                }
			    v.insert(v.end(), w.begin(), w.end());
		    }
	    }
	    return v;
    }

    vector<float> GetAllRootCounts()	const {
	    vector<float> v;
	    v.reserve(this->GetNcol() * Nstate);
        for (int j=0; j<this->GetNcol(); j++) {
            auto w = this->GetVal(0,j).GetRootCounts();
            for (auto x : w)    {
                if (isnan(x))   {
                    cerr << "nan in all root counts\n";
                    exit(1);
                }
                if (isinf(x))   {
                    cerr << "inf in all root counts\n";
                    exit(1);
                }
                if (x<0)   {
                    cerr << "negative root count\n";
                    exit(1);
                }
            }
            v.insert(v.end(), w.begin(), w.end());
        }
	    return v;
    }

  private:
    int Nstate;
};
