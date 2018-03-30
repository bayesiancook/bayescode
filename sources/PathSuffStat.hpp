
#ifndef PATHSUFFSTAT_H
#define PATHSUFFSTAT_H

#include "SuffStat.hpp"
#include "SubMatrix.hpp"
#include "Array.hpp"
#include "BidimArray.hpp"
#include "BranchArray.hpp"
#include "NodeArray.hpp"
#include <map>
// #include "BranchSitePath.hpp"
#include "PhyloProcess.hpp"


/**
 * \brief A general sufficient statistic for substitution histories, as a function of the substitution rate matrix.
 *
 * The probability of a set of detailed substitution histories (collectively denoted as S),
 * as a function of some rate matrix Q = (Q_ab), with equilibrium frequencies pi = (pi_a), can be written as:
 *
 * p(S | Q) propto (prod_a pi_a^u_a) (prod_a exp(t_a Q_aa)) (prod_ab Q_ab^v_ab),
 *
 * where u_a is the total number of times state a was seen at the root (root count statistic),
 * v_ab (pair is the total number of substitution events from a to b (pair count stat),
 * and t_a is the total waiting time in state a (waiting time stat) --
 * all this, across all substitution histories included in S.
 *
 * PathSuffStat implements this idea, by providing methods for gathering sufficient statistics across substitution histories
 * (see also BranchSitePath::AddPathSuffStat), adding them across sites and/or branches,
 * and calculating the log p(S | Q) for any matrix Q
 *
 * These path suffstats can be used for any Markovian substitution process (any Q).
 * In some cases (i.e. for Muse and Gaut codon models),
 * they can be further simplified, as a function of the nucleotide rate parameters or the omega parameter of the Q matrix,
 * leading to even more compact suff stats (see OmegaPathSuffStat and NucPathSuffStat).
 *
 * In terms of implementation, these suffstats are encoded as sparse data structures
 * (since a very small subset of all possible pairs of codons will typcially be visited
 * by the substitution history of a given site, for instance).
 * This sparse encoding is crucial for efficiency (both in terms of time and in terms of RAM usage).
 */

class PathSuffStat : public SuffStat	{

	public:

	PathSuffStat() {}
	~PathSuffStat() {}

    //! set suff stats to 0
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

    //! add path sufficient statistics from PhyloProcess (site-homogeneous case)
    void AddSuffStat(const PhyloProcess& process)   {
        process.AddPathSuffStat(*this);
    }

    void Add(const PathSuffStat& suffstat)  {
		for (std::map<int,int>::const_iterator i = suffstat.GetRootCountMap().begin(); i!= suffstat.GetRootCountMap().end(); i++)	{
			AddRootCount(i->first,i->second);
		}
		for (std::map<pair<int,int>, int>::const_iterator i = suffstat.GetPairCountMap().begin(); i!= suffstat.GetPairCountMap().end(); i++)	{
            AddPairCount(i->first.first,i->first.second,i->second);
		}
		for (std::map<int,double>::const_iterator i = suffstat.GetWaitingTimeMap().begin(); i!= suffstat.GetWaitingTimeMap().end(); i++)	{
            AddWaitingTime(i->first,i->second);
        }
	}


    PathSuffStat& operator+=(const PathSuffStat& from)  {
        Add(from);
        return *this;
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

    //! return log p(S | Q) as a function of the Q matrix given as the argument
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

    //! const access to the ordered map giving the root count stat (sparse data structure)
    const std::map<int,int>& GetRootCountMap() const {return rootcount;}
    //! const access to the ordered map giving the pair count stat (sparse data structure)
    const std::map<pair<int,int>,int>& GetPairCountMap() const {return paircount;}
    //! const access to the ordered map giving the waiting time stat (sparse data structure)
    const std::map<int,double>& GetWaitingTimeMap() const {return waitingtime;}

	private:

	std::map<int,int> rootcount;
	std::map<pair<int,int>,int> paircount;
	std::map<int,double> waitingtime;
};


/**
 * \brief An array of substitution path sufficient statistics
 *
 * This class provides an interface for dealing with cases where
 * each item (each site, or each component of a mixture) has a different rate matrix Q_i
 */

class PathSuffStatArray : public SimpleArray<PathSuffStat>	{

	public:

	PathSuffStatArray(int insize) : SimpleArray<PathSuffStat>(insize) {}
	~PathSuffStatArray() {}

    //! set all suff stats to 0
	void Clear()	{
		for (int i=0; i<GetSize(); i++)	{
			(*this)[i].Clear();
		}
	}

    //! add path sufficient statistics from PhyloProcess (site-heterogeneous case)
    void AddSuffStat(const PhyloProcess& process)   {
        process.AddPathSuffStat(*this);
    }

    //! return total log prob (summed over all items), given an array of rate matrices
	double GetLogProb(const Selector<SubMatrix>& matrixarray) const	{

		double total = 0;
		for (int i=0; i<GetSize(); i++)	{
			total += GetVal(i).GetLogProb(matrixarray.GetVal(i));
		}
		return total;
	}

    //! \brief add suffstatarray given as argument to this array based on the allocations provided as the second argument (mixture models)
    //!
    //! specifically, for each i=0..GetSize()-1, (*this)[alloc[i]] += suffstatarray[i]
    void Add(const Selector<PathSuffStat>& suffstatarray, const Selector<int>& alloc)   {
        for (int i=0; i<suffstatarray.GetSize(); i++) {
            (*this)[alloc.GetVal(i)] += suffstatarray.GetVal(i);
        }
    }
};

/**
 * \brief A NodeArray of substitution path sufficient statistics
 *
 * This class provides an interface for dealing with cases where
 * each branch has a different rate matrix Q_j
 */

class PathSuffStatNodeArray : public SimpleNodeArray<PathSuffStat>	{

	public:

	PathSuffStatNodeArray(const Tree& intree) : SimpleNodeArray<PathSuffStat>(intree) {}
	~PathSuffStatNodeArray() {}

    //! set all suff stats to 0
	void Clear()	{
		for (int i=0; i<GetNnode(); i++)	{
			(*this)[i].Clear();
		}
	}

    //! add path sufficient statistics from PhyloProcess (site-heterogeneous case)
    void AddSuffStat(const PhyloProcess& process)   {
        process.AddPathSuffStat(*this);
    }

    //! return total log prob (summed over all items), given an array of rate matrices
	double GetLogProb(const BranchSelector<SubMatrix>& matrixarray, const SubMatrix& rootmatrix) const	{

        double ret = RecursiveGetLogProb(GetTree().GetRoot(),matrixarray,rootmatrix);
        return ret;
    }

    double RecursiveGetLogProb(const Link* from, const BranchSelector<SubMatrix>& matrixarray, const SubMatrix& rootmatrix) const   {

        double total = 0;
        if (from->isRoot()) {
            total += GetVal(from->GetNode()->GetIndex()).GetLogProb(rootmatrix);
        }
        else    {
            total += GetVal(from->GetNode()->GetIndex()).GetLogProb(matrixarray.GetVal(from->GetBranch()->GetIndex()));
        }
        for (const Link* link=from->Next(); link!=from; link=link->Next())  {
            total += RecursiveGetLogProb(link->Out(),matrixarray,rootmatrix);
        }
		return total;
	}
};


/**
 * \brief A bi-dimensional array of substitution path sufficient statistics
 *
 * This class provides an interface for dealing with cases where
 * each site/condition pair has a different rate matrix Q_i
 * (essentially, the DiffSelModel).
 */

class PathSuffStatBidimArray : public SimpleBidimArray<PathSuffStat>	{

	public:

	PathSuffStatBidimArray(int inncol, int innrow) : SimpleBidimArray<PathSuffStat>(inncol,innrow,PathSuffStat()) {}
	~PathSuffStatBidimArray() {}

    //! set all suff stats to 0
	void Clear()	{
        for (int i=0; i<this->GetNrow(); i++)  {
            for (int j=0; j<this->GetNcol(); j++)   {
                (*this)(i,j).Clear();
            }
        }
	}

    //! add path sufficient statistics from PhyloProcess (site- and condition-heterogeneous case)
    void AddSuffStat(const PhyloProcess& process, const BranchAllocationSystem& branchalloc)   {
        process.AddPathSuffStat(*this,branchalloc);
    }


		    //! \brief add suffstatarray given as argument to this array based on the allocations provided as the second argument (mixture models)
		    //!
		    //! specifically, for each i=0..GetSize()-1, (*this)[alloc[i]] += suffstatarray[i]
				void Add(const BidimSelector<PathSuffStat>& suffstatarray, const Selector<int>& alloc)   {
					for (int i=0; i<this->GetNcol(); i++) {
						for (int j=0; j<this->GetNrow(); j++)	{
							(*this)(j,alloc.GetVal(i)) += suffstatarray.GetVal(j,i);
						}
					}
				}



    //! return total log prob (summed over all items), given a bi-dimensional array of rate matrices
    double GetLogProb(const BidimSelector<SubMatrix>& matrixarray) const  {
        double total = 0;
        for (int j=0; j<this->GetNcol(); j++)   {
            total += GetLogProb(j,matrixarray);
        }
        return total;
    }

    //! return log prob summed over a given column
    double GetLogProb(int j, const BidimSelector<SubMatrix>& matrixarray) const   {
        double total = 0;
        for (int i=0; i<this->GetNrow(); i++)  {
            total += GetVal(i,j).GetLogProb(matrixarray.GetVal(i,j));
        }
        return total;
    }

    //! return log prob summed over a given column (and only for items for which flag is non 0)
    double GetLogProb(int j, const vector<int>& flag, const BidimSelector<SubMatrix>& matrixarray) const   {
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
