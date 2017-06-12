
#ifndef CODONSUFFSTAT_H
#define CODONSUFFSTAT_H

#include "PathSuffStat.hpp"
#include "CodonSubMatrixArray.hpp"
#include <typeinfo>

/*
class NucPathSuffStat : public PathSuffStat	{

	public:
	NucPathSuffStat() : PathSuffStat() {}
	~NucPathSuffStat() {}

	// assumes pathsuffstat is 61x61
	// collect the 4x4 path suff stat out of codonpathsuffstat
	void AddSuffStat(const MGCodonSubMatrix& codonstatespace, const PathSuffStat& codonpathsuffstat);
	void AddSuffStat(const MGCodonSubMatrix& codonstatespace, const PathSuffStatArray& codonpathsuffstat);
};
*/

class OmegaSuffStat {

	public:

	OmegaSuffStat() {}
	~OmegaSuffStat() {}

	// assumes pathsuffstat is 61x61
	// tease out syn and non-syn substitutions and sum up count and beta stats  
	void AddSuffStat(const MGOmegaCodonSubMatrix& codonsubmatrix, const PathSuffStat& pathsuffstat)	{
		pathsuffstat.AddOmegaSuffStat(*this,codonsubmatrix);
	}

	// summing over all entries of an array
	void AddSuffStat(const MGOmegaCodonSubMatrix& codonsubmatrix, const PathSuffStatArray& pathsuffstatarray)	{
		for (int i=0; i<pathsuffstatarray.GetSize(); i++)	{
			AddSuffStat(codonsubmatrix,pathsuffstatarray.GetVal(i));
		}
	}

	void Clear()	{
		count = beta = 0;
	}

	void IncrementCount()	{
		count++;
	}

	void AddCount(int in)	{
		count += in;
	}

	void AddBeta(double in)	{
		beta += in;
	}

	void AddSuffStat(int incount, double inbeta)	{
		count += incount;
		beta += inbeta;
	}

	int GetCount() const {
		return count;
	}

	double GetBeta() const {
		return beta;
	}

	double GetLogProb(double rate) const {
		return count*log(rate) - beta*rate;
	}

	double GetMarginalLogProb(double shape, double scale) const {
		return shape*log(scale) - Random::logGamma(shape) - (shape + count)*log(scale + beta) + Random::logGamma(shape + count);
	}

	private:

	int count;
	double beta;
};

void PathSuffStat::AddOmegaSuffStat(OmegaSuffStat& omegasuffstat, const MGOmegaCodonSubMatrix& matrix) const {

    int ncodon = matrix.GetNstate();
    const CodonStateSpace* statespace = matrix.GetCodonStateSpace();

    double beta = 0;
    for (std::map<int,double>::iterator i = waitingtime.begin(); i!= waitingtime.end(); i++)	{
        double totnonsynrate = 0;
        int a = i->first;
        for (int b=0; b<ncodon; b++)	{
            if (b != a)	{
                if (matrix(a,b) != 0)	{
                    if (!statespace->Synonymous(a,b))	{
                        totnonsynrate += matrix(a,b);
                    }
                }
            }
        }
        beta += i->second * totnonsynrate;
    }
    beta /= matrix.GetOmega();

    int count = 0;
    for (std::map<pair<int,int>, int>::iterator i = paircount.begin(); i!= paircount.end(); i++)	{
        if (! statespace->Synonymous(i->first.first,i->first.second))	{
            count += i->second;
        }
    }
    omegasuffstat.AddSuffStat(count,beta);
}

class OmegaSuffStatArray : public SimpleArray<OmegaSuffStat>    {

	public:

	OmegaSuffStatArray(int insize) : SimpleArray<OmegaSuffStat>(insize) {}
	~OmegaSuffStatArray() {}

	void AddSuffStat(const MGOmegaHeterogeneousCodonSubMatrixArray& codonsubmatrixarray, const PathSuffStatArray& pathsuffstatarray)	{
		for (int i=0; i<pathsuffstatarray.GetSize(); i++)	{
            pathsuffstatarray.GetVal(i).AddOmegaSuffStat((*this)[i],codonsubmatrixarray.GetVal(i));
		}
	}

	void AddSuffStat(const MGOmegaHeterogeneousCodonSubMatrixArray& codonsubmatrixarray, const PathSuffStatArray& pathsuffstatarray, const ConstArray<int>& alloc)	{
		for (int i=0; i<pathsuffstatarray.GetSize(); i++)	{
            pathsuffstatarray.GetVal(i).AddOmegaSuffStat((*this)[i],codonsubmatrixarray.GetVal(alloc.GetVal(i)));
		}
	}

	void Clear()	{
		for (int i=0; i<GetSize(); i++)	{
			(*this)[i].Clear();
		}
	}

	double GetLogProb(const Array<double>* ratearray) const{
		double total = 0;
		for (int i=0; i<GetSize(); i++)	{
			total += GetVal(i).GetLogProb(ratearray->GetVal(i));
		}
		return total;
	}

	double GetMarginalLogProb(double shape, double scale)	const {
		double total = 0;
		/*
		for (int i=0; i<GetSize(); i++)	{
			total += GetVal(i).GetMarginalLogProb(shape,scale);
		}
		*/
		// factoring out prior factor
		for (int i=0; i<GetSize(); i++)	{
			int count = GetVal(i).GetCount();
			double beta = GetVal(i).GetBeta();
			total += -(shape+count)*log(scale+beta) + Random::logGamma(shape+count);
		}
		total += GetSize() * (shape*log(scale) - Random::logGamma(shape));
		return total;
	}
};

/*
class BranchOmegaSuffStatArray : public BranchPoissonSuffStatArray	{

	public:
	void AddSuffStat(CodonStateSpace* codonstatespace, BranchPathSuffStatArray* branchpathsuffstatarray)	{
		for (int i=0; i<GetNbranch(); i++)	{
			GetOmegaSuffStat(i).AddSuffStat(codonstatespace,branchpathsuffstatarray->GetVal(i));
		}
	}
};
*/

#endif
