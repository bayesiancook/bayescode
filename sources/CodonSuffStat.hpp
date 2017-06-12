
#ifndef CODONSUFFSTAT_H
#define CODONSUFFSTAT_H

#include "PathSuffStat.hpp"
#include "CodonSubMatrixArray.hpp"
#include "PoissonSuffStat.hpp"
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

class OmegaSuffStat : public PoissonSuffStat {

	public:

	OmegaSuffStat() {}
	~OmegaSuffStat() {}

	// assumes pathsuffstat is 61x61
	// tease out syn and non-syn substitutions and sum up count and beta stats  
	void AddSuffStat(const MGOmegaCodonSubMatrix& codonsubmatrix, const PathSuffStat& pathsuffstat)	{

        int ncodon = codonsubmatrix.GetNstate();
        const CodonStateSpace* statespace = codonsubmatrix.GetCodonStateSpace();

        const std::map<pair<int,int>,int>& paircount = pathsuffstat.GetPairCountMap();
        const std::map<int,double>& waitingtime = pathsuffstat.GetWaitingTimeMap();

        double tmpbeta = 0;
        for (std::map<int,double>::const_iterator i = waitingtime.begin(); i!= waitingtime.end(); i++)	{
            double totnonsynrate = 0;
            int a = i->first;
            for (int b=0; b<ncodon; b++)	{
                if (b != a)	{
                    if (codonsubmatrix(a,b) != 0)	{
                        if (!statespace->Synonymous(a,b))	{
                            totnonsynrate += codonsubmatrix(a,b);
                        }
                    }
                }
            }
            tmpbeta += i->second * totnonsynrate;
        }
        tmpbeta /= codonsubmatrix.GetOmega();

        int tmpcount = 0;
        for (std::map<pair<int,int>, int>::const_iterator i = paircount.begin(); i!= paircount.end(); i++)	{
            if (! statespace->Synonymous(i->first.first,i->first.second))	{
                tmpcount += i->second;
            }
        }

        PoissonSuffStat::AddSuffStat(tmpcount,tmpbeta);
    }
};

class OmegaSuffStatArray : public SimpleArray<OmegaSuffStat>, public Array<PoissonSuffStat>    {

	public:

	OmegaSuffStatArray(int insize) : SimpleArray<OmegaSuffStat>(insize) {}
	~OmegaSuffStatArray() {}

    int GetSize() const override {return array.size();}
    const OmegaSuffStat& GetVal(int i) const override {return array[i];}
    OmegaSuffStat& operator[](int i) override {return array[i];}

	void Clear()	{
		for (int i=0; i<GetSize(); i++)	{
			(*this)[i].Clear();
		}
	}

	void AddSuffStat(const ConstArray<MGOmegaCodonSubMatrix>& codonsubmatrixarray, const ConstArray<PathSuffStat>& pathsuffstatarray)	{
		for (int i=0; i<GetSize(); i++)	{
            (*this)[i].AddSuffStat(codonsubmatrixarray.GetVal(i),pathsuffstatarray.GetVal(i));
		}
	}

	double GetLogProb(const Array<double>* omegaarray) const{
		double total = 0;
		for (int i=0; i<GetSize(); i++)	{
			total += GetVal(i).GetLogProb(omegaarray->GetVal(i));
		}
		return total;
	}

	double GetMarginalLogProb(double shape, double scale)	const {
		double total = 0;
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

#endif
