#ifndef MULTINOMALLOC_H
#define MULTINOMALLOC_H

#include "Array.hpp"
#include "OccupancySuffStat.hpp"

class MultinomialAllocationVector : public SimpleArray<int> {

	public:

	MultinomialAllocationVector(int insize, const vector<double>& inweight) : SimpleArray<int>(insize), weight(inweight) {
            SampleAlloc();
	}

	~MultinomialAllocationVector() {}

	void SampleAlloc()  {
		for (int i=0; i<GetSize(); i++) {
			(*this)[i] = Random::DrawFromDiscreteDistribution(weight);
		}
	}
	
    void AddSuffStat(OccupancySuffStat& suffstat) const {
        for (int i=0; i<GetSize(); i++) {
            suffstat.Increment(GetVal(i));
        }
    }

    void GibbsResample(int i, const vector<double>& postprob)   {
        (*this)[i] = Random::DrawFromDiscreteDistribution(postprob);
    }

	void GibbsResample(const vector<vector<double> >& postprobarray)	{
		for (int i=0; i<GetSize(); i++) {
		    (*this)[i] = Random::DrawFromDiscreteDistribution(postprobarray[i]);
		}
	}

    void Permute(const Selector<int>& permut) override {
        if (permut.GetSize() != int(weight.size()))  {
            cerr << "error in MultinomialAllocationVector::Permute: non matching array size\n";
            exit(1);
        }
        vector<int> invpermut(permut.GetSize(),0);
        for (int k=0; k<permut.GetSize(); k++)  {
            invpermut[permut.GetVal(k)] = k;
        }
        for (int i=0; i<GetSize(); i++) {
            (*this)[i] = invpermut[(*this)[i]];
        }
    }

    void SwapComponents(int cat1, int cat2) {
        for (int i=0; i<GetSize(); i++) {
            if ((*this)[i] == cat1)  {
                (*this)[i] = cat2;
            }
            else if ((*this)[i] == cat2)    {
                (*this)[i] = cat1;
            }
        }
    }

	private:
	const vector<double>& weight;
};

#endif
