#ifndef MULTINOMALLOC_H
#define MULTINOMALLOC_H

#include "Array.hpp"

class MultinomialAllocationVector : public SimpleArray<int> {

	public:

	MultinomialAllocationVector(int insize, const vector<double>& inweight) : SimpleArray<int>(insize), weight(inweight), occupancy(inweight.size()) {
            SampleAlloc();
	}

	~MultinomialAllocationVector() {}

	void SampleAlloc()  {
		for (int i=0; i<GetSize(); i++) {
			(*this)[i] = Random::DrawFromDiscreteDistribution(weight);
		}
        UpdateOccupancies();
	}
	
    void GibbsResample(int i, const vector<double>& postprob)   {
        (*this)[i] = Random::DrawFromDiscreteDistribution(postprob);
    }

	void GibbsResample(const vector<vector<double> >& postprobarray)	{
		for (int i=0; i<GetSize(); i++) {
		    (*this)[i] = Random::DrawFromDiscreteDistribution(postprobarray[i]);
		}
        UpdateOccupancies();
	}

    const vector<int>& GetOccupancies() const   {
        return occupancy;
    }

	void UpdateOccupancies() const {
		if (occupancy.size() != weight.size())  {
			cerr << "error: non matching size\n";
			exit(1);
		}
		for (unsigned int k=0 ; k<occupancy.size(); k++)  {
			occupancy[k] = 0;
		}
		for (int i=0; i<GetSize(); i++) {
			occupancy[GetVal(i)]++;
		}
	}

	private:
	const vector<double>& weight;
    mutable vector<int> occupancy;
};

#endif
