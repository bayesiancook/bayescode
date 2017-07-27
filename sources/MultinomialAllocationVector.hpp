#ifndef MULTINOMALLOC_H
#define MULTINOMALLOC_H

#include "Array.hpp"

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
	
	void GibbsResample(const vector<vector<double> >& postprobarray)	{
		for (int i=0; i<GetSize(); i++) {
		    (*this)[i] = Random::DrawFromDiscreteDistribution(postprobarray[i]);
		}

	}

	/*
	void AddSuffStat(vector<int>& occupancy) const {
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
	*/

	private:
	const vector<double>& weight;
};

#endif
