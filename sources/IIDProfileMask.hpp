
#include "Array.hpp"
#include "Random.hpp"

#ifndef IIDPROFILEMASK_H
#define IIDPROFILEMASK_H

/**
 * \brief An array of IID 0/1 masks of a fixed dimension
 *
 * This class is used in AAMutSelSparseOmegaModel and DiffSelDoublySparseModel. In those two cases, the masks are over the 20 amino-acids, and the array if of size Nsite (i.e. the array implements site-specific masks over the alignment). The masks are meant to define a low/high fitness distribution over the 20 amino-acids.
 *
 * Each mask is made of dim iid Bernoulli(pi), conditional on at least one entry of the mask being equal to 1.
 */

class IIDProfileMask : public SimpleArray<vector<int> >     {

    public:

    //! constructor, parameterized by array size, mask dimension and Bernoulli probability parameter
    IIDProfileMask(int size, int indim, double pi) : SimpleArray(size,vector<int>(indim,1)), dim(indim), pi(0.1) {}

    //! return dimension of the masks
    int GetDim() const  {
        return dim;
    }

    //! set probability parameter of the Bernoulli to a new value
    void SetPi(double inpi) {
        pi = inpi;
    }

    //! return total log probability over entire array
    double GetLogProb() const   {
        double total = 0;
        for (int i=0; i<GetSize(); i++) {
            total += GetLogProb(i);
        }
        return total;
    }

    //! return log probability for entry i
    double GetLogProb(int i) const  {
        int naa = 0;
        const vector<int>& x = GetVal(i);
        for (int k=0; k<GetDim(); k++)  {
            naa += x[k];
        }
        if (! naa)  {
            cerr << "error in IIDProfileMask: all entries are null\n";
            exit(1);
        }
        // probability is conditional on at least one entry being 1
        return naa*log(pi) + (GetDim()-naa)*log(1.0-pi) - log(1.0 - exp(GetDim()*log(1.0-pi)));
    }

    //! return mean width (i.e. mean number of entries equal to 1) across the array
    double GetMeanWidth() const {
        double mean = 0;
        for (int i=0; i<GetSize(); i++) {
            const vector<int>& x = GetVal(i);
            for (int k=0; k<GetDim(); k++)  {
                mean += x[k];
            }
        }
        mean /= GetSize();
        return mean;
    }

    private:
    int dim;
    double pi;
};

#endif

