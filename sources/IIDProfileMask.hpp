
#include "Array.hpp"
#include "Random.hpp"

#ifndef IIDPROFILEMASK_H
#define IIDPROFILEMASK_H

class IIDProfileMask : public SimpleArray<vector<int> >     {

    public:

    IIDProfileMask(int size, int indim, double pi) : SimpleArray(size,vector<int>(indim,1)), dim(indim), pi(0.1) {}

    int GetDim() const  {
        return dim;
    }

    void SetPi(double inpi) {
        pi = inpi;
    }

    double GetLogProb() const   {
        double total = 0;
        for (int i=0; i<GetSize(); i++) {
            total += GetLogProb(i);
        }
        return total;
    }

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

