
#ifndef IIDMULTIBERN_H
#define IIDMULTIBERN_H

#include "Array.hpp"
#include "BidimArray.hpp"

class IIDMultiBernoulli : public SimpleArray<vector<int> >    {

    public:

    IIDMultiBernoulli(int insize, int indim, double inprob) :
        SimpleArray(insize,vector<int>(indim,1)), dim(indim), prob(inprob)    {
        Sample();
    }

    int GetDim() const {return dim;}

    void SetProb(double inprob) {
        prob = inprob;
    }

    void Sample()   {
        for (int i=0; i<GetSize(); i++) {
            Sample(i);
        }
    }

    void Sample(int i)  {
        vector<int>& x = (*this)[i];
        for (int k=0; k<GetDim(); k++) {
            x[k] = (Random::Uniform() < prob);
        }
    }

    int GetEventNumber(int i) const {
        int tot = 0;
        const vector<int>& x = GetVal(i);
        for (int k=0; k<GetDim(); k++)  {
            tot += x[k];
        }
        return tot;
    }

    double GetLogProb() const {
        return 0;
    }

    protected:

    int dim;
    double prob;
};

class BidimIIDMultiBernoulli : public SimpleBidimArray<vector<int> >    {

    public:

    // columns: sites
    // rows: conditions
    BidimIIDMultiBernoulli(int innrow, int inncol, int indim, const vector<double>& inprob) :
        SimpleBidimArray(innrow,inncol,vector<int>(indim,1)), dim(indim), prob(inprob)    {
        Sample();
    }

    int GetDim() const {return dim;}

    void Sample()   {
        for (int i=0; i<GetNrow(); i++)  {
            for (int j=0; j<GetNcol(); j++)   {
                Sample(i,j);
            }
        }
    }

    void Sample(int i, int j)   {
        vector<int>& x = (*this)(i,j);
        for (int k=0; k<GetDim(); k++) {
            x[k] = (Random::Uniform() < prob[i]);
        }
    }

    int GetEventNumber(int i, int j) const {
        int tot = 0;
        const vector<int>& x = GetVal(i,j);
        for (int k=0; k<GetDim(); k++)  {
            tot += x[k];
        }
        return tot;
    }

    int GetRowEventNumber(int j) const {
        int tot = 0;
        for (int i=0; i<GetNcol(); i++) {
            tot += GetEventNumber(i,j);
        }
        return tot;
    }

    double GetLogProb() const {
        return 0;
    }

    protected:

    int dim;
    const vector<double>& prob;
};

#endif


