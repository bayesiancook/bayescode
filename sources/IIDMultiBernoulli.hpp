
#ifndef IIDMULTIBERN_H
#define IIDMULTIBERN_H

#include "BidimArray.hpp"

class BidimIIDMultiBernoulli : public SimpleBidimArray<vector<int> >    {

    public:

    // columns: sites
    // rows: conditions
    BidimIIDMultiBernoulli(int inncol, int innrow, int indim, const vector<double>& inprob) :
    // BidimIIDMultiBernoulli(int inncol, int innrow, int indim, const Selector<double>& inprob) :
        SimpleBidimArray(inncol,innrow,vector<int>(indim,1)), dim(indim), prob(inprob)    {
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
            x[k] = (Random::Uniform() < prob[k]);
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
    /*
    double GetLogProb() const {
        double total = 0;
        for (int j=0; j<GetNcol(); j++)   {
            total += GetColumnLogProb(j);
        }
        return total;
    }

    double GetRowLogProb(int i) const {
        double total = 0;
        for (int j=0; j<GetNcol(); j++)   {
            total += GetLogProb(i,j);
        }
        return total;
    }

    double GetColumnLogProb(int j) const  {
        double total = 0;
        for (int i=0; i<GetNrow(); i++)  {
            total += GetLogProb(i,j);
        }
        return total;
    }

    double GetColumnLogProb(int j, const vector<int>& flag) const   {
        double total = 0;
        for (int i=0; i<GetNrow(); i++)  {
            if (flag[i])    {
                total += GetLogProb(i,j);
            }
        }
        return total;
    }

    double GetLogProb(int i, int j) const {
        const vector<double>& x = GetVal(i,j);
        double total = 0;
        for (int k=0; k<GetDim(); k++) {
            total += shape * log(shape/center[k]) - Random::logGamma(shape) + (shape-1)*log(x[k]) - shape/center[k]*x[k];
        }
        return total;
    }
    */

    protected:

    int dim;
    const vector<double>& prob;
};

#endif


