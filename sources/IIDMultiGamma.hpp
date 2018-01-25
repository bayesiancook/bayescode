
#ifndef IIDMULTIGAMMA_H
#define IIDMULTIGAMMA_H

#include "BidimArray.hpp"


class BidimIIDMultiGamma : public SimpleBidimArray<vector<double> >    {

    public:

    // columns: sites
    // rows: conditions
    BidimIIDMultiGamma(int innrow, int inncol, int indim, double inshape, const vector<double>& incenter) :
        SimpleBidimArray(innrow,inncol,vector<double>(indim,1.0/indim)), dim(indim), shape(inshape), center(incenter) {
        Sample();
    }

    void SetShape(double inshape)   {
        shape = inshape;
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
        vector<double>& x = (*this)(i,j);
        for (int k=0; k<GetDim(); k++) {
            x[k] = Random::sGamma(shape*center[k]);
            // x[k] = Random::Gamma(shape,shape/center[k]);
        }
    }

    double GetMeanRelVar(int k) const  {

        double mean = 0;
        for (int j=0; j<GetNcol(); j++) {
            mean += GetRelVar(k,j);
        }
        mean /= GetNcol();
        return mean;
    }

    double GetRelVar(int k, int j) const   {

        double mean = 0;
        double var = 0;
        const vector<double>& x = GetVal(k,j);
        for (int l=0; l<GetDim(); l++)  {
            mean += x[l];
            var += x[l]*x[l];
        }
        mean /= GetDim();
        var /= GetDim();
        var -= mean*mean;
        var /= mean*mean;
        return var;
    }

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
            double alpha = shape*center[k];
            total += - Random::logGamma(alpha) + (alpha-1)*log(x[k]) - x[k];
            // total += shape * log(shape/center[k]) - Random::logGamma(shape) + (shape-1)*log(x[k]) - shape/center[k]*x[k];
        }
        return total;
    }

    double GetLogProb(int i, int j, const vector<int>& toggle) const {
        const vector<double>& x = GetVal(i,j);
        double total = 0;
        for (int k=0; k<GetDim(); k++) {
            if (toggle[k])  {
                double alpha = shape*center[k];
                total += - Random::logGamma(alpha) + (alpha-1)*log(x[k]) - x[k];
                // total += shape * log(shape/center[k]) - Random::logGamma(shape) + (shape-1)*log(x[k]) - shape/center[k]*x[k];
            }
        }
        return total;
    }

    protected:

    int dim;
    double shape;
    const vector<double>& center;
};

#endif

