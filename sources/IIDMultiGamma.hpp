
#ifndef IIDMULTIGAMMA_H
#define IIDMULTIGAMMA_H

#include "BidimArray.hpp"


/**
 * \brief A BidimArray of Nrow * Ncol iid vectors of gamma random variables (each vector being itself of dimension dim)
 *
 * This class is used in DiffSelSparseModel.
 * In this context, rows are conditions, columns are sites, 
 * and (x_ijk)_k=1..20 is the vector of fitness parameters over the 20 amino-acids for site j under condition i.
 * Thus, for i=0..Nrow-1, j=0..Ncol-1, k=0..dim-1:
 *
 * x_ijk ~ Gamma(shape,center[k]).
 */

class BidimIIDMultiGamma : public SimpleBidimArray<vector<double> >    {

    public:

    //! constructor, parameterized by number of rows, of columns, dimension of the vectors, shape parameter and vector of means (center)
    BidimIIDMultiGamma(int innrow, int inncol, int indim, double inshape, const vector<double>& incenter) :
        SimpleBidimArray(innrow,inncol,vector<double>(indim,1.0/indim)), dim(indim), shape(inshape), center(incenter) {
        Sample();
    }

    //! set shape parameter to new value
    void SetShape(double inshape)   {
        shape = inshape;
    }

    //! return dimension of vectors
    int GetDim() const {return dim;}

    //! sample all entries from prior distribution
    void Sample()   {
        for (int i=0; i<GetNrow(); i++)  {
            for (int j=0; j<GetNcol(); j++)   {
                Sample(i,j);
            }
        }
    }

    //! sample entry i,j
    void Sample(int i, int j)   {
        vector<double>& x = (*this)(i,j);
        for (int k=0; k<GetDim(); k++) {
            x[k] = Random::sGamma(shape*center[k]);
            // x[k] = Random::Gamma(shape,shape/center[k]);
        }
    }

    //! get mean relative variance over row (condition) k (mean over all sites)
    double GetMeanRelVar(int k) const  {

        double mean = 0;
        for (int j=0; j<GetNcol(); j++) {
            mean += GetRelVar(k,j);
        }
        mean /= GetNcol();
        return mean;
    }

    //! get relative variance (i.e. heterogeneity) of fitness profile over the 20 amino-acids, for site j in condition k
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

    //! return total log prob, summed over all entries
    double GetLogProb() const {
        double total = 0;
        for (int j=0; j<GetNcol(); j++)   {
            total += GetColumnLogProb(j);
        }
        return total;
    }

    //! return total log prob for row i
    double GetRowLogProb(int i) const {
        double total = 0;
        for (int j=0; j<GetNcol(); j++)   {
            total += GetLogProb(i,j);
        }
        return total;
    }

    //! return total log prob for column j 
    double GetColumnLogProb(int j) const  {
        double total = 0;
        for (int i=0; i<GetNrow(); i++)  {
            total += GetLogProb(i,j);
        }
        return total;
    }

    //! return total log prob for column j, only for those entries that are flagged
    double GetColumnLogProb(int j, const vector<int>& flag) const   {
        double total = 0;
        for (int i=0; i<GetNrow(); i++)  {
            if (flag[i])    {
                total += GetLogProb(i,j);
            }
        }
        return total;
    }

    //! return log prob for entry i,j
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

    //! return logprob for entry i,j, only for those amino-acids that for which toggle[k] != 0
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

