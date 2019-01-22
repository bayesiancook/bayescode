
#ifndef IIDMVNORMAL_H
#define IIDMVNORMAL_H

#include "Array.hpp"
#include "BidimArray.hpp"

/**
 * \brief A BidimArray of iid vectors (of dimension dim) of normal random
 * variates.
 *
 * The Nrow*Ncol array is parameterized by a vector of variance parameters (var,
 * of dimension Nrow). Then, for i=0..Nrow-1, j=0..Ncol-1, k=0..dim-1:
 *
 * x_ijk ~ Normal(0,var[i])
 *
 * This class is used in DiffSelModel.
 * In this context, rows are conditions, columns are sites,
 * and (x_ijk)_k=1..20 is the vector of differential selection effects over the
 * 20 amino-acids, for site j and under condition i (compared to reference
 * condition).
 */

class BidimIIDMVNormal : public SimpleBidimArray<vector<double>> {
  public:
    //! \brief constructor, parameterized by number of columns, dimension of
    //! vectors and vector of variance parameters
    //!
    //! number of rows is defined implicitly: equal to dimension of the vector of
    //! variance parameters.
    BidimIIDMVNormal(int inncol, int indim, const Selector<double> &invar)
        : SimpleBidimArray(invar.GetSize(), inncol, vector<double>(indim, 0)),
          dim(indim),
          var(invar) {
        Sample();
    }

    //! return dimension of vectors
    int GetDim() const { return dim; }

    //! sample all entries from prior distribution
    void Sample() {
        for (int i = 0; i < GetNrow(); i++) {
            for (int j = 0; j < GetNcol(); j++) {
                Sample(i, j);
            }
        }
    }

    //! set all entries to 0
    void SetToZero() {
        for (int i = 0; i < GetNrow(); i++) {
            for (int j = 0; j < GetNcol(); j++) {
                vector<double> &x = (*this)(i, j);
                for (int k = 0; k < GetDim(); k++) {
                    x[k] = 0;
                }
            }
        }
    }

    //! sample entry i,j
    void Sample(int i, int j) {
        vector<double> &x = (*this)(i, j);
        double v = var.GetVal(i);
        for (int k = 0; k < GetDim(); k++) {
            x[k] = sqrt(v) * Random::sNormal();
        }
    }

    //! return mean (over columns, or sites) of variance (across rows, or
    //! conditions) for component (or amino-acid) k.
    double GetMeanVar(int k) const {
        double mean = 0;
        for (int j = 0; j < GetNcol(); j++) {
            mean += GetVar(k, j);
        }
        mean /= GetNcol();
        return mean;
    }

    //! return variance (across rows, or conditions) for component (or amino-acid)
    //! k.
    double GetVar(int k, int j) const {
        double mean = 0;
        double var = 0;
        const vector<double> &x = GetVal(k, j);
        for (int l = 0; l < GetDim(); l++) {
            mean += x[l];
            var += x[l] * x[l];
        }
        mean /= GetDim();
        var /= GetDim();
        var -= mean * mean;
        return var;
    }

    //! return total log prob over array
    double GetLogProb() const {
        double total = 0;
        for (int j = 0; j < GetNcol(); j++) {
            total += GetColumnLogProb(j);
        }
        return total;
    }

    //! return total log prob over row i
    double GetRowLogProb(int i) const {
        double total = 0;
        for (int j = 0; j < GetNcol(); j++) {
            total += GetLogProb(i, j);
        }
        return total;
    }

    //! return total log prob over column j
    double GetColumnLogProb(int j) const {
        double total = 0;
        for (int i = 0; i < GetNrow(); i++) {
            total += GetLogProb(i, j);
        }
        return total;
    }

    //! return log prob over column j, counting only those entries that are
    //! flagged
    double GetColumnLogProb(int j, const vector<int> &flag) const {
        double total = 0;
        for (int i = 0; i < GetNrow(); i++) {
            if (flag[i]) {
                total += GetLogProb(i, j);
            }
        }
        return total;
    }

    //! return log prob contribution of entry i,j
    double GetLogProb(int i, int j) const {
        double total = 0;
        const vector<double> &x = GetVal(i, j);
        double v = var.GetVal(i);
        for (int k = 0; k < GetDim(); k++) {
            total += x[k] * x[k];
        }
        return -0.5 * (total / v + GetDim() * log(2 * Pi * v));
    }

  protected:
    int dim;
    const Selector<double> &var;
};

/*
class IIDMVNormal : public SimpleArray<vector<double >    {

    public:

    IIDMVNormal(int insize, int indim, double invar) :
SimpleArray<vector<double> >(insize,vector<double>(indim,0)), dim(indim),
var(invar)  { Sample();
    }

    ~IIDMVNormal() {}

    void Sample()   {
        for (int i=0; i<GetSize(); i++) {
            for (int j=0; j<GetDim(); j++) {
                (*this)[i][j] = sqrt(var)*Random::sNormal();
            }
        }
    }

    double GetLogProb() {
        double total = 0;
        for (int i=0; i<GetSize(); i++) {
            total += GetLogProb(i);
        }
        return total;
    }

    double GetLogProb(int i)    {
        double total = 0;
        for (int j=0; j<GetDim(); j++)  {
            double tmp = GetVal(i)[j];
            total += tmp*tmp;
        }
        double ret = -0.5*(total/var + GetDim()*log(var));
    }

    int GetDim()    {
        return dim;
    }

    void SetVar(double invar)   {
        var = invar;
    }

    protected:

    int dim;

    double var;

};
*/

#endif
