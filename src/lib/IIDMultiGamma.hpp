
#ifndef IIDMULTIGAMMA_H
#define IIDMULTIGAMMA_H

#include "BidimArray.hpp"
#include "global/Random.hpp"

class IIDMultiGamma : public SimpleArray<std::vector<double>> {
  public:
    //! constructor, parameterized by number of rows, of columns, dimension of the
    //! vectors, shape parameter and center (frequency vector)
    IIDMultiGamma(int insize, int indim, double inshape, const std::vector<double> &incenter)
        : SimpleArray(insize, std::vector<double>(indim, 1.0 / indim)),
          dim(indim),
          shape(inshape),
          center(incenter) {
        Sample();
    }

    //! set shape parameter to new value
    void SetShape(double inshape) { shape = inshape; }

    //! return dimension of vectors
    int GetDim() const { return dim; }

    //! sample all entries from prior distribution
    void Sample() {
        for (int i = 0; i < GetSize(); i++) { Sample(i); }
    }

    void Sample(int i) {
        std::vector<double> &x = (*this)[i];
        for (int k = 0; k < GetDim(); k++) { x[k] = Random::sGamma(shape * center[k]); }
    }

    void PriorResample(const Selector<std::vector<int>> &mask, double min = 0) {
        for (int i = 0; i < GetSize(); i++) {
            std::vector<double> &x = (*this)[i];
            const std::vector<int> &s = mask.GetVal(i);
            for (int k = 0; k < GetDim(); k++) {
                if (!s[k]) {
                    x[k] = Random::sGamma(shape * center[k]);
                    if (x[k] < min) { x[k] = min; }
                }
            }
        }
    }

    void SetUniform() {
        for (int i = 0; i < GetSize(); i++) {
            std::vector<double> &x = (*this)[i];
            for (int k = 0; k < GetDim(); k++) { x[k] = 1.0; }
        }
    }

    //! return total log prob, summed over all entries
    double GetLogProb() const {
        double total = 0;
        for (int i = 0; i < GetSize(); i++) { total += GetLogProb(i); }
        return total;
    }

    //! return log prob for entry i
    double GetLogProb(int i) const {
        const std::vector<double> &x = GetVal(i);
        double total = 0;
        for (int k = 0; k < GetDim(); k++) {
            double alpha = shape * center[k];
            total += -Random::logGamma(alpha) + (alpha - 1) * log(x[k]) - x[k];
        }
        return total;
    }

    double GetLogProb(int i, const std::vector<int> &toggle) const {
        const std::vector<double> &x = GetVal(i);
        double total = 0;
        for (int k = 0; k < GetDim(); k++) {
            if (toggle[k]) {
                double alpha = shape * center[k];
                total += -Random::logGamma(alpha) + (alpha - 1) * log(x[k]) - x[k];
            }
        }
        return total;
    }

  protected:
    int dim;
    double shape;
    const std::vector<double> &center;
};

/**
 * \brief A BidimArray of iid vectors (of dimension dim) of gamma random
 * variables.
 *
 * The Nrow*Ncol array is parameterized by a (positive) shape parameter and a
 * center (a frequency vector of dimension dim). Then, for i=0..Nrow-1,
 * j=0..Ncol-1, k=0..dim-1:
 *
 * x_ijk ~ Gamma(shape,center[k]).
 *
 * This class is used in DiffSelSparseModel.
 * In this context, rows are conditions, columns are sites,
 * and (x_ijk)_k=1..20 is the vector of fitness parameters over the 20
 * amino-acids for site j under condition i.
 *
 * Note that renormalizing those vectors, i.e. defining:
 *
 * y_ijk = x_ijk / Z_ij, where Z_ij = sum_k y_ijk
 *
 * then, the y_ij vectors are iid from a Dirichlet distribution:
 *
 * y_ij ~ Dirichlet(shape*center)
 *
 * In other words, a 'MultiGamma(shape,center)' can be seen as an unnormalized
 * Dirichlet(shape*center), and similarly, a BidimIIDMultiGamma can be seen as a
 * BidimArray of unnormalized iid Dirichlet random variables. In the context of
 * DiffSelSparseModel, working with the unnormalized fitness profiles x_ij,
 * rather than with the normalized versions given by y_ij,
 * turns out to be more practical.
 */

class BidimIIDMultiGamma : public SimpleBidimArray<std::vector<double>> {
  public:
    //! constructor, parameterized by number of rows, of columns, dimension of the
    //! vectors, shape parameter and center (frequency vector)
    BidimIIDMultiGamma(
        int innrow, int inncol, int indim, double inshape, const std::vector<double> &incenter)
        : SimpleBidimArray(innrow, inncol, std::vector<double>(indim, 1.0 / indim)),
          dim(indim),
          shape(inshape),
          center(incenter) {
        Sample();
    }

    //! set shape parameter to new value
    void SetShape(double inshape) { shape = inshape; }

    //! return dimension of vectors
    int GetDim() const { return dim; }

    //! sample all entries from prior distribution
    void Sample() {
        for (int i = 0; i < GetNrow(); i++) {
            for (int j = 0; j < GetNcol(); j++) { Sample(i, j); }
        }
    }

    //! sample entry i,j
    void Sample(int i, int j) {
        std::vector<double> &x = (*this)(i, j);
        for (int k = 0; k < GetDim(); k++) {
            x[k] = Random::sGamma(shape * center[k]);
            // x[k] = Random::Gamma(shape,shape/center[k]);
        }
    }

    //! get mean relative variance for row k (i.e. mean of GetRelVar(k,j) over all
    //! columns j=0..Ncol-1)
    double GetMeanRelVar(int k) const {
        double mean = 0;
        for (int j = 0; j < GetNcol(); j++) { mean += GetRelVar(k, j); }
        mean /= GetNcol();
        return mean;
    }

    //! get relative variance (i.e. heterogeneity) of fitness profile over the 20
    //! amino-acids, for site j in condition k
    double GetRelVar(int k, int j) const {
        double mean = 0;
        double var = 0;
        const std::vector<double> &x = GetVal(k, j);
        for (int l = 0; l < GetDim(); l++) {
            mean += x[l];
            var += x[l] * x[l];
        }
        mean /= GetDim();
        var /= GetDim();
        var -= mean * mean;
        var /= mean * mean;
        return var;
    }

    //! return total log prob, summed over all entries
    double GetLogProb() const {
        double total = 0;
        for (int j = 0; j < GetNcol(); j++) { total += GetColumnLogProb(j); }
        return total;
    }

    //! return total log prob for row i
    double GetRowLogProb(int i) const {
        double total = 0;
        for (int j = 0; j < GetNcol(); j++) { total += GetLogProb(i, j); }
        return total;
    }

    //! return total log prob for column j
    double GetColumnLogProb(int j) const {
        double total = 0;
        for (int i = 0; i < GetNrow(); i++) { total += GetLogProb(i, j); }
        return total;
    }

    //! return total log prob for column j, only for those entries that are
    //! flagged
    double GetColumnLogProb(int j, const std::vector<int> &flag) const {
        double total = 0;
        for (int i = 0; i < GetNrow(); i++) {
            if (flag[i]) { total += GetLogProb(i, j); }
        }
        return total;
    }

    //! return log prob for entry i,j
    double GetLogProb(int i, int j) const {
        const std::vector<double> &x = GetVal(i, j);
        double total = 0;
        for (int k = 0; k < GetDim(); k++) {
            double alpha = shape * center[k];
            total += -Random::logGamma(alpha) + (alpha - 1) * log(x[k]) - x[k];
            // total += shape * log(shape/center[k]) - Random::logGamma(shape) +
            // (shape-1)*log(x[k]) - shape/center[k]*x[k];
        }
        return total;
    }

    //! return logprob for entry i,j, only for those amino-acids for which
    //! toggle[k] != 0
    double GetLogProb(int i, int j, const std::vector<int> &toggle) const {
        const std::vector<double> &x = GetVal(i, j);
        double total = 0;
        for (int k = 0; k < GetDim(); k++) {
            if (toggle[k]) {
                double alpha = shape * center[k];
                total += -Random::logGamma(alpha) + (alpha - 1) * log(x[k]) - x[k];
                // total += shape * log(shape/center[k]) - Random::logGamma(shape) +
                // (shape-1)*log(x[k]) - shape/center[k]*x[k];
            }
        }
        return total;
    }

  protected:
    int dim;
    double shape;
    const std::vector<double> &center;
};

#endif
