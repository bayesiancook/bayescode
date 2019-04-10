#pragma once

#include "BidimArray.hpp"
#include "global/Random.hpp"

/**
 * \brief A BidimArray of iid vectors (of dimension dim) of Bernoulli variables
 *
 * The Nrow*Ncol bidim array is parameterized by a vector of probability
 * parameters (prob, of dimension Nrow). Then, for i=0..Nrow-1, j=0..Ncol-1,
 * k=0..dim-1:
 *
 * a_ijk ~ Bernoulli(prob[i])
 *
 * This class is used in DiffSelSparseModel.
 * In this context, rows are conditions, columns are sites,
 * and (a_ijk)_k=1..20 is a vector of toggles:
 * a_ijk == 1 means that the fitness of amino-acid k at site j is modified (i.e.
 * resampled from gamma prior) upon undergoing a transition to condition i. This
 * occurs with condition-dependent probability prob[i].
 */

class BidimIIDMultiBernoulli : public SimpleBidimArray<std::vector<int>> {
  public:
    //! constructor, parameterized by number of rows, of columns, dimension of the
    //! vectors, and the probability vector
    BidimIIDMultiBernoulli(int innrow, int inncol, int indim, const std::vector<double> &inprob)
        : SimpleBidimArray(innrow, inncol, std::vector<int>(indim, 1)), dim(indim), prob(inprob) {
        Sample();
    }

    //! return dimension of vectors
    int GetDim() const { return dim; }

    //! sample all entries from prior distribution
    void Sample() {
        for (int i = 0; i < GetNrow(); i++) {
            for (int j = 0; j < GetNcol(); j++) { Sample(i, j); }
        }
    }

    //! set all entries to 0
    void Reset() {
        for (int i = 0; i < GetNrow(); i++) {
            for (int j = 0; j < GetNcol(); j++) {
                std::vector<int> &x = (*this)(i, j);
                for (int k = 0; k < GetDim(); k++) { x[k] = 0; }
            }
        }
    }

    //! sample entry i,j
    void Sample(int i, int j) {
        std::vector<int> &x = (*this)(i, j);
        for (int k = 0; k < GetDim(); k++) { x[k] = (Random::Uniform() < prob[i]); }
    }

    //! get number of shift events (number of 1's out of dim) for row i and column
    //! j
    int GetEventNumber(int i, int j) const {
        int tot = 0;
        const std::vector<int> &x = GetVal(i, j);
        for (int k = 0; k < GetDim(); k++) { tot += x[k]; }
        return tot;
    }

    //! get total number of shift events for row i
    int GetRowEventNumber(int i) const {
        int tot = 0;
        for (int j = 0; j < GetNcol(); j++) { tot += GetEventNumber(i, j); }
        return tot;
    }

    //! return total log prob (not yet implemented)
    double GetLogProb() const { return 0; }

  protected:
    int dim;
    const std::vector<double> &prob;
};