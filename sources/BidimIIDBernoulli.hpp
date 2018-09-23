#ifndef BIDIMIIDBERN_H
#define BIDIMIIDBERN_H

#include "Array.hpp"
#include "BidimArray.hpp"

class BidimIIDBernoulli : public SimpleBidimArray<int> {
  public:
    //! constructor, parameterized by number of rows, of columns, dimension of the
    //! vectors, and the probability vector
    BidimIIDBernoulli(int innrow, int inncol, const Selector<double>& inprob)
        : SimpleBidimArray(innrow, inncol, 0), prob(inprob) {
        Sample();
    }

    //! sample all entries from prior distribution
    void Sample() {
        for (int i = 0; i < GetNrow(); i++) {
            for (int j = 0; j < GetNcol(); j++) {
                Sample(i, j);
            }
        }
    }

    //! set all entries to 0
    void Reset() {
        for (int i = 0; i < GetNrow(); i++) {
            for (int j = 0; j < GetNcol(); j++) {
                (*this)(i, j) =  0;
            }
        }
    }

    //! sample entry i,j
    void Sample(int i, int j) {
        (*this)(i, j) = (Random::Uniform() < prob.GetVal(i));
    }

    //! return total log prob (not yet implemented)
    double GetLogProb() const {
        double total = 0;
        for (int i = 0; i < GetNrow(); i++) {
            for (int j = 0; j < GetNcol(); j++) {
                total += (GetVal(i,j) ? log(prob.GetVal(i)) : log(1.0 - prob.GetVal(i)));
            }
        }
        return total;
    }

  protected:
    const Selector<double> &prob;
};

#endif
