#ifndef BIDIMIIDGAMMA_H
#define BIDIMIIDGAMMA_H

#include "BidimArray.hpp"

class BidimIIDGamma : public SimpleBidimArray<double> {
  public:
    //! constructor, parameterized by number of rows, of columns, dimension of the
    //! vectors, shape parameter and center (frequency vector)
    BidimIIDGamma(int innrow, int inncol, const Selector<double>& inshape, const Selector<double>& inscale)   {
        : SimpleBidimArray<double>(innrow, inncol, 1.0), 
          shape(inshape),
          scale(inscale) {
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

    //! sample entry i,j
    void Sample(int i, int j) {
        (*this)(i, j) = Random::GammaSample(shape.GetVal(i),scale.GetVal(i));
    }

    //! return total log prob, summed over all entries
    double GetLogProb() const {
        double total = 0;
        for (int j = 0; j < GetNcol(); j++) {
            total += GetColumnLogProb(j);
        }
        return total;
    }

    //! return total log prob for row i
    double GetRowLogProb(int i) const {
        double total = 0;
        for (int j = 0; j < GetNcol(); j++) {
            total += GetLogProb(i, j);
        }
        return total;
    }

    //! return total log prob for column j
    double GetColumnLogProb(int j) const {
        double total = 0;
        for (int i = 0; i < GetNrow(); i++) {
            total += GetLogProb(i, j);
        }
        return total;
    }

    //! return log prob for entry i,j
    double GetLogProb(int i, int j) const {
        return Random::logGammaDensity(GetVal(i,j),shape.GetVal(i),scale.GetVal(i));
    }

  protected:
    const Selector<double>& shape;
    double Selector<double>& scale;
};

#endif
