#pragma once

#include "BidimArray.hpp"
#include "BidimProduct.hpp"

template<class Mean>
class GammaBidimArray : public SimpleBidimArray<double> {

  public:

    GammaBidimArray(const Mean& inmean, double ininvshape):
        SimpleBidimArray(inmean.GetNrow(), inmean.GetNcol(), 0), 
        mean(inmean), invshape(ininvshape)  {
            Sample();
    }

    ~GammaBidimArray() {}

    void SetInvShape(double ininvshape) {
        invshape = ininvshape;
    }

    //! return total log prob (over all genes and over all branches)
    double GetLogProb() const {
        double total = 0;
        for (int i=0; i<GetNrow(); i++) {
            for (int j=0; j<GetNcol(); j++)  {
                total += GetLogProb(i,j);
            }
        }
        return total;
    }

    double GetRowLogProb(int i) const    {
        double total = 0;
        for (int j=0; j<GetNcol(); j++) {
            total += GetLogProb(i,j);
        }
        return total;
    }

    double GetColLogProb(int j) const   {
        double total = 0;
        for (int i=0; i<GetNrow(); i++) {
            total += GetLogProb(i,j);
        }
        return total;
    }

    //! return log prob for one entry
    double GetLogProb(int gene, int branch) const {
        double shape = 1.0 / invshape;
        double scale = shape / mean.GetVal(gene, branch);
        return Random::logGammaDensity(GetVal(gene, branch), shape, scale);
    }

    void Sample()   {
        for (int i=0; i<GetNrow(); i++) {
            for (int j=0; j<GetNcol(); j++)  {
                Sample(i,j);
            }
        }
    }
    
    void Sample(int gene, int branch)   {
        double shape = 1.0 / invshape;
        double scale = shape / mean.GetVal(gene, branch);
        (*this)(gene,branch) = Random::Gamma(shape, scale);
    }

  private:
    const Mean& mean;
    // const BidimSelector<double>& mean;
    double invshape;
};

