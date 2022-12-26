
#ifndef CSMEANGAMMATREE_H
#define CSMEANGAMMATREE_H

#include "BidimArray.hpp"

class BidimProduct {
    public:
        BidimProduct(const Selector<double>& in_rows, const Selector<double>& in_cols) : 
            rows(in_rows), cols(in_cols) {}

        ~BidimProduct() {}

        int GetNrow() const {return rows.GetSize();}
        int GetNcol() const {return cols.GetSize();}

        double GetVal(int i, int j) const {return rows.GetVal(i) * cols.GetVal(j);}

    private:
        const Selector<double>& rows;
        const Selector<double>& cols;
};


class ConditionSpecificMeanGammaMixBidimArray : public SimpleBidimArray<double> {

  public:

    ConditionSpecificMeanGammaMixBidimArray(const BidimProduct& inmean, double ininvshape, double ininvshape_ratio, double inpi) : 
        SimpleBidimArray(inmean.GetNrow(), inmean.GetNcol(), 0), 
        mean(inmean), invshape(ininvshape), invshape_ratio(ininvshape_ratio), pi(inpi) {
            Sample();
    }

    ~ConditionSpecificMeanGammaMixBidimArray() {}

    void SetParams(double ininvshape, double ininvshape_ratio, double inpi) {
        invshape = ininvshape;
        invshape_ratio = ininvshape_ratio;
        pi = inpi;
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
        if (pi) {
            double shape1 = 1.0 / invshape;
            double shape2 = shape1 / invshape_ratio;
            double scale1 = shape1 / mean.GetVal(gene, branch);
            double scale2 = shape2 / mean.GetVal(gene, branch);
            double logl1 = Random::logGammaDensity(GetVal(gene, branch), shape1, scale1);
            double logl2 = Random::logGammaDensity(GetVal(gene, branch), shape2, scale2);
            double max = (logl1 > logl2) ? logl1 : logl2;
            double l1 = exp(logl1-max);
            double l2 = exp(logl2-max);
            double logl = log((1-pi)*l1 + pi*l2) + max;
            return logl;
        }
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
        if (pi) {
            double shape1 = 1.0 / invshape;
            double shape2 = shape1 / invshape_ratio;
            double scale1 = shape1 / mean.GetVal(gene, branch);
            double scale2 = shape2 / mean.GetVal(gene, branch);
            if (Random::Uniform() < pi) {
                (*this)(gene,branch) = Random::Gamma(shape2, scale2);
            }
            else    {
                (*this)(gene,branch) = Random::Gamma(shape1, scale1);
            }
        }
        else    {
            double shape = 1.0 / invshape;
            double scale = shape / mean.GetVal(gene, branch);
            (*this)(gene,branch) = Random::Gamma(shape, scale);
        }
    }

  private:
    const BidimProduct& mean;
    double invshape;
    double invshape_ratio;
    double pi;
};

#endif
