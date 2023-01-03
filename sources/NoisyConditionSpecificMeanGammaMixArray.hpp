#pragma once

#include "BidimArray.hpp"
#include "BidimProduct.hpp"

class NoisyConditionSpecificMeanGammaMixBidimArray : public SimpleBidimArray<double> {

  public:

    NoisyConditionSpecificMeanGammaMixBidimArray(const BidimProduct& inmean, const Selector<double>& intimescale, double ininvshape, double inmean2, double ininvshape2, double inpi, int inrelative) : 
        SimpleBidimArray(inmean.GetNrow(), inmean.GetNcol(), 0), 
        mean(inmean), timescale(intimescale), invshape(ininvshape),
        mean2(inmean2), invshape2(ininvshape2), pi(inpi), relative(inrelative) {
            Sample();
    }

    ~NoisyConditionSpecificMeanGammaMixBidimArray() {}

    void SetParams(double ininvshape,
            double inmean2, double ininvshape2, double inpi) {
        invshape = ininvshape;
        mean2 = inmean2;
        invshape2 = ininvshape2;
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
            double scale1 = shape1 / mean.GetVal(gene, branch);

            double shape2 = 1.0 / invshape2;
            double scale2 = shape2 / mean2;
            if (relative)   {
                scale2 *= timescale.GetVal(branch);
            }

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
        double ret = Random::logGammaDensity(GetVal(gene, branch), shape, scale);
        if (std::isinf(ret) || std::isnan(ret)) {
            cerr << "error in rescaled gamma array get log prob\n";
            cerr << invshape << '\t' << shape << '\t' << scale << '\t' << mean.GetVal(gene, branch) << '\t' << GetVal(gene, branch) << '\n';
            exit(1);
        }
        return ret;
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
            double scale1 = shape1 / mean.GetVal(gene, branch);

            double shape2 = 1.0 / invshape2;
            double scale2 = shape2 / mean2;
            if (relative)   {
                scale2 *= timescale.GetVal(branch);
            }

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
    const Selector<double>& timescale;
    double invshape;
    double mean2;
    double invshape2;
    double pi;
    int relative;
};

