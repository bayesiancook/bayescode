
#ifndef CSMEANGAMMATREE_H
#define CSMEANGAMMATREE_H

#include "MPIBuffer.hpp"
#include "PoissonSuffStat.hpp"
#include "Product.hpp"
#include "Random.hpp"

class ConditionSpecificMeanGammaMixArray : public SimpleArray<double> {
  public:
    ConditionSpecificMeanGammaMixArray(const Selector<double> &inmean, double ininvshape, double ininvshape_ratio, double inpi)
        : SimpleArray<double>(inmean.GetSize()), 
        mean(inmean), invshape(ininvshape),
        invshape_ratio(ininvshape_ratio), pi(inpi) {
        Sample();
    }

    ~ConditionSpecificMeanGammaMixArray() {}

    void SetParams(double ininvshape, double ininvshape_ratio, double inpi) { 
        invshape = ininvshape; 
        invshape_ratio = ininvshape_ratio; 
        pi = inpi;
    }

    //! sample all entries from prior
    void Sample() {
        for (int i = 0; i < GetSize(); i++) {
            if (! invshape) {
                (*this)[i] = mean.GetVal(i);
            }
            else    {
                if (pi && Random::Uniform() < pi) {
                    double shape = 1.0 / invshape / invshape_ratio;
                    double scale = shape / mean.GetVal(i);
                    (*this)[i] = Random::GammaSample(shape, scale);
                }
                else    {
                    double shape = 1.0 / invshape;
                    double scale = shape / mean.GetVal(i);
                    (*this)[i] = Random::GammaSample(shape, scale);
                }
            }
        }
    }

    //! resample entries based on a Array of PoissonSuffStat
    /*
    void GibbsResample(const Selector<PoissonSuffStat> &suffstatarray) {
        if (invshape) {
            for (int i = 0; i < GetSize(); i++) {
                const PoissonSuffStat &suffstat = suffstatarray.GetVal(i);
                if (pi) {
                    double scale = shape / mean.GetVal(i);
                    double shape1 = 1.0 / invshape;
                    double logl1 = suffstat.GetMarginalLogProb(shape1, scale);

                    double shape2 = shape1 / invshape_ratio;
                    double logl2 = suffstat.GetMarginalLogProb(shape2, scale);

                    double max = (logl1 > logl2) ? logl1 : logl2;
                    double l1 = exp(logl1-max);
                    double l2 = exp(logl2-max);
                    double p1 = (1-pi) * l1;
                    double p2 = pi * l2;
                    double tot = p1 + p2;
                    // double logl = log(tot) + max;
                    p1 /= tot;
                    p2 /= tot;

                    if (Random::Uniform() < p2) {
                        (*this)[i] =
                            Random::GammaSample(shape2 + suffstat.GetCount(), scale + suffstat.GetBeta());
                    }
                    else    {
                        (*this)[i] =
                            Random::GammaSample(shape1 + suffstat.GetCount(), scale + suffstat.GetBeta());
                    }
                }
                else    {
                    double shape = 1.0 / invshape;
                    double scale = shape / mean.GetVal(i);
                    (*this)[i] =
                        Random::GammaSample(shape + suffstat.GetCount(), scale + suffstat.GetBeta());
                }
            }
        }
        else    {
            for (int i = 0; i < GetSize(); i++) {
                (*this)[i] = mean.GetVal(i);
            }
        }
    }
    */

    //! return total log prob summed over all entries
    double GetLogProb() const {
        double total = 0;
        if (invshape)   {
            for (int i = 0; i < GetSize(); i++) {
                total += GetLogProb(i);
            }
        }
        return total;
    }

    //! return log prob for one entry
    double GetLogProb(int index) const {
        if (pi) {
            double shape1 = 1.0 / invshape;
            double shape2 = shape1 / invshape_ratio;
            double scale1 = shape1 / mean.GetVal(index);
            double scale2 = shape2 / mean.GetVal(index);
            double logl1 = Random::logGammaDensity(GetVal(index), shape1, scale1);
            double logl2 = Random::logGammaDensity(GetVal(index), shape2, scale2);
            double max = (logl1 > logl2) ? logl1 : logl2;
            double l1 = exp(logl1-max);
            double l2 = exp(logl2-max);
            double logl = log((1-pi)*l1 + pi*l2) + max;
            return logl;
        }
        double shape = 1.0 / invshape;
        double scale = shape / mean.GetVal(index);
        return Random::logGammaDensity(GetVal(index), shape, scale);
    }

    double GetMean() const {
        double m1 = 0;
        for (int i = 0; i < GetSize(); i++) {
            m1 += GetVal(i);
        }
        return m1 / GetSize();
    }

    double GetVar() const {
        double m1 = 0;
        double m2 = 0;
        for (int i = 0; i < GetSize(); i++) {
            m1 += GetVal(i);
            m2 += GetVal(i) * GetVal(i);
        }
        m1 /= GetSize();
        m2 /= GetSize();
        m2 -= m1 * m1;
        return m2;
    }

  protected:
    const Selector<double> &mean;
    double invshape;
    double invshape_ratio;
    double pi;
};

class ConditionSpecificMeanGammaMixBidimArray : public Array<ConditionSpecificMeanGammaMixArray> {
  public:
    //! constructor: parameterized by the number of genes, the tree, the means
    //! over branches and the shape parameter
    ConditionSpecificMeanGammaMixBidimArray(const ProductArray &inmean, double ininvshape, double ininvshape_ratio, double inpi)
        : mean(inmean),
          invshape(ininvshape), invshape_ratio(ininvshape_ratio), pi(inpi),
          array(inmean.GetSize(), (ConditionSpecificMeanGammaMixArray *)0) {
        for (int gene = 0; gene < GetSize(); gene++) {
            array[gene] = new ConditionSpecificMeanGammaMixArray(mean.GetVal(gene), invshape, invshape_ratio, pi);
        }
    }

    ~ConditionSpecificMeanGammaMixBidimArray() {
        for (int gene = 0; gene < GetSize(); gene++) {
            delete array[gene];
        }
    }

    //! set the shape parameter (should be called whenever the shape parameter has
    //! changed during the MCMC)
    void SetParams(double ininvshape, double ininvshape_ratio, double inpi) {
        invshape = ininvshape;
        invshape_ratio = ininvshape_ratio;
        pi = inpi;
        for (int gene = 0; gene < GetSize(); gene++) {
            array[gene]->SetParams(invshape, invshape_ratio, pi);
        }
    }

    //! return total number of entries (number of genes)
    int GetSize() const { return mean.GetSize(); }

    //! return total number of conditions
    int GetNcond() const { return array[0]->GetSize(); }

    //! const access for the given gene
    const ConditionSpecificMeanGammaMixArray &GetVal(int gene) const { return *array[gene]; }

    //! non-const access for the given gene
    ConditionSpecificMeanGammaMixArray &operator[](int gene) { return *array[gene]; }

    //! return mean tree length over genes
    double GetMean() const {
        double tot = 0;
        for (int j = 0; j < GetNcond(); j++) {
            double mean = 0;
            for (int g = 0; g < GetSize(); g++) {
                double tmp = array[g]->GetVal(j);
                mean += tmp;
            }
            mean /= GetSize();
            tot += mean;
        }
        return tot;
    }

    //! return variance of tree length across genes
    double GetVar() const {
        double tot = 0;
        for (int j = 0; j < GetNcond(); j++) {
            double mean = 0;
            double var = 0;
            for (int g = 0; g < GetSize(); g++) {
                double tmp = array[g]->GetVal(j);
                mean += tmp;
                var += tmp * tmp;
            }
            mean /= GetSize();
            var /= GetSize();
            var -= mean * mean;
            tot += var;
        }
        tot /= GetSize();
        return tot;
    }

    //! return total log prob (over all genes and over all branches)
    double GetLogProb() const {
        double total = 0;
        for (int gene = 0; gene < GetSize(); gene++) {
            total += array[gene]->GetLogProb();
        }
        return total;
    }

  private:
    const ProductArray &mean;
    double invshape;
    double invshape_ratio;
    double pi;
    vector<ConditionSpecificMeanGammaMixArray *> array;
};

#endif
