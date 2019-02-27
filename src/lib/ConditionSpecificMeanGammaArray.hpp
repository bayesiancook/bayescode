#pragma once

#include "PoissonSuffStat.hpp"
#include "Product.hpp"
#include "global/Random.hpp"

class ConditionSpecificMeanGammaArray : public SimpleArray<double> {
  public:
    ConditionSpecificMeanGammaArray(const Selector<double> &inblmean, double ininvshape)
        : SimpleArray<double>(inblmean.GetSize()), blmean(inblmean), invshape(ininvshape) {
        Sample();
    }

    ~ConditionSpecificMeanGammaArray() {}

    double GetInvShape() const { return invshape; }

    void SetInvShape(double ininvshape) { invshape = ininvshape; }

    //! sample all entries from prior
    void Sample() {
        for (int i = 0; i < GetSize(); i++) {
            double shape = 1.0 / invshape;
            double scale = shape / blmean.GetVal(i);
            (*this)[i] = Random::GammaSample(shape, scale);
        }
    }

    //! resample entries based on a Array of PoissonSuffStat
    void GibbsResample(const Selector<PoissonSuffStat> &suffstatarray) {
        for (int i = 0; i < GetSize(); i++) {
            const PoissonSuffStat &suffstat = suffstatarray.GetVal(i);
            double shape = 1.0 / invshape;
            double scale = shape / blmean.GetVal(i);
            (*this)[i] =
                Random::GammaSample(shape + suffstat.GetCount(), scale + suffstat.GetBeta());
        }
    }

    //! return total log prob summed over all entries
    double GetLogProb() {
        double total = 0;
        for (int i = 0; i < GetSize(); i++) { total += GetLogProb(i); }
        return total;
    }

    //! return log prob for one entry
    double GetLogProb(int index) const {
        double shape = 1.0 / invshape;
        double scale = shape / blmean.GetVal(index);
        return Random::logGammaDensity(GetVal(index), shape, scale);
    }

    double GetMean() const {
        double m1 = 0;
        for (int i = 0; i < GetSize(); i++) { m1 += GetVal(i); }
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
    const Selector<double> &blmean;
    double invshape;
};

class ConditionSpecificMeanGammaBidimArray : public Array<ConditionSpecificMeanGammaArray> {
  public:
    //! constructor: parameterized by the number of genes, the tree, the means
    //! over branches and the shape parameter
    ConditionSpecificMeanGammaBidimArray(const ProductArray &inmean, double ininvshape)
        : mean(inmean),
          invshape(ininvshape),
          array(inmean.GetSize(), (ConditionSpecificMeanGammaArray *)0) {
        for (int gene = 0; gene < GetSize(); gene++) {
            array[gene] = new ConditionSpecificMeanGammaArray(mean.GetVal(gene), invshape);
        }
    }

    ~ConditionSpecificMeanGammaBidimArray() {
        for (int gene = 0; gene < GetSize(); gene++) { delete array[gene]; }
    }

    //! set the shape parameter (should be called whenever the shape parameter has
    //! changed during the MCMC)
    void SetInvShape(double ininvshape) {
        invshape = ininvshape;
        for (int gene = 0; gene < GetSize(); gene++) { array[gene]->SetInvShape(invshape); }
    }

    //! return total number of entries (number of genes)
    int GetSize() const { return mean.GetSize(); }

    //! return total number of conditions
    int GetNcond() const { return array[0]->GetSize(); }

    //! const access for the given gene
    const ConditionSpecificMeanGammaArray &GetVal(int gene) const { return *array[gene]; }

    //! non-const access for the given gene
    ConditionSpecificMeanGammaArray &operator[](int gene) { return *array[gene]; }

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
        for (int gene = 0; gene < GetSize(); gene++) { total += array[gene]->GetLogProb(); }
        return total;
    }

  private:
    const ProductArray &mean;
    double invshape;
    vector<ConditionSpecificMeanGammaArray *> array;
};