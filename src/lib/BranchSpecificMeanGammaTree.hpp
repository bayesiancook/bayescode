#pragma once

#include "BranchProduct.hpp"
#include "PoissonSuffStat.hpp"
#include "global/Random.hpp"

class BranchSpecificMeanGammaTree : public SimpleBranchArray<double> {
  public:
    BranchSpecificMeanGammaTree(const BranchSelector<double> &inblmean, double ininvshape)
        : SimpleBranchArray<double>(inblmean.GetTree()), blmean(inblmean), invshape(ininvshape) {
        Sample();
    }

    ~BranchSpecificMeanGammaTree() {}

    double GetInvShape() const { return invshape; }

    void SetInvShape(double ininvshape) { invshape = ininvshape; }

    //! sample all entries from prior
    void Sample() {
        for (int i = 0; i < GetNbranch(); i++) {
            double shape = 1.0 / invshape;
            double scale = shape / blmean.GetVal(i);
            (*this)[i] = Random::GammaSample(shape, scale);
        }
    }

    //! resample entries based on a BranchArray of PoissonSuffStat
    void GibbsResample(const BranchSelector<PoissonSuffStat> &suffstatarray) {
        for (int i = 0; i < GetNbranch(); i++) {
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
        for (int i = 0; i < GetNbranch(); i++) { total += GetLogProb(i); }
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
        for (int i = 0; i < GetNbranch(); i++) { m1 += GetVal(i); }
        return m1 / GetNbranch();
    }

    double GetVar() const {
        double m1 = 0;
        double m2 = 0;
        for (int i = 0; i < GetNbranch(); i++) {
            m1 += GetVal(i);
            m2 += GetVal(i) * GetVal(i);
        }
        m1 /= GetNbranch();
        m2 /= GetNbranch();
        m2 -= m1 * m1;
        return m2;
    }

  protected:
    const BranchSelector<double> &blmean;
    double invshape;
};

class BranchSpecificMeanGammaTreeArray : public Array<BranchSpecificMeanGammaTree> {
  public:
    //! constructor: parameterized by the number of genes, the tree, the means
    //! over branches and the shape parameter
    BranchSpecificMeanGammaTreeArray(const BranchProductArray &inmean, double ininvshape)
        : mean(inmean),
          invshape(ininvshape),
          array(inmean.GetSize(), (BranchSpecificMeanGammaTree *)0) {
        for (int gene = 0; gene < GetSize(); gene++) {
            array[gene] = new BranchSpecificMeanGammaTree(mean.GetVal(gene), invshape);
        }
    }

    ~BranchSpecificMeanGammaTreeArray() {
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

    //! return total number of branches of the underlying tree
    int GetNbranch() const { return array[0]->GetNbranch(); }

    //! const access for the given gene
    const BranchSpecificMeanGammaTree &GetVal(int gene) const { return *array[gene]; }

    //! non-const access for the given gene
    BranchSpecificMeanGammaTree &operator[](int gene) { return *array[gene]; }

    //! return mean tree length over genes
    double GetMean() const {
        double tot = 0;
        for (int j = 0; j < GetNbranch(); j++) {
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
        for (int j = 0; j < GetNbranch(); j++) {
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
        tot /= GetNbranch();
        return tot;
    }

    //! return total log prob (over all genes and over all branches)
    double GetLogProb() const {
        double total = 0;
        for (int gene = 0; gene < GetSize(); gene++) { total += array[gene]->GetLogProb(); }
        return total;
    }

  private:
    const BranchProductArray &mean;
    double invshape;
    vector<BranchSpecificMeanGammaTree *> array;
};