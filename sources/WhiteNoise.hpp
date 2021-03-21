#ifndef WHITENOISE_H
#define WHITENOISE_H

#include "BranchArray.hpp"
#include "PoissonSuffStat.hpp"

/**
 * \brief A tree-structured branch-wise array of Gamma variables, with
 * branch-specific means but same shape parameter
 *
 * One should be careful about the fact that the shape parameter is given by
 * copy (not by ref) to the array. Thus, each time the shape parameter is
 * modified during the MCMC, the new value should be given to the array (using
 * the SetShape method).
 */
class GammaWhiteNoise : public SimpleBranchArray<double> {
  public:

    // mode 0   : ugam
    // mode 1/2 : wn of mean 1 and variance per unit of time = 1/shape
    // mode 1   : realized branch variable: integral of wn over time interval = blmean
    // mode 2   : realized branch variable: mean     of wn over time interval = blmean

    GammaWhiteNoise(const Tree &intree, const BranchSelector<double> &inblmean, double inshape, int inmode = 0)
        : SimpleBranchArray<double>(intree), blmean(inblmean), shape(inshape) {
        mode = inmode;
        Sample();
    }

    ~GammaWhiteNoise() {}

    void SetMode(int inmode)    {
        mode = inmode;
    }

    double GetShape() const { return shape; }

    void SetShape(double inshape) { shape = inshape; }

    double GetAlpha(int i) const {
        if (mode == 2)  {
            return shape*blmean.GetVal(i);
        }
        if (mode == 1)   {
            return shape*blmean.GetVal(i);
        }
        return shape;
    }

    double GetBeta(int i) const {
        if (mode == 2)  {
            return shape*blmean.GetVal(i);
        }
        if (mode == 1)   {
            return shape;
        }
        return shape / blmean.GetVal(i);
    }
        
    //! sample all entries from prior
    void Sample() {
        for (int i = 0; i < GetNbranch(); i++) {
            (*this)[i] = Random::GammaSample(GetAlpha(i), GetBeta(i));
        }
    }

    //! resample entries based on a BranchArray of PoissonSuffStat
    void GibbsResample(const PoissonSuffStatBranchArray &suffstatarray) {
        for (int i = 0; i < GetNbranch(); i++) {
            const PoissonSuffStat &suffstat = suffstatarray.GetVal(i);
            (*this)[i] = Random::GammaSample(GetAlpha(i) + suffstat.GetCount(), GetBeta(i) + suffstat.GetBeta());
            if ((*this)[i] == 0)    {
                cerr << "gibbs null bl : " << GetAlpha(i) << '\t' << GetBeta(i) << '\t' << shape << '\t' << blmean.GetVal(i) << '\n';
                (*this)[i] = 0.001;
            }
        }
    }

    //! resample entries based on a BranchArray of PoissonSuffStat
    void ResampleEmptyBranches(const PoissonSuffStatBranchArray &suffstatarray) {
        for (int i = 0; i < GetNbranch(); i++) {
            if (suffstatarray.GetVal(i).GetBeta() == 0)  {
                (*this)[i] = Random::GammaSample(GetAlpha(i), GetBeta(i));
                if ((*this)[i] == 0)    {
                    cerr << "empty null bl : " << GetAlpha(i) << '\t' << GetBeta(i) << '\t' << shape << '\t' << blmean.GetVal(i) << '\n';
                    (*this)[i] = 0.001;
                }
            }
        }
    }

    //! return total log prob summed over all entries
    double GetLogProb() {
        double total = 0;
        for (int i = 0; i < GetNbranch(); i++) {
            total += GetLogProb(i);
        }
        return total;
    }

    //! return log prob for one entry
    double GetLogProb(int index) const {
        return Random::logGammaDensity(GetVal(index), GetAlpha(index), GetBeta(index));
    }

    double GetTotalLength() const {
        double m1 = 0;
        for (int i = 0; i < GetNbranch(); i++) {
            m1 += GetVal(i);
        }
        return m1;
    }

  protected:
    const BranchSelector<double> &blmean;
    double shape;
    int mode;
};

/**
 * \brief An array of GammaWhiteNoise
 *
 * This is an Array of BranchArray's. It is used for implementing modulations of
 * branch-lengths across genes, while implementing some basic shrinkage.
 * Specifically, for a given branch j,
 * the branch-lengths across genes are gamma of mean blmean[j] and shape
 * parameter uniform across all branches.
 */

class GammaWhiteNoiseArray : public Array<GammaWhiteNoise> {
  public:
    //! constructor: parameterized by the number of genes, the tree, the means
    //! over branches and the shape parameter
    GammaWhiteNoiseArray(int inNgene, const Tree &intree, const BranchSelector<double> &inblmean,
                         double inshape)
        : Ngene(inNgene),
          tree(intree),
          blmean(inblmean),
          shape(inshape),
          blarray(Ngene, (GammaWhiteNoise *)0) {
        for (int gene = 0; gene < Ngene; gene++) {
            blarray[gene] = new GammaWhiteNoise(tree, blmean, shape);
        }
    }

    ~GammaWhiteNoiseArray() {
        for (int gene = 0; gene < Ngene; gene++) {
            delete blarray[gene];
        }
    }

    //! set the shape parameter (should be called whenever the shape parameter has
    //! changed during the MCMC)
    void SetShape(double inshape) {
        shape = inshape;
        for (int gene = 0; gene < Ngene; gene++) {
            blarray[gene]->SetShape(shape);
        }
    }

    //! return total number of entries (number of genes)
    int GetSize() const { return Ngene; }

    //! return total number of genes
    int GetNgene() const { return Ngene; }

    //! return total number of branches of the underlying tree
    int GetNbranch() const { return blarray[0]->GetNbranch(); }

    //! const access to the GammaWhiteNoise BranchArray for the given gene
    const GammaWhiteNoise &GetVal(int gene) const { return *blarray[gene]; }

    //! non-const access to the GammaWhiteNoise BranchArray for the given gene
    GammaWhiteNoise &operator[](int gene) { return *blarray[gene]; }

    //! return mean tree length over genes
    double GetMeanLength() const {
        double tot = 0;
        for (int j = 0; j < GetNbranch(); j++) {
            double mean = 0;
            for (int g = 0; g < GetNgene(); g++) {
                double tmp = blarray[g]->GetVal(j);
                mean += tmp;
            }
            mean /= GetNgene();
            tot += mean;
        }
        return tot;
    }

    //! return variance of tree length across genes
    double GetVarLength() const {
        double tot = 0;
        for (int j = 0; j < GetNbranch(); j++) {
            double mean = 0;
            double var = 0;
            for (int g = 0; g < GetNgene(); g++) {
                double tmp = blarray[g]->GetVal(j);
                mean += tmp;
                var += tmp * tmp;
            }
            mean /= GetNgene();
            var /= GetNgene();
            var -= mean * mean;
            tot += var;
        }
        tot /= GetNbranch();
        return tot;
    }

    //! return total log prob (over all genes and over all branches)
    double GetLogProb() const {
        double total = 0;
        for (int gene = 0; gene < GetNgene(); gene++) {
            total += blarray[gene]->GetLogProb();
        }
        return total;
    }

  private:
    int Ngene;
    const Tree &tree;
    const BranchSelector<double> &blmean;
    double shape;
    vector<GammaWhiteNoise *> blarray;
};

#endif
