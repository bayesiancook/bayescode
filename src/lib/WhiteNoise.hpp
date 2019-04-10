#pragma once

#include "BranchArray.hpp"
#include "PoissonSuffStat.hpp"
#include "global/Random.hpp"

/**
 * \brief A tree-structured branch-wise array of Gamma variables, with
 * branch-specific means but same invshape parameter
 */
class GammaWhiteNoise : public SimpleBranchArray<double> {
  public:
    // mode 0 : ugam
    // mode 1 : wn
    GammaWhiteNoise(const Tree &intree, const BranchSelector<double> &inblmean,
        const double &ininvshape, int inmode = 0)
        : SimpleBranchArray<double>(intree), blmean(inblmean), invshape(ininvshape) {
        mode = inmode;
        Sample();
    }

    ~GammaWhiteNoise() {}

    GammaWhiteNoise &operator=(const GammaWhiteNoise &from) {
        array = from.array;
        return *this;
    }

    void SetMode(int inmode) { mode = inmode; }

    double GetShape() const { return 1.0 / invshape; }

    //! return the shape (alpha) parameter for branch i
    double GetAlpha(int i) const {
        if (mode) { return blmean.GetVal(i) / invshape; }
        return 1.0 / invshape;
    }

    //! return the scale (beta) parameter for branch i
    double GetBeta(int i) const {
        if (mode) { return 1.0 / invshape; }
        return 1.0 / invshape / blmean.GetVal(i);
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
            (*this)[i] = Random::GammaSample(
                GetAlpha(i) + suffstat.GetCount(), GetBeta(i) + suffstat.GetBeta());
            if ((*this)[i] == 0) {
                std::cerr << "gibbs null bl : " << GetAlpha(i) << '\t' << GetBeta(i) << '\t'
                          << invshape << '\t' << blmean.GetVal(i) << '\n';
                (*this)[i] = 0.001;
            }
        }
    }

    //! resample entries based on a BranchArray of PoissonSuffStat
    void ResampleEmptyBranches(const PoissonSuffStatBranchArray &suffstatarray) {
        for (int i = 0; i < GetNbranch(); i++) {
            if (suffstatarray.GetVal(i).GetBeta() == 0) {
                (*this)[i] = Random::GammaSample(GetAlpha(i), GetBeta(i));
                if ((*this)[i] == 0) {
                    std::cerr << "empty null bl : " << GetAlpha(i) << '\t' << GetBeta(i) << '\t'
                              << invshape << '\t' << blmean.GetVal(i) << '\n';
                    (*this)[i] = 0.001;
                }
            }
        }
    }

    //! return total log prob summed over all entries
    double GetLogProb() const {
        double total = 0;
        for (int i = 0; i < GetNbranch(); i++) { total += GetLogProb(i); }
        return total;
    }

    //! return log prob for one entry
    double GetLogProb(int index) const {
        return Random::logGammaDensity(GetVal(index), GetAlpha(index), GetBeta(index));
    }

    double GetTotalLength() const {
        double m1 = 0;
        for (int i = 0; i < GetNbranch(); i++) { m1 += GetVal(i); }
        return m1;
    }

  protected:
    const BranchSelector<double> &blmean;
    const double &invshape;
    int mode;
};

template <>
struct has_custom_serialization<GammaWhiteNoise> : std::true_type {};

/**
 * \brief An array of GammaWhiteNoise
 *
 * This is an Array of BranchArray's. It is used for implementing modulations of
 * branch-lengths across genes, while implementing some basic shrinkage.
 * Specifically, for a given branch j,
 * the branch-lengths across genes are gamma of mean blmean[j] and shape
 * parameter uniform across all branches.
 */

class GammaWhiteNoiseArray : public SimpleArray<GammaWhiteNoise> {
  public:
    //! constructor: parameterized by the number of genes, the tree, the means
    //! over branches and the shape parameter
    GammaWhiteNoiseArray(int inNgene, const Tree &intree, const BranchSelector<double> &inblmean,
        const double &ininvshape)
        : SimpleArray(inNgene, GammaWhiteNoise(intree, inblmean, ininvshape)) {}

    ~GammaWhiteNoiseArray() {}

    //! return total number of genes
    int GetNgene() const { return size(); }

    //! return total number of branches of the underlying tree
    int GetNbranch() const { return array[0].GetNbranch(); }

    //! return mean tree length over genes
    double GetMeanLength() const {
        double tot = 0;
        for (int j = 0; j < GetNbranch(); j++) {
            double mean = 0;
            for (int g = 0; g < GetNgene(); g++) {
                double tmp = array[g].GetVal(j);
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
                double tmp = array[g].GetVal(j);
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
        for (int gene = 0; gene < GetNgene(); gene++) { total += array[gene].GetLogProb(); }
        return total;
    }
};

template <>
struct has_custom_serialization<GammaWhiteNoiseArray> : std::true_type {};
