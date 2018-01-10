
#ifndef GAMMASUFFSTAT_H
#define GAMMASUFFSTAT_H

#include "IIDGamma.hpp"

/**
 * \brief A sufficient statistic for a collection of gamma variates, as a function of the shape and scale parameters
 *
 * Suppose you have x = (x_i)_i=1..N, iid Gamma(shape,scale).
 * Then, p(X | shape,scale) can be expressed as a function of compact sufficient statistics: sum x_i's, sum log(x_i)'s and N.
 * GammaSuffStat implements this idea, by providing methods for collecting these suff stats and returning the log prob for a
 * given
 * value for the shape and scale parameters.
 */

class GammaSuffStat : public SuffStat, public tc::Component {
  public:
    GammaSuffStat() {}
    ~GammaSuffStat() {}

    //! set suff stats to 0
    void Clear() {
        sum = 0;
        sumlog = 0;
        n = 0;
    }

    //! add the contribution of one gamma variate (x) to this suffstat
    void AddSuffStat(double x, double logx, int c = 1) {
        sum += x;
        sumlog += logx;
        n += c;
    }

    //! (*this) += from
    void Add(const GammaSuffStat& from) {
        sum += from.GetSum();
        sumlog += from.GetSumLog();
        n += from.GetN();
    }

    //! (*this) += from, operator version
    GammaSuffStat& operator+=(const GammaSuffStat& from) {
        Add(from);
        return *this;
    }

    //! get suff stats from an IIDGamma array
    void AddSuffStat(const IIDGamma& array) {
        for (int i = 0; i < array.GetSize(); i++) {
            AddSuffStat(array.GetVal(i), log(array.GetVal(i)));
        }
    }

    //! get suff stats from entries of an IIDGamma array such that occupancy[i] != 0
    void AddSuffStat(const IIDGamma& array, const Selector<int>& occupancy) {
        for (int i = 0; i < array.GetSize(); i++) {
            if (occupancy.GetVal(i)) {
                AddSuffStat(array.GetVal(i), log(array.GetVal(i)));
            }
        }
    }

    //! get suff stats from entries of an IIDGamma array such that probarray[i] != 0
    void AddSuffStat(const IIDGamma& array, const Selector<double>& probarray) {
        for (int i = 0; i < array.GetSize(); i++) {
            if (probarray.GetVal(i)) {
                AddSuffStat(array.GetVal(i), log(array.GetVal(i)));
            }
        }
    }

    //! get suff stats from a BranchIIDGamma array
    void AddSuffStat(const BranchIIDGamma& array) {
        for (int i = 0; i < array.GetNbranch(); i++) {
            AddSuffStat(array.GetVal(i), log(array.GetVal(i)));
        }
    }

    //! get suff stats from GammaWhiteNoise (!! valid only if all blmean[i]'s are the same across all branches)
    void AddSuffStat(const GammaWhiteNoise& array) {
        for (int i = 0; i < array.GetNbranch(); i++) {
            AddSuffStat(array.GetVal(i), log(array.GetVal(i)));
        }
    }

    //! return object size, when put into an MPI buffer
    unsigned int GetMPISize() const { return 3; }

    //! put object into MPI buffer
    void MPIPut(MPIBuffer& buffer) const { buffer << sum << sumlog << n; }

    //! read object from MPI buffer
    void MPIGet(const MPIBuffer& buffer) { buffer >> sum >> sumlog >> n; }

    //! read a GammaSuffStat from MPI buffer and add it to this
    void Add(const MPIBuffer& buffer) {
        double temp;
        buffer >> temp;
        sum += temp;
        buffer >> temp;
        sumlog += temp;

        int tmp;
        buffer >> tmp;
        n += tmp;
    }

    //! return log prob, as a function of the given shape and scale parameters
    double GetLogProb(double shape, double scale) const {
        return n * (shape * log(scale) - Random::logGamma(shape)) + (shape - 1) * sumlog - scale * sum;
    }

    //! return sum x_i's
    double GetSum() const { return sum; }
    //! return sum log x_i's
    double GetSumLog() const { return sumlog; }
    //! return N, total number of gamma variates contributing to the suff stat
    int GetN() const { return n; }

  private:
    double sum;
    double sumlog;
    int n;
};

/**
 * \brief A tree-structured branch-wise array of gamma sufficient statistics
 *
 * Useful for gene-specific branch lengths (see for instance MultiGeneCodonM2aModel).
 * In this model,
 * for a given branch of the tree,
 * genes have differing lengths, which are iid Gamma, of shape and scale that are branch-specific.
 * Thus, it is useful to collect suff stats across genes, and do this branchwise (storing the result into a branch-wise
 * array).
 * This array of suffstats can then be used to do fast MCMC moves on the hyperparameters of the distribution of
 * branch lengths across genes.
 */

class GammaSuffStatBranchArray : public SimpleBranchArray<GammaSuffStat> {
  public:
    GammaSuffStatBranchArray(const Tree& intree) : SimpleBranchArray<GammaSuffStat>(intree) {}
    ~GammaSuffStatBranchArray() {}

    //! member-wise addition between the two arrays ((*this) += from)
    void Add(const GammaSuffStatBranchArray& from) {
        for (int i = 0; i < GetNbranch(); i++) {
            (*this)[i].Add(from.GetVal(i));
        }
    }

    //! member-wise addition between the two arrays, operator version
    GammaSuffStatBranchArray& operator+=(const GammaSuffStatBranchArray& from) {
        Add(from);
        return *this;
    }

    //! get suff stats from a GammaWhiteNoise
    void AddSuffStat(const GammaWhiteNoise& array) {
        for (int i = 0; i < array.GetNbranch(); i++) {
            (*this)[i].AddSuffStat(array.GetVal(i), log(array.GetVal(i)));
        }
    }

    //! get suff stats from a GammaWhiteNoiseArray
    void AddSuffStat(GammaWhiteNoiseArray& array) {
        for (int gene = 0; gene < array.GetNgene(); gene++) {
            AddSuffStat(array.GetVal(gene));
        }
    }

    //! return array size, when put into an MPI buffer
    unsigned int GetMPISize() const { return 3 * GetNbranch(); }

    //! put array into MPI buffer
    void MPIPut(MPIBuffer& buffer) const {
        for (int i = 0; i < GetNbranch(); i++) {
            buffer << GetVal(i);
        }
    }

    //! read array from MPI buffer
    void MPIGet(const MPIBuffer& buffer) {
        for (int i = 0; i < GetNbranch(); i++) {
            buffer >> (*this)[i];
        }
    }

    //! read from MPI buffer and add to current array
    void Add(const MPIBuffer& buffer) {
        for (int i = 0; i < GetNbranch(); i++) {
            (*this)[i].Add(buffer);
        }
    }

    //! set all suff stats to 0
    void Clear() {
        for (int i = 0; i < GetNbranch(); i++) {
            (*this)[i].Clear();
        }
    }

    //! get total log prob
    double GetLogProb(const BranchSelector<double>& blmean, double invshape) const {
        double total = 0;
        for (int i = 0; i < GetNbranch(); i++) {
            double shape = 1.0 / invshape;
            double scale = 1.0 / blmean.GetVal(i);
            total += GetVal(i).GetLogProb(shape, scale);
        }
        return total;
    }
};

#endif
