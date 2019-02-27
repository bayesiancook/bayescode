#pragma once

#include <cmath>
#include "Array.hpp"
#include "BranchArray.hpp"
#include "PhyloProcess.hpp"
#include "SuffStat.hpp"

/**
 * \brief A Poisson-like sufficient statistic
 *
 * This sufficient statistic deals with all cases where the probability of the
 * variable(s) of interest (say X), as a function of some rate parameter r
 * (positive real number), can be written as P(X | r) \propto r^count
 * exp(-beta*r), for some suff stats count (integer) and beta (positive real
 * number). When several independent variables are entailed by the generic X
 * mentioned above, the suff stats just have to be summed over all items
 * (separately for count and beta).
 *
 * As an example, the probability of the substitution histories over all sites,
 * for a given branch of the tree, as a function of the branch length l, can be
 * written as follows: p(S | l) = K l^count exp(-beta*l), where count is the
 * total number of substitution events along the branch (across all sites), beta
 * is the average rate away from the current state (averaged over all paths) and
 * K is some normalization constant (that may depend on other parameters than
 * l). Thus, conditional on the current substitution histories, the sufficient
 * statistics count and beta can be first computed, and then MCMC moves on
 * branch lengths can be done based on the knowledge of these two numbers. Other
 * examples include the probability of substitution histories as a function of
 * omega = dN/dS (in that case, the count suff stat is the number of
 * non-synonymous substitutions only). In these two cases, the suffstats have to
 * be calculated based on the substitution histories (see
 * BranchSitePath::AddLengthSuffStat) and then summed over branches, sites, or
 * any other more subtle pattern, depending on the exact structure of the model.
 *
 * PoissonSuffStat implements this general idea:
 * collecting count and beta suff stats across all relevant items
 * and providing methods for calculating the probability, as a function of the
 * rate parameter.
 */

class PoissonSuffStat : public SuffStat {
  public:
    PoissonSuffStat() { Clear(); }
    ~PoissonSuffStat() {}

    //! set count and beta to 0
    void Clear() { count = beta = 0; }

    void IncrementCount() { count++; }

    void AddCount(int in) { count += in; }

    void AddBeta(double in) { beta += in; }

    void AddSuffStat(int incount, double inbeta) {
        count += incount;
        beta += inbeta;
    }

    void Add(const PoissonSuffStat &from) {
        count += from.GetCount();
        beta += from.GetBeta();
    }

    PoissonSuffStat &operator+=(const PoissonSuffStat &from) {
        Add(from);
        return *this;
    }

    //! write structure into generic output stream
    void ToStream(std::ostream &os) const { os << count << '\t' << beta << '\n'; }

    //! read structure from generic input stream
    void FromStream(std::istream &is) { is >> count >> beta; }

    int GetCount() const { return count; }

    double GetBeta() const { return beta; }

    //! return the log probability as a function of the rate: essentially
    //! count*log(rate) - beta*rate
    double GetLogProb(double rate) const { return count * log(rate) - beta * rate; }

    //! return the log of the marginal probability when the rate is from a gamma
    //! distribution
    double GetMarginalLogProb(double mean, double invshape) const {
        double shape = 1.0 / invshape;
        double scale = shape / mean;
        return shape * log(scale) - Random::logGamma(shape) - (shape + count) * log(scale + beta) +
               Random::logGamma(shape + count);
    }

    template <class T>
    void serialization_interface(T &x) {
        x.add(count, beta);
    }

    // protected:
    int count;
    double beta;
};

template <>
struct has_custom_serialization<PoissonSuffStat> : std::true_type {};

inline std::ostream &operator<<(std::ostream &os, const PoissonSuffStat &suffstat) {
    suffstat.ToStream(os);
    return os;
}

inline std::istream &operator>>(std::istream &is, PoissonSuffStat &suffstat) {
    suffstat.FromStream(is);
    return is;
}

/**
 * \brief An array of Poisson sufficient statistics
 *
 * This class provides an interface for dealing with cases where
 * each item has a different rate parameter (e.g. rates r_i's or omega_i's
 * across sites). In that case, we need a distinct PoissonSuffStat for each
 * item. Each suff stat might still be summing over other subcases (e.g. across
 * all branches for a given site)
 */

class PoissonSuffStatArray : public SimpleArray<PoissonSuffStat> {
  public:
    //! constructor with array size
    PoissonSuffStatArray(int insize) : SimpleArray<PoissonSuffStat>(insize) {}
    ~PoissonSuffStatArray() {}

    //! set suff stats to 0
    void Clear() {
        for (int i = 0; i < GetSize(); i++) { (*this)[i].Clear(); }
    }

    //! add path sufficient statistics for resampling branch lengths from
    //! PhyloProcess
    void AddRatePathSuffStat(const PhyloProcess &process) { process.AddRateSuffStat(*this); }

    //! member-wise addition: essentially (*this)[i] += from[i], for
    //! i=0..GetSize()-1
    void Add(const PoissonSuffStatArray &from) {
        for (int i = 0; i < GetSize(); i++) { (*this)[i].Add(from.GetVal(i)); }
    }

    //! member-wise addition, operator version
    PoissonSuffStatArray &operator+=(const PoissonSuffStatArray &from) {
        Add(from);
        return *this;
    }

    //! \brief get logprob, based on an array of rates (of same size)
    //!
    //! specifically, total log prob such as returned by the function is sum_i log
    //! p((*this)[i] | ratearray[i])
    double GetLogProb(const Array<double> &ratearray) const {
        double total = 0;
        for (int i = 0; i < GetSize(); i++) { total += GetVal(i).GetLogProb(ratearray.GetVal(i)); }
        return total;
    }

    //! get marginal log prob, based on an array of rates that are iid from a
    //! gamma of given shape and scale parameters
    double GetMarginalLogProb(double mean, double invshape) const {
        double total = 0;
        for (int i = 0; i < GetSize(); i++) {
            total += GetVal(i).GetMarginalLogProb(mean, invshape);
        }
        return total;
    }
};

/**
 * \brief A tree-structured branch-wise array of Poisson sufficient statistics
 *
 * This class provides an interface for dealing with cases where
 * each branch has a different rate (or, for that matter, length) parameter
 * In that case, we need a distinct PoissonSuffStat for each item (for each
 * branch). Each suff stat might still be summing over subcases (e.g. across all
 * sites for a given branch)
 */

class PoissonSuffStatBranchArray : public SimpleBranchArray<PoissonSuffStat> {
  public:
    //! constructor parameterized by underlying treee
    PoissonSuffStatBranchArray(const Tree &intree) : SimpleBranchArray<PoissonSuffStat>(intree) {}
    ~PoissonSuffStatBranchArray() {}

    //! set all suff stats to 0
    void Clear() {
        for (int i = 0; i < GetNbranch(); i++) { (*this)[i].Clear(); }
    }

    //! get total log prob, based on an array of branch-specific rates (or
    //! lengths)
    double GetLogProb(const BranchSelector<double> &ratearray) const {
        double total = 0;
        for (int i = 0; i < GetNbranch(); i++) {
            total += GetVal(i).GetLogProb(ratearray.GetVal(i));
        }
        return total;
    }

    //! member-wise addition (*this)[i] = from[i] for all i=0..GetNbranch()-1
    void Add(const PoissonSuffStatBranchArray &from) {
        for (int i = 0; i < GetNbranch(); i++) { (*this)[i].Add(from.GetVal(i)); }
    }

    //! member-wise addition: operator version
    PoissonSuffStatBranchArray &operator+=(const PoissonSuffStatBranchArray &from) {
        Add(from);
        return *this;
    }

    //! add path sufficient statistics for resampling branch lengths from
    //! PhyloProcess
    void AddLengthPathSuffStat(const PhyloProcess &process) { process.AddLengthSuffStat(*this); }

    //! get total (summed) marginal log prob integrated over branch-specific rates
    //! (or lengths) iid from a gamma(shape,scale)
    double GetMarginalLogProb(double mean, double invshape) const {
        double total = 0;
        for (int i = 0; i < GetNbranch(); i++) {
            total += GetVal(i).GetMarginalLogProb(mean, invshape);
        }
        return total;
    }

    //! push two C-style arrays of counts and betas into this array (erase current
    //! values)
    void Push(int *count, double *beta) const {
        for (int i = 0; i < GetNbranch(); i++) {
            count[i] = GetVal(i).GetCount();
            beta[i] = GetVal(i).GetBeta();
        }
    }

    //! push-and-add two C-style arrays of counts and betas into this array
    //! (member-wise addition)
    void Add(const int *count, const double *beta) {
        for (int i = 0; i < GetNbranch(); i++) {
            (*this)[i].AddCount(count[i]);
            (*this)[i].AddBeta(beta[i]);
        }
    }
};

template <>
struct has_custom_serialization<PoissonSuffStatBranchArray> : std::true_type {};

class PoissonSuffStatTreeArray : public Array<PoissonSuffStatBranchArray> {
  public:
    //! constructor, parameterized by underlying tree and size (number of genes)
    PoissonSuffStatTreeArray(const Tree &intree, int insize)
        : tree(intree), size(insize), array(insize, (PoissonSuffStatBranchArray *)0) {
        for (int i = 0; i < GetSize(); i++) { array[i] = new PoissonSuffStatBranchArray(tree); }
    }

    ~PoissonSuffStatTreeArray() {
        for (int i = 0; i < GetSize(); i++) { delete[] array[i]; }
    }

    int GetSize() const override { return size; }

    const PoissonSuffStatBranchArray &GetVal(int i) const override { return *array[i]; }

    PoissonSuffStatBranchArray &operator[](int i) override { return *array[i]; }

    //! clear all suff stats
    void Clear() {
        for (int i = 0; i < GetSize(); i++) { array[i]->Clear(); }
    }

  private:
    const Tree &tree;
    int size;
    std::vector<PoissonSuffStatBranchArray *> array;
};
