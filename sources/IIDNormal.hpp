
#ifndef IIDMORMAL_H
#define IIDNORMAL_H

#include "Array.hpp"
#include "BranchArray.hpp"
#include "MPIBuffer.hpp"
#include "Random.hpp"
#include "SuffStat.hpp"

/**
 * \brief An array of IID normal random variables
 *
 * Note that the mean and variance parameters are given by copy (not by ref) to
 * the array. Thus, each time these parameters are modified during
 * the MCMC, the new values should be given to the array (using the SetMean and SetVar methods);
 */

class IIDNormal : public SimpleArray<double> {
  public:
    //! constructor specifies array size and initial mean and var parameters
    IIDNormal(int insize, double inmean, double invar)
        : SimpleArray<double>(insize), mean(inmean), var(invar) {
        Sample();
    }

    ~IIDNormal() {}

    //! return mean parameter
    double GetMean() const { return mean; }

    //! return var parameter
    double GetVar() const { return var; }

    //! set mean parameter to a new value
    void SetMean(double inmean) { mean = inmean; }

    //! set var parameter to a new value
    void SetVar(double invar) { var = invar; }

    //! sample all entries, given current mean and var params
    void Sample() {
        for (int i = 0; i < GetSize(); i++) {
            (*this)[i] = Random::NormalSample(mean, var);
        }
    }

    //! return total log prob (summed over the array) given current mean and
    //! var params
    double GetLogProb() const {
        double total = 0;
        for (int i = 0; i < GetSize(); i++) {
            total += GetLogProb(i);
        }
        return total;
    }

    //! return log prob of one specific entry
    double GetLogProb(int index) const {
        return Random::logNormalDensity(GetVal(index), mean, var);
    }

    //! resample all entries for which occupancy[i] == 0 from the prior (from a
    //! Normal(mean,var))
    void PriorResample(const Selector<int> &occupancy) {
        for (int i = 0; i < GetSize(); i++) {
            if (!occupancy.GetVal(i)) {
                (*this)[i] = Random::NormalSample(mean, var);
            }
        }
    }

  protected:
    double mean;
    double var;
};

/**
 * \brief A tree-structured branch-wise array of iid Normal variables
 * (tree-stuctured version of IIDNormal)
 *
 * One should be careful about the fact that the mean and var parameters are
 * given by copy (not by ref) to the array. Thus, each time the mean and var
 * parameters are modified during the MCMC, the new values should be given to
 * the array (using the SetMean and SetVar methods).
 */

class BranchIIDNormal : public SimpleBranchArray<double> {
  public:
    BranchIIDNormal(const Tree &intree, double inmean, double invar)
        : SimpleBranchArray<double>(intree), mean(inmean), var(invar) {
        Sample();
    }

    ~BranchIIDNormal() {}

    double GetMean() const { return mean; }
    double GetVar() const { return var; }

    void SetMean(double inmean) { mean = inmean; }

    void SetVar(double invar) { var = invar; }

    //! set all entries equal to inval
    void SetAllBranches(double inval) {
        for (int i = 0; i < GetNbranch(); i++) {
            (*this)[i] = inval;
        }
    }

    //! sample all entries from prior
    void Sample() {
        for (int i = 0; i < GetNbranch(); i++) {
            (*this)[i] = Random::NormalSample(mean, var);
        }
    }

    //! get total log prob summed over all branches
    double GetLogProb() {
        double total = 0;
        for (int i = 0; i < GetNbranch(); i++) {
            total += GetLogProb(i);
        }
        return total;
    }

    //! get log prob for a given branch
    double GetLogProb(int index) { return Random::logNormalDensity(GetVal(index), mean, var); }

  protected:
    double mean;
    double var;
};

class NormalSuffStat : public SuffStat {
  public:
    NormalSuffStat() {}
    ~NormalSuffStat() {}

    //! set suff stats to 0
    void Clear() {
        m1 = 0;
        m2 = 0;
        n = 0;
    }

    //! add the contribution of one normal variate (x) to this suffstat
    void AddSuffStat(double x)  {
        m1 += x;
        m2 += x*x;
        n ++;
    }

    //! (*this) += from
    void Add(const NormalSuffStat &from) {
        m1 += from.m1;
        m2 += from.m2;
        n += from.n;
    }

    //! (*this) += from, operator version
    NormalSuffStat &operator+=(const NormalSuffStat &from) {
        Add(from);
        return *this;
    }

    //! get suff stats from an IIDNormal array
    void AddSuffStat(const IIDNormal &array) {
        for (int i = 0; i < array.GetSize(); i++) {
            AddSuffStat(array.GetVal(i));
        }
    }

    //! get suff stats from entries of an IIDNormal array such that occupancy[i] !=
    //! 0
    void AddSuffStat(const IIDNormal &array, const Selector<int> &occupancy) {
        for (int i = 0; i < array.GetSize(); i++) {
            if (occupancy.GetVal(i)) {
                AddSuffStat(array.GetVal(i));
            }
        }
    }

    //! get suff stats from a BranchIIDNormal array
    void AddSuffStat(const BranchIIDNormal &array) {
        for (int i = 0; i < array.GetNbranch(); i++) {
            AddSuffStat(array.GetVal(i));
        }
    }

    //! return object size, when put into an MPI buffer
    unsigned int GetMPISize() const { return 3; }

    //! put object into MPI buffer
    void MPIPut(MPIBuffer &buffer) const { buffer << m1 << m2 << n; }

    //! read object from MPI buffer
    void MPIGet(const MPIBuffer &buffer) { buffer >> m1 >> m2 >> n; }

    //! read a NormalSuffStat from MPI buffer and add it to this
    void Add(const MPIBuffer &buffer) {
        double temp;
        buffer >> temp;
        m1 += temp;
        buffer >> temp;
        m2 += temp;

        int tmp;
        buffer >> tmp;
        n += tmp;
    }

    //! return log prob, as a function of the given shape and scale parameters
    double GetLogProb(double mean, double var) const    {
        return -0.5 * (n*log(2*Pi*var) + (m2 - 2*mean*m1 + n*mean*mean)/var);
    }

    //! return sum of normal variates contributing to the suff stat
    double GetSum() const { return m1; }
    //! return sum of squares of normal variates contributing to the suff stat
    double GetSumOfSquares() const { return m2; }
    //! return N, total number of normal variates contributing to the suff stat
    int GetN() const { return n; }

  private:
    double m1;
    double m2;
    int n;
};
#endif
