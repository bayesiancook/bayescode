
#ifndef IIDMORMAL_H
#define IIDNORMAL_H

#include "Array.hpp"
#include "BranchArray.hpp"
#include "MPIBuffer.hpp"
#include "Random.hpp"

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

#endif
