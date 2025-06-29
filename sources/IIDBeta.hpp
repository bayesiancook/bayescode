
#pragma once

#include "Array.hpp"
#include "CodonSuffStat.hpp"
#include "MPIBuffer.hpp"
#include "PoissonSuffStat.hpp"
#include "Random.hpp"

/**
 * \brief A sufficient statistic for a collection of beta variates, as a
 * function of the parameters of the beta distribution
 *
 * Suppose you have x = (x_i)_i=1..N, iid Beta(a,b).
 * Then, p(X | a,b) can be expressed as a function of compact sufficient
 * statistics: sum log(x_i)'s, sum log(1-x_i)'s and N. BetaSuffStat implements
 * this idea, by providing methods for collecting these suff stats and returning
 * the log prob for a given value for the a and b parameters.
 */

class BetaSuffStat : public SuffStat {
  public:
    BetaSuffStat() {}
    ~BetaSuffStat() {}

    //! set suff stats to 0
    void Clear() {
        sumlog0 = 0;
        sumlog1 = 0;
        n = 0;
    }

    //! add the contribution of one beta variate (x) to this suffstat
    void AddSuffStat(double log0, double log1, int c = 1) {
        sumlog0 += log0;
        sumlog1 += log1;
        n += c;
    }

    //! (*this) += from
    void Add(const BetaSuffStat &from) {
        sumlog0 += from.GetSumLog0();
        sumlog1 += from.GetSumLog1();
        n += from.GetN();
    }

    //! (*this) += from, operator version
    BetaSuffStat &operator+=(const BetaSuffStat &from) {
        Add(from);
        return *this;
    }

    //! return object size, when put into an MPI buffer
    unsigned int GetMPISize() const { return 3; }

    //! put object into MPI buffer
    void MPIPut(MPIBuffer &buffer) const { buffer << sumlog0 << sumlog1 << n; }

    //! read object from MPI buffer
    void MPIGet(const MPIBuffer &buffer) { buffer >> sumlog0 >> sumlog1 >> n; }

    //! read a BetaSuffStat from MPI buffer and add it to this
    void Add(const MPIBuffer &buffer) {
        double temp;
        buffer >> temp;
        sumlog0 += temp;
        buffer >> temp;
        sumlog1 += temp;

        int tmp;
        buffer >> tmp;
        n += tmp;
    }

    //! return log prob of suff stats, as a function of mean and inverse
    //! concentration
    double GetMeanInvConcLogProb(double mean, double invconc) const {
        double alpha = mean / invconc;
        double beta = (1 - mean) / invconc;
        return GetLogProb(alpha, beta);
    }

    //! return log prob of suff stats, as a function of alpha and beta parameters
    double GetLogProb(double alpha, double beta) const {
        return n * (Random::logGamma(alpha + beta) - Random::logGamma(alpha) -
                    Random::logGamma(beta)) +
               (alpha - 1) * sumlog0 + (beta - 1) * sumlog1;
    }

    //! return N (number of beta variates included in this suff stat)
    int GetN() const { return n; }
    //! return sum log(x_i)'s
    double GetSumLog0() const { return sumlog0; }
    //! return sum log(1-x_i)'s
    double GetSumLog1() const { return sumlog1; }

  private:
    double sumlog0;
    double sumlog1;
    int n;
};

/**
 * \brief An array of IID beta random variables
 *
 * Note that the parameters of the beta distribution are given by copy (not by
 * ref) to the array. Thus, each time these parameters are modified during the
 * MCMC, the new values should be given to the array (using SetAlpha, SetBeta or
 * SetMeanInvConc).
 */

class IIDBeta : public SimpleArray<double> {
  public:
    //! constructor parameterized by array size, and alpha and beta parameters
    IIDBeta(int insize, double inalpha, double inbeta)
        : SimpleArray<double>(insize), alpha(inalpha), beta(inbeta) {
        Sample();
    }

    ~IIDBeta() {}

    //! return value of the alpha parameter
    double GetAlpha() const { return alpha; }
    //! return value of the beta parameter
    double GetBeta() const { return beta; }

    //! return expected mean (not actual mean)
    double GetMean() const { return alpha / (alpha + beta); }
    //! return inverse concentration
    double GetInvConcentration() const { return 1.0 / (alpha + beta); }

    //! get mean over the array
    double GetTrueMean() const {
        double m1 = 0;
        for (int i = 0; i < GetSize(); i++) {
            m1 += GetVal(i);
        }
        m1 /= GetSize();
        return m1;
    }

    //! get variance over the array
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

    //! set mean and inverse concentration (i.e. set alpha and beta, based on
    //! these two values)
    void SetMeanInvConc(double mean, double invconc) {
        double alpha = mean / invconc;
        double beta = (1 - mean) / invconc;
        SetAlpha(alpha);
        SetBeta(beta);
    }

    //! set alpha parameter to new value
    void SetAlpha(double inalpha) { alpha = inalpha; }
    //! set beta parameter to new value
    void SetBeta(double inbeta) { beta = inbeta; }

    //! sample from prior
    void Sample() {
        for (int i = 0; i < GetSize(); i++) {
            (*this)[i] = Random::BetaSample(alpha, beta);
        }
    }

    //! return log prior density of current values across array
    double GetLogProb() const {
        double total = 0;
        for (int i = 0; i < GetSize(); i++) {
            total += GetLogProb(i);
        }
        return total;
    }

    //! return log prior for entry i
    double GetLogProb(int i) const { return Random::logBetaDensity(GetVal(i), alpha, beta); }

    //! add current values stored in array to BetaSuffStat given as argument
    void AddSuffStat(BetaSuffStat &suffstat) {
        for (int i = 0; i < GetSize(); i++) {
            suffstat.AddSuffStat(log(GetVal(i)), log(1 - GetVal(i)), 1);
        }
    }

  protected:
    double alpha;
    double beta;
};


class BranchIIDBeta : public SimpleBranchArray<double> {
  public:
    BranchIIDBeta(const Tree &intree, double inalpha, double inbeta)
        : SimpleBranchArray<double>(intree), alpha(inalpha), beta(inbeta) {
        Sample();
    }

    ~BranchIIDBeta() {}

    //! return value of the alpha parameter
    double GetAlpha() const { return alpha; }
    //! return value of the beta parameter
    double GetBeta() const { return beta; }

    //! return expected mean (not actual mean)
    double GetMean() const { return alpha / (alpha + beta); }
    //! return inverse concentration
    double GetInvConcentration() const { return 1.0 / (alpha + beta); }

    //! set mean and inverse concentration (i.e. set alpha and beta, based on
    //! these two values)
    void SetMeanInvConc(double mean, double invconc) {
        double alpha = mean / invconc;
        double beta = (1 - mean) / invconc;
        SetAlpha(alpha);
        SetBeta(beta);
    }

    //! set alpha parameter to new value
    void SetAlpha(double inalpha) { alpha = inalpha; }
    //! set beta parameter to new value
    void SetBeta(double inbeta) { beta = inbeta; }

    //! set all entries equal to inval
    void SetAllBranches(double inval) {
        for (int i = 0; i < GetNbranch(); i++) {
            (*this)[i] = inval;
        }
    }

    //! sample all entries from prior
    void Sample() {
        for (int i = 0; i < GetNbranch(); i++) {
            (*this)[i] = Random::BetaSample(alpha, beta);
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
    double GetLogProb(int index) { return Random::logBetaDensity(GetVal(index), alpha, beta); }

    //! add current values stored in array to BetaSuffStat given as argument
    void AddSuffStat(BetaSuffStat &suffstat) {
        for (int i = 0; i < GetNbranch(); i++) {
            suffstat.AddSuffStat(log(GetVal(i)), log(1 - GetVal(i)), 1);
        }
    }
    //! get sum over all entries (name is rather specialized... could change..)
    double GetTotalLength() const {
        double m1 = 0;
        for (int i = 0; i < GetNbranch(); i++) {
            m1 += GetVal(i);
        }
        return m1;
    }

    //! get mean over the array
    double GetTrueMean() const {
        double m1 = 0;
        for (int i = 0; i < GetNbranch(); i++) {
            m1 += GetVal(i);
        }
        m1 /= GetNbranch();
        return m1;
    }

    //! get variance over the array
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
    double alpha;
    double beta;
};

