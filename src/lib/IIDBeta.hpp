
#pragma once

#include "Array.hpp"
#include "CodonSuffStat.hpp"
#include "PoissonSuffStat.hpp"
#include "global/Random.hpp"

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

    //! return log prob of suff stats, as a function of mean and inverse
    //! concentration
    double GetLogProb(double mean, double invconc) const {
        double alpha = mean / invconc;
        double beta = (1 - mean) / invconc;
        return GetLogProb(alpha, beta);
    }

    //! return N (number of beta variates included in this suff stat)
    int GetN() const { return n; }
    //! return sum log(x_i)'s
    double GetSumLog0() const { return sumlog0; }
    //! return sum log(1-x_i)'s
    double GetSumLog1() const { return sumlog1; }

    template <class T>
    void serialization_interface(T &x) {
        x.add(sumlog0, sumlog1, n);
    }

  private:
    double sumlog0;
    double sumlog1;
    int n;
};

template <>
struct has_custom_serialization<BetaSuffStat> : std::true_type {};

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
    IIDBeta(int insize, const double &inmean, const double &ininvconc)
        : SimpleArray<double>(insize), mean(inmean), invconc(ininvconc) {
        Sample();
    }

    ~IIDBeta() {}

    //! return value of the alpha parameter
    double GetAlpha() const { return mean / invconc; }
    //! return value of the beta parameter
    double GetBeta() const { return (1.0 - mean) / invconc; }

    //! sample from prior
    void Sample() {
        for (int i = 0; i < GetSize(); i++) {
            (*this)[i] = Random::BetaSample(GetAlpha(), GetBeta());
        }
    }

    //! return log prior density of current values across array
    double GetLogProb() const {
        double total = 0;
        for (int i = 0; i < GetSize(); i++) { total += GetLogProb(i); }
        return total;
    }

    //! return log prior for entry i
    double GetLogProb(int i) const {
        return Random::logBetaDensity(GetVal(i), GetAlpha(), GetBeta());
    }

    //! add current values stored in array to BetaSuffStat given as argument
    void AddSuffStat(BetaSuffStat &suffstat) {
        for (int i = 0; i < GetSize(); i++) {
            suffstat.AddSuffStat(log(GetVal(i)), log(1 - GetVal(i)), 1);
        }
    }

  protected:
    const double &mean;
    const double &invconc;
};

template <>
struct has_custom_serialization<IIDBeta> : std::true_type {};
