
#pragma once

#include "Array.hpp"
#include "CodonSuffStat.hpp"
#include "PoissonSuffStat.hpp"
#include "global/Random.hpp"

/**
 * \brief A sufficient statistic for a collection of random variables, from a
 * mixture of a point mass at 0 and a Beta distribution (see IIDBernoulliBeta)
 */

class BernoulliBetaSuffStat : public SuffStat {
  public:
    BernoulliBetaSuffStat() {}
    ~BernoulliBetaSuffStat() {}

    //! set suff stats to 0
    void Clear() {
        sumlog0 = 0;
        sumlog1 = 0;
        n0 = 0;
        n1 = 0;
    }

    //! add the contribution of one beta variate equal to 0 to this suffstat
    void AddNullSuffStat(int c = 1) { n0 += c; }

    //! add the contribution of one beta variate > 0 to this suffstat
    void AddPosSuffStat(double log0, double log1, int c = 1) {
        sumlog0 += log0;
        sumlog1 += log1;
        n1 += c;
    }

    //! (*this) += from
    void Add(const BernoulliBetaSuffStat& from) {
        sumlog0 += from.GetSumLog0();
        sumlog1 += from.GetSumLog1();
        n0 += from.GetN0();
        n1 += from.GetN1();
    }

    //! (*this) += from, operator version
    BernoulliBetaSuffStat& operator+=(const BernoulliBetaSuffStat& from) {
        Add(from);
        return *this;
    }

    //! return log prob of suff stats, as a function of parameters of the mixture
    //! and the Beta distribution
    double GetLogProb(double pi, double mean, double invconc) const {
        double alpha = mean / invconc;
        double beta = (1.0 - mean) / invconc;
        double logbern = n0 * log(1 - pi) + n1 * log(pi);
        double logbeta = n1 * (Random::logGamma(alpha + beta) - Random::logGamma(alpha) -
                                  Random::logGamma(beta)) +
                         (alpha - 1) * sumlog0 + (beta - 1) * sumlog1;
        return logbern + logbeta;
    }

    //! get number of contributing variables that are == 0
    int GetN0() const { return n0; }
    //! get number of contributing variables that are > 0
    int GetN1() const { return n1; }
    //! get total sum of log(1-x_i), for all x_i's that are > 0
    double GetSumLog0() const { return sumlog0; }
    //! get total sum of log(x_i), for all x_i's that are > 0
    double GetSumLog1() const { return sumlog1; }

    template <class T>
    void serialization_interface(T& x) {
        x.add(sumlog0, sumlog1, n0, n1);
    }

  private:
    double sumlog0;
    double sumlog1;
    int n0;
    int n1;
};

template <>
struct has_custom_serialization<BernoulliBetaSuffStat> : std::true_type {};

/**
 * \brief An array of IID random variables, from a mixture of a point mass at 0
 * (weight 1-pi) and a Beta distribution (weight pi)
 */

class IIDBernoulliBeta : public SimpleArray<double> {
  public:
    //! constructor, parameterized by size, weight parameter pi of the mixture,
    //! and parameters of the Beta (alpha and beta)
    IIDBernoulliBeta(int insize, const double& inpi, const double& inmean, const double& ininvconc)
        : SimpleArray<double>(insize), pi(inpi), mean(inmean), invconc(ininvconc) {
        Sample();
    }

    ~IIDBernoulliBeta() {}

    //! return alpha parameter of the Beta distribution
    double GetAlpha() const { return mean / invconc; }
    //! return beta parameter of the Beta distribution
    double GetBeta() const { return (1 - mean) / invconc; }

    //! sample all entries from prior distribution
    void Sample() {
        for (int i = 0; i < GetSize(); i++) {
            if (Random::Uniform() < pi) {
                (*this)[i] = Random::BetaSample(GetAlpha(), GetBeta());
            } else {
                (*this)[i] = 0;
            }
        }
    }

    //! get number of entries that are equal to 0
    int GetNullSet() const {
        int tot = 0;
        for (int i = 0; i < GetSize(); i++) {
            if (!GetVal(i)) { tot++; }
        }
        return tot;
    }

    //! return log probability summed over all entries
    double GetLogProb() const {
        double total = 0;
        for (int i = 0; i < GetSize(); i++) { total += GetLogProb(i); }
        return total;
    }

    //! return log probability for entry i
    double GetLogProb(int i) const {
        double ret = 0;
        if (!GetVal(i)) {
            ret += log(1 - pi);
        } else {
            ret += log(pi);
            ret += Random::logBetaDensity(GetVal(i), GetAlpha(), GetBeta());
        }
        return ret;
    }

    //! add all entries to sufficient statistic
    void AddSuffStat(BernoulliBetaSuffStat& suffstat) {
        for (int i = 0; i < GetSize(); i++) {
            if (!GetVal(i)) {
                suffstat.AddNullSuffStat(1);
            } else {
                suffstat.AddPosSuffStat(log(GetVal(i)), log(1 - GetVal(i)), 1);
            }
        }
    }

    //! get mean of all of those entries that are not 0
    double GetPosMean() const {
        int tot = 0;
        double m1 = 0;
        for (int i = 0; i < GetSize(); i++) {
            if (GetVal(i)) {
                m1 += GetVal(i);
                tot++;
            }
        }
        if (!tot) { return 0; }
        m1 /= tot;
        return m1;
    }

  protected:
    const double& pi;
    const double& mean;
    const double& invconc;
};

template <>
struct has_custom_serialization<IIDBernoulliBeta> : std::true_type {};
