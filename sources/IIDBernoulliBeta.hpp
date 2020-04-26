
#ifndef BERNBETA_H
#define BERNBETA_H

#include "Array.hpp"
#include "CodonSuffStat.hpp"
#include "MPIBuffer.hpp"
#include "PoissonSuffStat.hpp"
#include "Random.hpp"

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
    void Add(const BernoulliBetaSuffStat &from) {
        sumlog0 += from.GetSumLog0();
        sumlog1 += from.GetSumLog1();
        n0 += from.GetN0();
        n1 += from.GetN1();
    }

    //! (*this) += from, operator version
    BernoulliBetaSuffStat &operator+=(const BernoulliBetaSuffStat &from) {
        Add(from);
        return *this;
    }

    //! return object size, when put into an MPI buffer
    unsigned int GetMPISize() const { return 4; }

    //! put object into MPI buffer
    void MPIPut(MPIBuffer &buffer) const { buffer << sumlog0 << sumlog1 << n0 << n1; }

    //! read object from MPI buffer
    void MPIGet(const MPIBuffer &buffer) { buffer >> sumlog0 >> sumlog1 >> n0 >> n1; }

    //! read a BernoulliBetaSuffStat from MPI buffer and add it to this
    void Add(const MPIBuffer &buffer) {
        double temp;
        buffer >> temp;
        sumlog0 += temp;
        buffer >> temp;
        sumlog1 += temp;

        int tmp;
        buffer >> tmp;
        n0 += tmp;
        buffer >> tmp;
        n1 += tmp;
    }

    //! return log prob of suff stats, as a function of parameters of the mixture
    //! and the Beta distribution
    double GetLogProb(double pi, double alpha, double beta) const {
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

  private:
    double sumlog0;
    double sumlog1;
    int n0;
    int n1;
};

/**
 * \brief An array of IID random variables, from a mixture of a point mass at 0
 * (weight 1-pi) and a Beta distribution (weight pi)
 */

class IIDBernoulliBeta : public SimpleArray<double> {
  public:
    //! constructor, parameterized by size, weight parameter pi of the mixture,
    //! and parameters of the Beta (alpha and beta)
    IIDBernoulliBeta(int insize, double inpi, double inalpha, double inbeta)
        : SimpleArray<double>(insize), pi(inpi), alpha(inalpha), beta(inbeta) {
        Sample();
    }

    ~IIDBernoulliBeta() {}

    //! return weight of the Beta component of the mixture
    double GetPi() const { return pi; }
    //! return alpha parameter of the Beta distribution
    double GetAlpha() const { return alpha; }
    //! return beta parameter of the Beta distribution
    double GetBeta() const { return beta; }

    //! set weight of the Beta component of the mixture to new value
    void SetPi(double inpi) { pi = inpi; }
    //! set alpha parameter of the Beta distribution to new value
    void SetAlpha(double inalpha) { alpha = inalpha; }
    //! set beta parameter of the Beta distribution to new value
    void SetBeta(double inbeta) { beta = inbeta; }

    //! sample all entries from prior distribution
    void Sample() {
        for (int i = 0; i < GetSize(); i++) {
            if (Random::Uniform() < pi) {
                (*this)[i] = Random::BetaSample(alpha, beta);
            } else {
                (*this)[i] = 0;
            }
        }
    }

    void OffsetZeros(double in) {
        for (int i=0; i<GetSize(); i++) {
            if (!(*this)[i])    {
                (*this)[i] = 0.1;
            }
        }
    }
    //! get number of entries that are equal to 0
    int GetNullSet() const {
        int tot = 0;
        for (int i = 0; i < GetSize(); i++) {
            if (!GetVal(i)) {
                tot++;
            }
        }
        return tot;
    }

    //! return log probability summed over all entries
    double GetLogProb() const {
        double total = 0;
        for (int i = 0; i < GetSize(); i++) {
            total += GetLogProb(i);
        }
        return total;
    }

    //! return log probability for entry i
    double GetLogProb(int i) const {
        double ret = 0;
        if (!GetVal(i)) {
            ret += log(1 - pi);
        } else {
            ret += log(pi);
            ret += Random::logBetaDensity(GetVal(i), alpha, beta);
        }
        return ret;
    }

    //! add all entries to sufficient statistic
    void AddSuffStat(BernoulliBetaSuffStat &suffstat) {
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
        if (!tot) {
            return 0;
        }
        m1 /= tot;
        return m1;
    }

  protected:
    double pi;
    double alpha;
    double beta;
};

#endif
