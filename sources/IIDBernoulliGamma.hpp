
#ifndef BERNGAMMA_H
#define BERNGAMMA_H

#include "Array.hpp"
#include "CodonSuffStat.hpp"
#include "MPIBuffer.hpp"
#include "PoissonSuffStat.hpp"
#include "Random.hpp"

/**
 * \brief A sufficient statistic for a collection of random variables, from a
 * mixture of a point mass at 0 and a Gamma distribution (see IIDBernoulliGamma)
 */

class BernoulliGammaSuffStat : public SuffStat {
  public:
    BernoulliGammaSuffStat() {}
    ~BernoulliGammaSuffStat() {}

    //! set suff stats to 0
    void Clear() {
        sum = 0;
        sumlog = 0;
        n0 = 0;
        n1 = 0;
    }

    //! add the contribution of one beta variate equal to 0 to this suffstat
    void AddNullSuffStat(int c = 1) { n0 += c; }

    //! add the contribution of one beta variate > 0 to this suffstat
    void AddPosSuffStat(double x, double logx, int c = 1) {
        sum += x;
        sumlog += logx;
        n1 += c;
    }

    //! (*this) += from
    void Add(const BernoulliGammaSuffStat &from) {
        sum += from.sum;
        sumlog += from.sumlog;
        n0 += from.n0;
        n1 += from.n1;
    }

    //! (*this) += from, operator version
    BernoulliGammaSuffStat &operator+=(const BernoulliGammaSuffStat &from) {
        Add(from);
        return *this;
    }

    //! return object size, when put into an MPI buffer
    unsigned int GetMPISize() const { return 4; }

    //! put object into MPI buffer
    void MPIPut(MPIBuffer &buffer) const { buffer << sum << sumlog << n0 << n1; }

    //! read object from MPI buffer
    void MPIGet(const MPIBuffer &buffer) { buffer >> sum >> sumlog >> n0 >> n1; }

    //! read a BernoulliGammaSuffStat from MPI buffer and add it to this
    void Add(const MPIBuffer &buffer) {
        double temp;
        buffer >> temp;
        sum += temp;
        buffer >> temp;
        sumlog += temp;

        int tmp;
        buffer >> tmp;
        n0 += tmp;
        buffer >> tmp;
        n1 += tmp;
    }

    //! return log prob of suff stats, as a function of parameters of the mixture
    //! and the Gamma distribution
    double GetLogProb(double pi, double alpha, double beta) const {
        double logbern = n0 * log(1 - pi) + n1 * log(pi);
        double loggamma =
            n1 * (alpha * log(beta) - Random::logGamma(alpha) + (alpha - 1) * sumlog - beta * sum);
        return logbern + loggamma;
    }

    //! get number of contributing variables that are == 0
    int GetN0() const { return n0; }
    //! get number of contributing variables that are > 0
    int GetN1() const { return n1; }
    //! get total sum of x_i, for all x_i's that are > 0
    double GetSum() const { return sum; }
    //! get total sum of log(x_i), for all x_i's that are > 0
    double GetSumLog() const { return sumlog; }

  private:
    double sum;
    double sumlog;
    int n0;
    int n1;
};

/**
 * \brief An array of IID random variables, from a mixture of a point mass at 0
 * (weight 1-pi) and a Gamma distribution (weight pi)
 */

class IIDBernoulliGamma : public SimpleArray<double> {
  public:
    //! constructor, parameterized by size, weight parameter pi of the mixture,
    //! and parameters of the Gamma (alpha and beta)
    IIDBernoulliGamma(int insize, double inpi, double inalpha, double inbeta)
        : SimpleArray<double>(insize), pi(inpi), alpha(inalpha), beta(inbeta) {
        Sample();
    }

    ~IIDBernoulliGamma() {}

    //! return weight of the Gamma component of the mixture
    double GetPi() const { return pi; }
    //! return alpha parameter of the Gamma distribution
    double GetShape() const { return alpha; }
    //! return beta parameter of the Gamma distribution
    double GetScale() const { return beta; }

    //! set weight of the Gamma component of the mixture to new value
    void SetPi(double inpi) { pi = inpi; }
    //! set alpha parameter of the Gamma distribution to new value
    void SetShape(double inalpha) { alpha = inalpha; }
    //! set beta parameter of the Gamma distribution to new value
    void SetScale(double inbeta) { beta = inbeta; }

    //! sample all entries from prior distribution
    void Sample() {
        for (int i = 0; i < GetSize(); i++) {
            if (Random::Uniform() < pi) {
                (*this)[i] = Random::Gamma(alpha, beta);
            } else {
                (*this)[i] = 0;
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
            ret += Random::logGammaDensity(GetVal(i), alpha, beta);
        }
        return ret;
    }

    //! add all entries to sufficient statistic
    void AddSuffStat(BernoulliGammaSuffStat &suffstat) {
        for (int i = 0; i < GetSize(); i++) {
            if (!GetVal(i)) {
                suffstat.AddNullSuffStat(1);
            } else {
                suffstat.AddPosSuffStat(GetVal(i), log(GetVal(i)), 1);
            }
        }
    }

    //! get number of positive entries
    int GetNpos() const {
        int n = 0;
        for (int i = 0; i < GetSize(); i++) {
            if (GetVal(i)) {
                n++;
            }
        }
        return n;
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
