
#ifndef BERNBETA_H
#define BERNBETA_H

#include "Array.hpp"
#include "CodonSuffStat.hpp"
#include "MPIBuffer.hpp"
#include "PoissonSuffStat.hpp"
#include "Random.hpp"

class BernoulliBetaSuffStat : public SuffStat {
  public:
    BernoulliBetaSuffStat() {}
    ~BernoulliBetaSuffStat() {}

    void Clear() {
        sumlog0 = 0;
        sumlog1 = 0;
        n0 = 0;
        n1 = 0;
    }

    void AddNullSuffStat(int c = 1) { n0 += c; }

    void AddPosSuffStat(double log0, double log1, int c = 1) {
        sumlog0 += log0;
        sumlog1 += log1;
        n1 += c;
    }

    void Add(const BernoulliBetaSuffStat& from) {
        sumlog0 += from.GetSumLog0();
        sumlog1 += from.GetSumLog1();
        n0 += from.GetN0();
        n1 += from.GetN1();
    }

    BernoulliBetaSuffStat& operator+=(const BernoulliBetaSuffStat& from) {
        Add(from);
        return *this;
    }

    unsigned int GetMPISize() const { return 4; }

    void MPIPut(MPIBuffer& buffer) const { buffer << sumlog0 << sumlog1 << n0 << n1; }

    void MPIGet(const MPIBuffer& buffer) { buffer >> sumlog0 >> sumlog1 >> n0 >> n1; }

    void Add(const MPIBuffer& buffer) {
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

    double GetLogProb(double pi, double alpha, double beta) const {
        double logbern = n0 * log(1 - pi) + n1 * log(pi);
        double logbeta = n1 * (Random::logGamma(alpha + beta) - Random::logGamma(alpha) - Random::logGamma(beta)) +
                         (alpha - 1) * sumlog0 + (beta - 1) * sumlog1;
        return logbern + logbeta;
    }

    int GetN0() const { return n0; }
    int GetN1() const { return n1; }
    double GetSumLog0() const { return sumlog0; }
    double GetSumLog1() const { return sumlog1; }

  private:
    double sumlog0;
    double sumlog1;
    int n0;
    int n1;
};

class IIDBernoulliBeta : public SimpleArray<double> {
  public:
    IIDBernoulliBeta(int insize, double inpi, double inalpha, double inbeta)
        : SimpleArray<double>(insize), pi(inpi), alpha(inalpha), beta(inbeta) {
        Sample();
    }

    ~IIDBernoulliBeta() {}

    double GetPi() const { return pi; }
    double GetAlpha() const { return alpha; }
    double GetBeta() const { return beta; }

    void SetPi(double inpi) { pi = inpi; }
    void SetAlpha(double inalpha) { alpha = inalpha; }
    void SetBeta(double inbeta) { beta = inbeta; }

    void Sample() {
        for (int i = 0; i < GetSize(); i++) {
            if (Random::Uniform() < pi) {
                (*this)[i] = Random::BetaSample(alpha, beta);
            } else {
                (*this)[i] = 0;
            }
        }
    }

    int GetNullSet() const {
        int tot = 0;
        for (int i = 0; i < GetSize(); i++) {
            if (!GetVal(i)) {
                tot++;
            }
        }
        return tot;
    }

    double GetLogProb() const {
        double total = 0;
        for (int i = 0; i < GetSize(); i++) {
            total += GetLogProb(i);
        }
        return total;
    }

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

    void AddSuffStat(BernoulliBetaSuffStat& suffstat) {
        for (int i = 0; i < GetSize(); i++) {
            if (!GetVal(i)) {
                suffstat.AddNullSuffStat(1);
            } else {
                suffstat.AddPosSuffStat(log(GetVal(i)), log(1 - GetVal(i)), 1);
            }
        }
    }

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
