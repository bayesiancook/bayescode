
#ifndef BETA_H
#define BETA_H

#include "Array.hpp"
#include "CodonSuffStat.hpp"
#include "MPIBuffer.hpp"
#include "PoissonSuffStat.hpp"
#include "Random.hpp"

class BetaSuffStat : public SuffStat {
  public:
    BetaSuffStat() {}
    ~BetaSuffStat() {}

    void Clear() {
        sumlog0 = 0;
        sumlog1 = 0;
        n = 0;
    }

    void AddSuffStat(double log0, double log1, int c = 1) {
        sumlog0 += log0;
        sumlog1 += log1;
        n += c;
    }

    void Add(const BetaSuffStat& from) {
        sumlog0 += from.GetSumLog0();
        sumlog1 += from.GetSumLog1();
        n += from.GetN();
    }

    BetaSuffStat& operator+=(const BetaSuffStat& from) {
        Add(from);
        return *this;
    }

    unsigned int GetMPISize() const { return 3; }

    void MPIPut(MPIBuffer& buffer) const { buffer << sumlog0 << sumlog1 << n; }

    void MPIGet(const MPIBuffer& buffer) { buffer >> sumlog0 >> sumlog1 >> n; }

    void Add(const MPIBuffer& buffer) {
        double temp;
        buffer >> temp;
        sumlog0 += temp;
        buffer >> temp;
        sumlog1 += temp;

        int tmp;
        buffer >> tmp;
        n += tmp;
    }

    double GetMeanInvConcLogProb(double mean, double invconc) const {
        double alpha = mean / invconc;
        double beta = (1 - mean) / invconc;
        return GetLogProb(alpha, beta);
    }

    double GetLogProb(double alpha, double beta) const {
        return n * (Random::logGamma(alpha + beta) - Random::logGamma(alpha) - Random::logGamma(beta)) +
               (alpha - 1) * sumlog0 + (beta - 1) * sumlog1;
    }

    int GetN() const { return n; }
    double GetSumLog0() const { return sumlog0; }
    double GetSumLog1() const { return sumlog1; }

  private:
    double sumlog0;
    double sumlog1;
    int n;
};

class IIDBeta : public SimpleArray<double> {
  public:
    IIDBeta(int insize, double inalpha, double inbeta) : SimpleArray<double>(insize), alpha(inalpha), beta(inbeta) {
        Sample();
    }

    ~IIDBeta() {}

    double GetAlpha() const { return alpha; }
    double GetBeta() const { return beta; }

    double GetMean() const { return alpha / (alpha + beta); }
    double GetInvConcentration() const { return 1.0 / (alpha + beta); }

    void SetMeanInvConc(double mean, double invconc) {
        double alpha = mean / invconc;
        double beta = (1 - mean) / invconc;
        SetAlpha(alpha);
        SetBeta(beta);
    }

    void SetAlpha(double inalpha) { alpha = inalpha; }
    void SetBeta(double inbeta) { beta = inbeta; }

    void Sample() {
        for (int i = 0; i < GetSize(); i++) {
            (*this)[i] = Random::BetaSample(alpha, beta);
        }
    }

    double GetLogProb() const {
        double total = 0;
        for (int i = 0; i < GetSize(); i++) {
            total += GetLogProb(i);
        }
        return total;
    }

    double GetLogProb(int i) const { return Random::logBetaDensity(GetVal(i), alpha, beta); }

    void AddSuffStat(BetaSuffStat& suffstat) {
        for (int i = 0; i < GetSize(); i++) {
            suffstat.AddSuffStat(log(GetVal(i)), log(1 - GetVal(i)), 1);
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
    double alpha;
    double beta;
};

#endif
