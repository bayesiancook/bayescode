#pragma once

#include "SuffStat.hpp"

class MeanPoissonSuffStat : public SuffStat {
  public:
    MeanPoissonSuffStat() {}
    MeanPoissonSuffStat(double incount, double inbeta) : count(incount), beta(inbeta) {}
    MeanPoissonSuffStat(const MeanPoissonSuffStat& from) : count(from.count), beta(from.beta) {}

    ~MeanPoissonSuffStat() {}

    //! set count and beta to 0
    void Clear() { count = beta = 0; }

    void IncrementCount() { count++; }

    void AddCount(double in) { count += in; }

    void AddBeta(double in) { beta += in; }

    void AddSuffStat(double incount, double inbeta) {
        count += incount;
        beta += inbeta;
    }

    void Add(const MeanPoissonSuffStat &from) {
        count += from.GetCount();
        beta += from.GetBeta();
    }

    MeanPoissonSuffStat &operator+=(const MeanPoissonSuffStat &from) {
        Add(from);
        return *this;
    }

    //! return size when put into an MPI buffer
    unsigned int GetMPISize() const { return 2; }

    //! put current value of count and beta into an MPI buffer
    void MPIPut(MPIBuffer &buffer) const { buffer << beta << count; }

    //! get value from MPI buffer
    void MPIGet(const MPIBuffer &buffer) { buffer >> beta >> count; }

    //! get a PoissonSuffStat from MPI buffer and then add it to this object
    void Add(const MPIBuffer &buffer) {
        double temp;
        buffer >> temp;
        beta += temp;

        double tmp;
        buffer >> tmp;
        count += tmp;
    }

    double GetCount() const { return count; }

    double GetBeta() const { return beta; }

    //! return the log probability as a function of the rate: essentially
    //! count*log(rate) - beta*rate
    double GetLogProb(double rate) const { return count * log(rate) - beta * rate; }

    //! return the log of the marginal probability when the rate is from a gamma
    //! distribution
    double GetMarginalLogProb(double shape, double scale) const {
        return shape * log(scale) - Random::logGamma(shape) - (shape + count) * log(scale + beta) +
               Random::logGamma(shape + count);
    }

    double count;
    double beta;
};

