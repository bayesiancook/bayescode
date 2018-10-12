
#ifndef MULTIGAMMASUFFSTAT_H
#define MULTIGAMMASUFFSTAT_H

#include "IIDMultiGamma.hpp"
#include "SuffStat.hpp"

/**
 * \brief A sufficient statistic for a collection of MultiGamma variates, as a
 * function of the shape and center parameters
 *
 * see BidimIIDMultiGamma.
 */

class MultiGammaSuffStat : public SuffStat {
  public:
    //! constructor, parameterized by dimension of the multi-gamma
    MultiGammaSuffStat(int indim) : dim(indim), sum(indim), sumlog(indim), n(indim) {}
    ~MultiGammaSuffStat() {}

    //! return dimension of the multi-gamma
    int GetDim() const { return dim; }

    //! set suff stats to 0
    void Clear() {
        for (int k = 0; k < GetDim(); k++) {
            sum[k] = 0;
            sumlog[k] = 0;
            n[k] = 0;
        }
    }

    //! add the contribution of one gamma variate (x) to this suffstat
    void AddSuffStat(const std::vector<double> &x) {
        for (int k = 0; k < GetDim(); k++) {
            sum[k] += x[k];
            sumlog[k] += log(x[k]);
            n[k]++;
        }
    }

    //! add the contribution of one gamma variate (x) to this suffstat
    void AddSuffStat(const std::vector<double> &x, const std::vector<int> &t) {
        for (int k = 0; k < GetDim(); k++) {
            if (t[k]) {
                sum[k] += x[k];
                sumlog[k] += log(x[k]);
                n[k]++;
            }
        }
    }

    //! add the contribution of one gamma variate (x) to this suffstat
    void AddSuffStat(
        const std::vector<double> &x, const std::vector<int> &t1, const std::vector<int> &t2) {
        for (int k = 0; k < GetDim(); k++) {
            if (t1[k] && t2[k]) {
                sum[k] += x[k];
                sumlog[k] += log(x[k]);
                n[k]++;
            }
        }
    }

    //! (*this) += from
    void Add(const MultiGammaSuffStat &from) {
        for (int k = 0; k < GetDim(); k++) {
            sum[k] += from.GetSum(k);
            sumlog[k] += from.GetSumLog(k);
            n[k] += from.GetN(k);
        }
    }

    //! (*this) += from, operator version
    MultiGammaSuffStat &operator+=(const MultiGammaSuffStat &from) {
        Add(from);
        return *this;
    }

    //! get suff stats from an IIDGamma array
    void AddSuffStat(
        const BidimIIDMultiGamma &array, const BidimSelector<std::vector<int>> &toggle) {
        for (int i = 0; i < array.GetNrow(); i++) {
            if (!i) {
                for (int j = 0; j < array.GetNcol(); j++) { AddSuffStat(array.GetVal(i, j)); }
            } else {
                for (int j = 0; j < array.GetNcol(); j++) {
                    AddSuffStat(array.GetVal(i, j), toggle.GetVal(i - 1, j));
                }
            }
        }
    }

    //! get suff stats from an IIDGamma array, with a double system of masks
    void AddSuffStat(const BidimIIDMultiGamma &array, const Selector<std::vector<int>> &mask,
        const BidimSelector<std::vector<int>> &toggle) {
        for (int i = 0; i < array.GetNrow(); i++) {
            if (!i) {
                for (int j = 0; j < array.GetNcol(); j++) {
                    AddSuffStat(array.GetVal(i, j), mask.GetVal(j));
                }
            } else {
                for (int j = 0; j < array.GetNcol(); j++) {
                    AddSuffStat(array.GetVal(i, j), mask.GetVal(j), toggle.GetVal(i - 1, j));
                }
            }
        }
    }

    void AddSuffStat(const IIDMultiGamma &array) {
        for (int i = 0; i < array.GetSize(); i++) { AddSuffStat(array.GetVal(i)); }
    }

    void AddSuffStat(const IIDMultiGamma &array, const Selector<std::vector<int>> &mask) {
        for (int i = 0; i < array.GetSize(); i++) { AddSuffStat(array.GetVal(i), mask.GetVal(i)); }
    }

    //! return object size, when put into an MPI buffer
    unsigned int GetMPISize() const { return 3 * dim; }

    //! put object into MPI buffer
    void MPIPut(MPIBuffer &buffer) const {
        for (int k = 0; k < GetDim(); k++) { buffer << sum[k] << sumlog[k] << n[k]; }
    }

    //! read object from MPI buffer
    void MPIGet(const MPIBuffer &buffer) {
        for (int k = 0; k < GetDim(); k++) { buffer >> sum[k] >> sumlog[k] >> n[k]; }
    }

    //! read a MultiGammaSuffStat from MPI buffer and add it to this
    void Add(const MPIBuffer &buffer) {
        double temp;
        int tmp;
        for (int k = 0; k < GetDim(); k++) {
            buffer >> temp;
            sum[k] += temp;
            buffer >> temp;
            sumlog[k] += temp;
            buffer >> tmp;
            n[k] += tmp;
        }
    }

    //! return log prob, as a function of the given shape and scale parameters
    double GetLogProb(double shape, const std::vector<double> &center) const {
        double tot = 0;
        for (int k = 0; k < GetDim(); k++) {
            double alpha = shape * center[k];
            tot += n[k] * (-Random::logGamma(alpha)) + (alpha - 1) * sumlog[k] - sum[k];
            /*
            double scale = shape / center[k];
                    tot += n[k]*(shape*log(scale) - Random::logGamma(shape)) +
            (shape-1)*sumlog[k] - scale*sum[k];
            */
        }
        return tot;
    }

    //! return sum x_i's
    double GetSum(int k) const { return sum[k]; }
    //! return sum log x_i's
    double GetSumLog(int k) const { return sumlog[k]; }
    //! return N, total number of gamma variates contributing to the suff stat
    int GetN(int k) const { return n[k]; }

  private:
    int dim;
    std::vector<double> sum;
    std::vector<double> sumlog;
    std::vector<int> n;
};

#endif
