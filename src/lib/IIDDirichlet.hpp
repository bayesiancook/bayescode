
#ifndef IIDDIR_H
#define IIDDIR_H

#include "Array.hpp"
#include "MPIBuffer.hpp"
#include "Random.hpp"
#include "SuffStat.hpp"

/**
 * \brief A sufficient statistic for a collection of Dirichlet random variables,
 * as a function of the parameters of the Dirichlet distribution
 *
 * Consider a Dirichlet distribution, over vectors of dimension K, and with
 * parameter alpha = (alpha_k)_k=1..K. Consider N Dirichlet random variables
 * (x_i)_i=1..N ~ Dirichlet(alpha). Then, up to an additive constant (not
 * depending on alpha):
 *
 * sum_i log p(x_i | alpha) = N (log Gamma(|alpha|) - sum_k log Gamma(alpha_k))
 * + sum_k S_k ^ (alpha_k-1)
 *
 * where S_k = sum_i log(x_ik).
 */

class DirichletSuffStat : public SuffStat {
  public:
    DirichletSuffStat(int indim) : sumlog(indim, 0), n(0) {}
    ~DirichletSuffStat() {}

    //! set suff stats to 0
    void Clear() {
        for (unsigned int i = 0; i < sumlog.size(); i++) { sumlog[i] = 0; }
        n = 0;
    }

    //! get dimension of underlying Dirichlet distribution
    int GetDim() { return (int)sumlog.size(); }

    //! add the contribution of one variate (x) to this suffstat
    void AddSuffStat(const std::vector<double> &pi) {
        for (unsigned int i = 0; i < sumlog.size(); i++) {
            if (pi[i] <= 0) {
                std::cerr << "error: negative pi in DirichletSuffStat: " << pi[i] << '\n';
                for (unsigned int j = 0; j < sumlog.size(); j++) { std::cerr << pi[j] << '\t'; }
                std::cerr << '\n';
                exit(1);
            }
            sumlog[i] += log(pi[i]);
        }
        n++;
    }

    //! add the contribution of one variate (x) to this suffstat
    void AddSuffStat(const double *insumlog, int d) {
        for (unsigned int i = 0; i < sumlog.size(); i++) { sumlog[i] += insumlog[i]; }
        n += d;
    }

    //! (*this) += from
    void Add(const DirichletSuffStat &from) {
        for (unsigned int i = 0; i < sumlog.size(); i++) { sumlog[i] += from.GetSumLog(i); }
        n += from.GetN();
    }

    //! (*this) += from, operator version
    DirichletSuffStat &operator+=(const DirichletSuffStat &from) {
        Add(from);
        return *this;
    }

    //! return S_k
    double GetSumLog(int i) const { return sumlog[i]; }

    //! return N (number of Dirichlet variables contributing to the sufficient
    //! statistic)
    int GetN() const { return n; }

    //! return logprob, as a function of center (renormalized alpha) and
    //! concentration (sum of the alpha vector)
    double GetLogProb(const std::vector<double> &center, double concentration) const {
        double tot = n * Random::logGamma(concentration);
        for (unsigned int i = 0; i < sumlog.size(); i++) {
            tot += -n * Random::logGamma(concentration * center[i]) +
                   (concentration * center[i] - 1) * sumlog[i];
        }
        return tot;
    }

    //! return object size, when put into an MPI buffer
    unsigned int GetMPISize() const { return sumlog.size() + 1; }

    //! put object into MPI buffer
    void MPIPut(MPIBuffer &buffer) const {
        for (unsigned int i = 0; i < sumlog.size(); i++) { buffer << sumlog[i]; }
        buffer << n;
    }

    //! read object from MPI buffer
    void MPIGet(const MPIBuffer &buffer) {
        for (unsigned int i = 0; i < sumlog.size(); i++) { buffer >> sumlog[i]; }
        buffer >> n;
    }

    //! read a DirichletSuffStat from MPI buffer and add it to this
    void Add(const MPIBuffer &buffer) {
        double tmp;
        for (unsigned int i = 0; i < sumlog.size(); i++) {
            buffer >> tmp;
            sumlog[i] += tmp;
        }
        int temp;
        buffer >> temp;
        n += temp;
    }

  private:
    std::vector<double> sumlog;
    int n;
};

/**
 * \brief A SimpleArray of DirichletSuffStat
 */

class DirichletSuffStatArray : public SimpleArray<DirichletSuffStat> {
  public:
    //! constructor, parameterized by array size and dimension of the Dirichlet
    //! distribution
    DirichletSuffStatArray(int insize, int indim)
        : SimpleArray<DirichletSuffStat>(insize, DirichletSuffStat(indim)), dim(indim) {}
    ~DirichletSuffStatArray() {}

    //! return dimension of the Dirichlet distribution
    int GetDim() const { return dim; }

    //! set suff stat to 0
    void Clear() {
        for (int i = 0; i < GetSize(); i++) { (*this)[i].Clear(); }
    }

    //! element-wise addition of array given as argument to this array
    void Add(const DirichletSuffStatArray &from) {
        for (int i = 0; i < GetSize(); i++) { (*this)[i].Add(from.GetVal(i)); }
    }

    //! element-wise addition, operator version
    DirichletSuffStatArray &operator+=(const DirichletSuffStatArray &from) {
        Add(from);
        return *this;
    }

    //! return object size, when put into an MPI buffer
    unsigned int GetMPISize() const { return GetSize() * GetVal(0).GetMPISize(); }

    //! put object into MPI buffer
    void MPIPut(MPIBuffer &buffer) const {
        for (int i = 0; i < GetSize(); i++) { buffer << GetVal(i); }
    }

    //! read object from MPI buffer
    void MPIGet(const MPIBuffer &buffer) {
        for (int i = 0; i < GetSize(); i++) { buffer >> (*this)[i]; }
    }

    //! read a DirichletSuffStatArray from MPI buffer and add it to this
    void Add(const MPIBuffer &buffer) {
        for (int i = 0; i < GetSize(); i++) { (*this)[i] += buffer; }
    }

  private:
    int dim;
};

/**
 * \brief An array of IID Dirichlet random variables
 */

class IIDDirichlet : public SimpleArray<std::vector<double>> {
  public:
    //! constructor, parameterized by array size and parameters of the Dirichlet
    //! distribution (center and concentration)
    IIDDirichlet(int insize, const std::vector<double> &incenter, double inconcentration)
        : SimpleArray<std::vector<double>>(insize),
          center(incenter),
          concentration(inconcentration) {
        for (int i = 0; i < GetSize(); i++) { (*this)[i].assign(center.size(), 0); }
        Sample();
    }

    ~IIDDirichlet() {}

    //! set center of the Dirichlet distribution
    void SetCenter(const std::vector<double> &incenter) { center = incenter; }

    //! set concentration of the Dirichlet distribution
    void SetConcentration(double inconcentration) { concentration = inconcentration; }

    //! set all entries equal to uniform vector
    void SetUniform() {
        int dim = GetDim();
        for (int i = 0; i < GetSize(); i++) {
            for (int k = 0; k < dim; k++) { (*this)[i][k] = 1.0 / dim; }
        }
    }

    //! get dimension of Dirichlet distribution
    int GetDim() const { return center.size(); }

    //! sample from prior
    void Sample() {
        for (int i = 0; i < GetSize(); i++) {
            Random::DirichletSample((*this)[i], center, concentration);
        }
    }

    //! get log probability density of all entries, given current values of center
    //! and concentration parameters
    double GetLogProb() const {
        double total = 0;
        for (int i = 0; i < GetSize(); i++) { total += GetLogProb(i); }
        return total;
    }

    //! get log probability of entry i
    double GetLogProb(int i) const {
        return Random::logDirichletDensity(GetVal(i), center, concentration);
    }

    //! add this array to sufficient statistic given as argument
    void AddSuffStat(DirichletSuffStat &suffstat) const {
        for (int i = 0; i < GetSize(); i++) { suffstat.AddSuffStat(GetVal(i)); }
    }

    //! add this array to sufficient statistic given as argument -- only entries
    //! with non zero occupancy
    void AddSuffStat(DirichletSuffStat &suffstat, const Selector<int> &occupancy) const {
        for (int i = 0; i < GetSize(); i++) {
            if (occupancy.GetVal(i)) { suffstat.AddSuffStat(GetVal(i)); }
        }
    }

    //! resample from prior those entries for which occupancy[i] == 0
    void PriorResample(const Selector<int> &occupancy) {
        for (int i = 0; i < GetSize(); i++) {
            if (!occupancy.GetVal(i)) {
                Random::DirichletSample((*this)[i], center, concentration);
            }
        }
    }

    //! get mean entropy over all elements of the array
    double GetMeanEntropy() const {
        double mean = 0;
        for (int i = 0; i < GetSize(); i++) { mean += Random::GetEntropy(GetVal(i)); }
        mean /= GetSize();
        return mean;
    }

    //! get mean of component k of all elements of the array
    double GetMean(int k) const {
        double m1 = 0;
        for (int i = 0; i < GetSize(); i++) { m1 += GetVal(i)[k]; }
        m1 /= GetSize();
        return m1;
    }

    //! get variance of component k of all elements of the array
    double GetVar(int k) const {
        double m1 = 0;
        double m2 = 0;
        for (int i = 0; i < GetSize(); i++) {
            m1 += GetVal(i)[k];
            m2 += GetVal(i)[k] * GetVal(i)[k];
        }
        m1 /= GetSize();
        m2 /= GetSize();
        m2 -= m1 * m1;
        return m2;
    }

  protected:
    std::vector<double> center;
    double concentration;
};

/**
 * \brief An array of Dirichlet random variables, each with its own center and
 * concentration parameters
 */

class MultiDirichlet : public SimpleArray<std::vector<double>> {
  public:
    //! constructor, parameterized by arrays of center and concentration
    //! parameters (both of same size, which will also be the size of this array)
    MultiDirichlet(const Selector<std::vector<double>> *incenterarray,
        const Selector<double> *inconcentrationarray)
        : SimpleArray<std::vector<double>>(incenterarray->GetSize()),
          dim(incenterarray->GetVal(0).size()),
          centerarray(incenterarray),
          concentrationarray(inconcentrationarray),
          weightarray(0) {
        if (centerarray->GetSize() != concentrationarray->GetSize()) {
            std::cerr << "error in multi dirichlet: center and concentration arrays "
                         "should have same size\n";
            exit(1);
        }

        for (int i = 0; i < GetSize(); i++) { (*this)[i].assign(dim, 0); }
        Sample();
    }

    //! constructor, parameterized by arrays of center and concentration
    //! parameters (both of same size, which will also be the size of this array)
    MultiDirichlet(const Selector<std::vector<double>> *inweightarray)
        : SimpleArray<std::vector<double>>(inweightarray->GetSize()),
          dim(inweightarray->GetVal(0).size()),
          centerarray(0),
          concentrationarray(0),
          weightarray(inweightarray) {
        for (int i = 0; i < GetSize(); i++) { (*this)[i].assign(dim, 0); }
        Sample();
    }

    ~MultiDirichlet() {}

    //! return dimension of Dirichlet distribution
    int GetDim() const { return dim; }

    //! sample all entries from prior
    void Sample() {
        if (weightarray) {
            for (int i = 0; i < GetSize(); i++) {
                Random::DirichletSample((*this)[i], weightarray->GetVal(i));
            }
        } else {
            for (int i = 0; i < GetSize(); i++) {
                Random::DirichletSample(
                    (*this)[i], centerarray->GetVal(i), concentrationarray->GetVal(i));
            }
        }
    }

    //! get total log prob (sum over all array)
    double GetLogProb() const {
        double total = 0;
        for (int i = 0; i < GetSize(); i++) { total += GetLogProb(i); }
        return total;
    }

    //! get log prob for entry i
    double GetLogProb(int i) const {
        double ret = 0;
        if (weightarray) {
            ret = Random::logDirichletDensity(GetVal(i), weightarray->GetVal(i));
        } else {
            ret = Random::logDirichletDensity(
                GetVal(i), centerarray->GetVal(i), concentrationarray->GetVal(i));
        }
        return ret;
    }

    //! \brief add entries of this array to an array DirichletSuffStat, based on
    //! specified allocation vector
    void AddSuffStat(Array<DirichletSuffStat> &suffstatarray, const Selector<int> &alloc) {
        for (int i = 0; i < GetSize(); i++) {
            suffstatarray[alloc.GetVal(i)].AddSuffStat(GetVal(i));
        }
    }

    //! resample entries for which occupancy[i] == 0
    void PriorResample(const Selector<int> &occupancy) {
        for (int i = 0; i < GetSize(); i++) {
            if (!occupancy.GetVal(i)) {
                if (weightarray) {
                    Random::DirichletSample((*this)[i], weightarray->GetVal(i));
                } else {
                    Random::DirichletSample(
                        (*this)[i], centerarray->GetVal(i), concentrationarray->GetVal(i));
                }
            }
        }
    }

    //! get mean entropy over all elements of the array
    double GetMeanEntropy() const {
        double mean = 0;
        for (int i = 0; i < GetSize(); i++) { mean += Random::GetEntropy(GetVal(i)); }
        mean /= GetSize();
        return mean;
    }

    //! get mean of component k of all elements of the array
    double GetMean(int k) const {
        double m1 = 0;
        for (int i = 0; i < GetSize(); i++) { m1 += GetVal(i)[k]; }
        m1 /= GetSize();
        return m1;
    }

    //! get variance of component k of all elements of the array
    double GetVar(int k) const {
        double m1 = 0;
        double m2 = 0;
        for (int i = 0; i < GetSize(); i++) {
            m1 += GetVal(i)[k];
            m2 += GetVal(i)[k] * GetVal(i)[k];
        }
        m1 /= GetSize();
        m2 /= GetSize();
        m2 -= m1 * m1;
        return m2;
    }

    void Flatten(){
        for (int cat = 0; cat < GetSize(); cat++) {
            std::fill((*this)[cat].begin(), (*this)[cat].end(), 1.0 / (*this)[cat].size());
        }
    }

  protected:
    int dim;
    const Selector<std::vector<double>> *centerarray;
    const Selector<double> *concentrationarray;
    const Selector<std::vector<double>> *weightarray;
};

#endif
