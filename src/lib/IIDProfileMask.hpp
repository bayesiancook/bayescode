#pragma once

#include <cmath>
#include "Array.hpp"

/**
 * \brief An array of IID 0/1 masks of a fixed dimension
 *
 * This class is used in AAMutSelSparseOmegaModel and DiffSelDoublySparseModel.
 * In those two cases, the masks are over the 20 amino-acids, and the array if
 * of size Nsite (i.e. the array implements site-specific masks over the
 * alignment). The masks are meant to define a low/high fitness distribution
 * over the 20 amino-acids.
 *
 * Each mask is made of dim iid Bernoulli(pi), conditional on at least one entry
 * of the mask being equal to 1.
 */

class IIDProfileMask : public SimpleArray<std::vector<int>> {
  public:
    //! constructor, parameterized by array size, mask dimension and Bernoulli
    //! probability parameter
    IIDProfileMask(int size, int indim, double pi)
        : SimpleArray(size, std::vector<int>(indim, 1)), dim(indim), pi(pi) {}

    //! return dimension of the masks
    int GetDim() const { return dim; }

    //! set probability parameter of the Bernoulli to a new value
    void SetPi(double inpi) { pi = inpi; }

    //! return total log probability over entire array
    double GetLogProb() const {
        double total = 0;
        for (int i = 0; i < GetSize(); i++) { total += GetLogProb(i); }
        return total;
    }

    //! return log probability for entry i
    double GetLogProb(int i) const {
        int naa = 0;
        const std::vector<int> &x = GetVal(i);
        for (int k = 0; k < GetDim(); k++) { naa += x[k]; }
        if (!naa) {
            std::cerr << "error in IIDProfileMask: all entries are null\n";
            exit(1);
        }
        // probability is conditional on at least one entry being 1
        return naa * log(pi) + (GetDim() - naa) * log(1.0 - pi) -
               log(1.0 - exp(GetDim() * log(1.0 - pi)));
    }

    //! return mean width (i.e. mean number of entries equal to 1) across the
    //! array
    double GetMeanWidth() const {
        double mean = 0;
        for (int i = 0; i < GetSize(); i++) {
            const std::vector<int> &x = GetVal(i);
            for (int k = 0; k < GetDim(); k++) { mean += x[k]; }
        }
        mean /= GetSize();
        return mean;
    }

  private:
    int dim;
    double pi;
};

class ProfileMask : public SimpleArray<std::vector<int>> {
  public:
    //! constructor, parameterized by array size, mask dimension and Bernoulli
    //! probability parameter
    ProfileMask(int size, const std::vector<double> &inpi)
        : SimpleArray(size, std::vector<int>(inpi.size(), 1)), dim(inpi.size()), pi(inpi) {}

    //! return dimension of the masks
    int GetDim() const { return dim; }

    //! return total log probability over entire array
    double GetLogProb() const {
        double total = 0;
        for (int i = 0; i < GetSize(); i++) { total += GetLogProb(i); }
        return total;
    }

    //! return log probability for entry i
    double GetLogProb(int i) const {
        int naa = 0;
        double ret = 0;
        double z = 1.0;
        const std::vector<int> &x = GetVal(i);
        for (int k = 0; k < GetDim(); k++) {
            if (x[k]) {
                ret += log(pi[k]);
            } else {
                ret += log(1 - pi[k]);
            }
            naa += x[k];
            z *= (1.0 - pi[k]);
        }
        if (!naa) {
            std::cerr << "error in IIDProfileMask: all entries are null\n";
            exit(1);
        }
        // probability is conditional on at least one entry being 1
        ret -= log(1.0 - z);
        return ret;
    }

    //! return mean width (i.e. mean number of entries equal to 1) across the
    //! array
    double GetMeanWidth() const {
        double mean = 0;
        for (int i = 0; i < GetSize(); i++) {
            const std::vector<int> &x = GetVal(i);
            for (int k = 0; k < GetDim(); k++) { mean += x[k]; }
        }
        mean /= GetSize();
        return mean;
    }

  private:
    int dim;
    const std::vector<double> &pi;
};