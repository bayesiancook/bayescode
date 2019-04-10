#ifndef MULTINOMALLOC_H
#define MULTINOMALLOC_H

#include "Array.hpp"
#include "OccupancySuffStat.hpp"
#include "global/Random.hpp"

/**
 * \brief A multinomial allocation vector
 *
 * MultinomialAllocationVector specifies the probability distribution of the
 * multinomial allocation of N items to the K components of the finite mixture.
 * The Distribution is parameterized by a vector of K weight.
 */

class MultinomialAllocationVector : public SimpleArray<int> {
  public:
    //! Constructor (parameterized by number of items (N) and mixture weights
    //! (vector of dimension K)
    MultinomialAllocationVector(int insize, const std::vector<double> &inweight)
        : SimpleArray<int>(insize), weight(inweight) {
        SampleAlloc();
    }

    ~MultinomialAllocationVector() {}

    //! sample allocations based on current weights (i.e. from the prior)
    void SampleAlloc() {
        for (int i = 0; i < GetSize(); i++) {
            (*this)[i] = Random::DrawFromDiscreteDistribution(weight);
        }
    }

    //! resample allocation of item i based on an externally given vector of
    //! posterior probabilities
    void GibbsResample(int i, const std::vector<double> &postprob) {
        (*this)[i] = Random::DrawFromDiscreteDistribution(postprob);
    }

    //! resample all allocations based on an externally given array of
    //! item-specific posterior probabilities
    void GibbsResample(const std::vector<std::vector<double>> &postprobarray) {
        for (int i = 0; i < GetSize(); i++) {
            (*this)[i] = Random::DrawFromDiscreteDistribution(postprobarray[i]);
        }
    }

    //! apply permutation to allocation vector
    void Permute(const Selector<int> &permut) override {
        if (permut.GetSize() != int(weight.size())) {
            std::cerr << "error in MultinomialAllocationVector::Permute: non matching "
                         "array size\n";
            exit(1);
        }
        std::vector<int> invpermut(permut.GetSize(), 0);
        for (int k = 0; k < permut.GetSize(); k++) { invpermut[permut.GetVal(k)] = k; }
        for (int i = 0; i < GetSize(); i++) { (*this)[i] = invpermut[(*this)[i]]; }
    }

    //! apply single-exchange permutation to allocation vector
    void SwapComponents(int cat1, int cat2) {
        for (int i = 0; i < GetSize(); i++) {
            if ((*this)[i] == cat1) {
                (*this)[i] = cat2;
            } else if ((*this)[i] == cat2) {
                (*this)[i] = cat1;
            }
        }
    }

  private:
    const std::vector<double> &weight;
};

#endif
