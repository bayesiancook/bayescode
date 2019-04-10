#pragma once

#include "Array.hpp"
#include "BidimArray.hpp"
#include "global/Random.hpp"
#include "model/probnode_utils.hpp"

/**
 * \brief An array of site- and condition-specific fitness profiles -- doubly
 * sparse version
 *
 * this array is used in DiffSelDoublySparseModel; It implements two levels of
 * sparsity, as follows:
 * - F_0ia = G_0ia * m_ia + epsilon * (1-m_ia)
 * - F_kia = F_0ia^(1-d_kia) * G_kia^(d_kia), for k=1..K-1
 *
 * where
 * - G_kia is the input fitness matrix (for amino-acid a, site i and condition
 * k=0..K-1)
 * - m_ia  is an array of masks (for amino-acid a and site i)
 * - d_kia is an array of toggles (for amino-acid a, site i and condition
 * k=1..K-1)
 * - F_kia is the actual fitness of amino acid a at site i and under condition k
 * (for k=0..K-1)
 *
 * In words, and compared to DiffSelSparseFitnessArray, the baseline profiles
 * are first masked for some amino-acids (which are then considered as
 * low-fitness amino-acids across all conditions). Then, for non-baseline
 * conditions k=1..K, if d_kia == 0, the fitness of amino-acid a at site i under
 * condition k is just the (masked) baseline F_0ia; otherwise, it is a 'new'
 * fitness parameter, such as defined by G_kia.
 *
 * Note that, when Nlevel == 2, the relation between G, d and F unfolds over two
 * levels:
 * - F_0ia = G_0ia * m_ia + epsilon * (1-m_ia)
 * - F_1ia = F_0ia^(1-d_1ia) * G_1ia^(d_1ia);
 * - F_kia = F_1ia^(1-d_kia) * G_kia^(d_kia), for k=2..K-1
 */

class DiffSelDoublySparseFitnessArray : public SimpleBidimArray<std::vector<double>> {
  public:
    //! constructor, parameterized by input fitness array, toggle array and Nlevel
    DiffSelDoublySparseFitnessArray(const BidimSelector<std::vector<double>> &infitness,
        const Selector<std::vector<int>> &inmask,
        std::function<indicator_t(int, int, int)> intoggle, int inNlevel, const double &inepsilon)
        : SimpleBidimArray<std::vector<double>>(infitness.GetNrow(), infitness.GetNcol(),
              std::vector<double>(infitness.GetVal(0, 0).size(), 0)),
          fitness(infitness),
          mask(inmask),
          get_toggle(intoggle),
          epsilon(inepsilon),
          Nlevel(inNlevel) {
        Update();
    }

    //! return dimension of fitness profiles (should normally be 20)
    int GetDim() const { return GetVal(0, 0).size(); }

    //! full update of the array
    void Update() {
        for (int i = 0; i < GetNrow(); i++) {
            for (int j = 0; j < GetNcol(); j++) { Update(i, j); }
        }
    }

    //! update of column j (i.e. site j)
    void UpdateColumn(int j) {
        for (int i = 0; i < GetNrow(); i++) { Update(i, j); }
    }

    //! update of column j (i.e. site j) and condition i
    void Update(int i, int j) {
        std::vector<double> &x = (*this)(i, j);
        double total = 0;
        for (int k = 0; k < GetDim(); k++) {
            if (mask.GetVal(j)[k]) {
                int l = 0;
                if (i > 0) {
                    if (get_toggle(i, j, k)) {
                        l = i;
                    } else {
                        if ((Nlevel == 2) && (get_toggle(1, j, k))) { l = 1; }
                    }
                }
                x[k] = fitness.GetVal(l, j)[k];
            } else {
                x[k] = epsilon;
            }
            total += x[k];
        }
        for (int k = 0; k < GetDim(); k++) { x[k] /= total; }
    }

  protected:
    const BidimSelector<std::vector<double>> &fitness;
    const Selector<std::vector<int>> &mask;
    std::function<indicator_t(int, int, int)> get_toggle;
    // const BidimSelector<std::vector<int>> &toggle;
    const double &epsilon;
    int Nlevel;
};
