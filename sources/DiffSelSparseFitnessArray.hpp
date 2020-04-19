
#ifndef DIFFSELSPARSEFIT_H
#define DIFFSELSPARSEFIT_H

#include "Array.hpp"
#include "BidimArray.hpp"

/**
 * \brief An array of site- and condition-specific fitness profiles -- sparse
 * version
 *
 * this array is used in DiffSelSparseModel; It implements the following
 * deterministic relation:
 * - F_0ia = G_0ia;
 * - F_kia = F_0ia^(1-d_kia) * G_kia^(d_kia), for k=1..K-1
 *
 * where
 * - G_kia is the input fitness matrix (for amino-acid a, site i and condition
 * k=0..K-1)
 * - d_kia is an array of toggles (for amino-acid a, site i and condition
 * k=1..K-1)
 * - F_kia is the actual fitness of amino acid a at site i and under condition k
 * (for k=0..K-1)
 *
 * In words, if d_kia == 0, then F_kia, the fitness of amino-acid a at site i
 * under condition k, is just the baseline G_0ia; otherwise, it is a 'new'
 * fitness parameter, such as defined by G_kia. Note that, when Nlevel == 2, the
 * relation between G, d and F unfolds over two levels:
 * - F_0ia = G_0ia;
 * - F_1ia = G_0ia^(1-d_1ia) * G_1ia^(d_1ia);
 * - F_kia = F_1ia^(1-d_kia) * G_kia^(d_kia), for k=2..K-1
 */

class DiffSelSparseFitnessArray : public SimpleBidimArray<vector<double>> {
  public:
    //! constructor, parameterized by input fitness array, toggle array and Nlevel
    DiffSelSparseFitnessArray(const BidimSelector<vector<double>> &infitness,
                              const BidimSelector<vector<int>> &intoggle, int inNlevel)
        : SimpleBidimArray<vector<double>>(infitness.GetNrow(), infitness.GetNcol(),
                                           vector<double>(infitness.GetVal(0, 0).size(), 0)),
          fitness(infitness),
          toggle(intoggle),
          Nlevel(inNlevel) {
        Update();
    }

    //! return dimension of fitness profiles (should normally be 20)
    int GetDim() const { return GetVal(0, 0).size(); }

    //! full update of the array
    void Update() {
        for (int i = 0; i < GetNrow(); i++) {
            for (int j = 0; j < GetNcol(); j++) {
                Update(i, j);
            }
        }
    }

    //! update of column j (i.e. site j)
    void UpdateColumn(int j) {
        for (int i = 0; i < GetNrow(); i++) {
            Update(i, j);
        }
    }

    //! update of column j (i.e. site j) and condition i
    void Update(int i, int j) {
        vector<double> &x = (*this)(i, j);
        double total = 0;
        for (int k = 0; k < GetDim(); k++) {
            int l = 0;
            if (i > 0) {
                if (toggle.GetVal(i - 1, j)[k]) {
                    l = i;
                } else {
                    if ((Nlevel == 2) && (toggle.GetVal(0, j)[k])) {
                        l = 1;
                    }
                }
            }
            x[k] = fitness.GetVal(l, j)[k];
            total += x[k];
        }
        for (int k = 0; k < GetDim(); k++) {
            x[k] /= total;
        }
    }

  protected:
    const BidimSelector<vector<double>> &fitness;
    const BidimSelector<vector<int>> &toggle;
    int Nlevel;
};

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
 * low-fitness amino-acids across all conditions). Then, dor non-baseline
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

class DiffSelDoublySparseFitnessArray : public SimpleBidimArray<vector<double>> {
  public:
    //! constructor, parameterized by input fitness array, toggle array and Nlevel
    DiffSelDoublySparseFitnessArray(const BidimSelector<vector<double>> &infitness,
                                    const Selector<vector<int>> &inmask,
                                    const BidimSelector<vector<int>> &intoggle, int inNlevel,
                                    double inepsilon)
        : SimpleBidimArray<vector<double>>(infitness.GetNrow(), infitness.GetNcol(),
                                           vector<double>(infitness.GetVal(0, 0).size(), 0)),
          fitness(infitness),
          mask(inmask),
          toggle(intoggle),
          Nlevel(inNlevel),
          epsilon(inepsilon) {
        Update();
    }

    //! return dimension of fitness profiles (should normally be 20)
    int GetDim() const { return GetVal(0, 0).size(); }

    //! notify new value for epsilon parameter
    void SetEpsilon(double ineps) { epsilon = ineps; }

    //! full update of the array
    void Update() {
        for (int i = 0; i < GetNrow(); i++) {
            for (int j = 0; j < GetNcol(); j++) {
                Update(i, j);
            }
        }
    }

    //! update of column j (i.e. site j)
    void UpdateColumn(int j) {
        for (int i = 0; i < GetNrow(); i++) {
            Update(i, j);
        }
    }

    //! update of column j (i.e. site j) and condition i
    void Update(int i, int j) {
        vector<double> &x = (*this)(i, j);
        double total = 0;
        for (int k = 0; k < GetDim(); k++) {
            if (mask.GetVal(j)[k]) {
                int l = 0;
                if (i > 0) {
                    if (toggle.GetVal(i - 1, j)[k]) {
                        l = i;
                    } else {
                        if ((Nlevel == 2) && (toggle.GetVal(0, j)[k])) {
                            l = 1;
                        }
                    }
                }
                x[k] = fitness.GetVal(l, j)[k];
            } else {
                x[k] = epsilon;
            }
            total += x[k];
        }
        for (int k = 0; k < GetDim(); k++) {
            x[k] /= total;
        }
    }

  protected:
    const BidimSelector<vector<double>> &fitness;
    const Selector<vector<int>> &mask;
    const BidimSelector<vector<int>> &toggle;
    int Nlevel;
    double epsilon;
};

/**
 * \brief A sparse array of site-specific fitness profiles
 *
 * This array is used in AAMutSelSparseOmegaModel. It implements sparse fitness
 * profiles across sites by masking some amino-acids, according to an array of
 * binary variables (masks), at each site. Specifically, the fitness at site i
 * for amino-acid a is given by:
 * - F_ia = G_ia * m_ia + epsilon * (1-m_ia)
 *
 * where
 * - G_ia is the input fitness array (for amino-acid a and site i)
 * - m_ia  is an array of masks (for amino-acid a and site i)
 * - F_ia is the actual fitness of amino acid a at site i
 */

class MutSelSparseFitnessArray : public SimpleArray<vector<double>> {
  public:
    //! constructor, parameterized by input fitness array, mask array and epsilon
    //! (background fitness of low-fitness amino-acids)
    MutSelSparseFitnessArray(const Selector<vector<double>> &infitness,
                             const Selector<vector<int>> &inmask, double inepsilon)
        : SimpleArray<vector<double>>(infitness.GetSize(),
                                      vector<double>(infitness.GetVal(0).size(), 0)),
          fitness(infitness),
          mask(inmask),
          epsilon(inepsilon) {
        Update();
    }

    //! returns dimension of fitness profiles (should normally be 20)
    int GetDim() const { return GetVal(0).size(); }

    //! notify new value for epsilon parameter
    void SetEpsilon(double ineps) { epsilon = ineps; }

    //! full update of the array
    void Update() {
        for (int i = 0; i < GetSize(); i++) {
            Update(i);
        }
    }

    //! update site i
    void Update(int i) {
        vector<double> &x = (*this)[i];
        double total1 = 0;
        int n1 = 0;
        for (int k = 0; k < GetDim(); k++) {
            if (mask.GetVal(i)[k]) {
                n1++;
                total1 += fitness.GetVal(i)[k];
            }
	    }
        double total2 = 0;
        for (int k = 0; k < GetDim(); k++) {
            if (mask.GetVal(i)[k]) {
                x[k] = fitness.GetVal(i)[k] / total1 * n1;
            } else {
                x[k] = epsilon;
            }
            total2 += x[k];
        }
        for (int k = 0; k < GetDim(); k++) {
            x[k] /= total2;
        }
    }

    double GetMeanEntropy() const {
        double mean = 0;
        for (int i = 0; i < GetSize(); i++) {
            mean += Random::GetEntropy(GetVal(i));
        }
        mean /= GetSize();
        return mean;
    }

  protected:
    const Selector<vector<double>> &fitness;
    const Selector<vector<int>> &mask;
    double epsilon;
};

#endif
