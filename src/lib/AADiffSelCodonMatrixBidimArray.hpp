
#ifndef AADIFFSELCODONMATRIXARRAY_H
#define AADIFFSELCODONMATRIXARRAY_H

#include "AAMutSelOmegaCodonSubMatrix.hpp"
#include "BidimArray.hpp"

/**
 * \brief A BidimArray of AAMutSelOmegaCodonSubMatrix (used in DiffSelModel and
 * DiffSelSparseModel)
 *
 * The constructor takes as arguments a single nucleotide matrix and a bidim
 * array of amino-acid fitness profiles (across sites and conditions). It then
 * constructs an array of AAMutSelOmegaCodonSubMatrix of same size as the array
 * of fitness profiles.
 */

class AADiffSelCodonMatrixBidimArray : public BidimArray<SubMatrix>,
                                       public BidimArray<AAMutSelOmegaCodonSubMatrix> {
  public:
    // two possible ways to do it:
    // AAMutSelNeCodonMatrixBidimArray(const Selector<vector<double> >&
    // infitnessarray, const Selector<double>& inNe, const CodonStateSpace&
    // incodonstatespace, const SubMatrix& innucmatrix);
    // AAMutSelNeCodonMatrixBidimArray(const Selector<vector<double> >&
    // infitnessarray, const vector<double>& inNe, const CodonStateSpace&
    // incodonstatespace, const SubMatrix& innucmatrix);

    //! constructor parameterized by a bidim array of fitness profiles, a codon
    //! state space and a single nucleotide matrix.
    AADiffSelCodonMatrixBidimArray(const BidimSelector<std::vector<double>> &infitnessarray,
        const CodonStateSpace &incodonstatespace, const SubMatrix &innucmatrix)
        : fitnessarray(infitnessarray),
          codonstatespace(incodonstatespace),
          nucmatrix(innucmatrix),
          matrixarray(infitnessarray.GetNrow(),
              std::vector<AAMutSelOmegaCodonSubMatrix *>(
                  infitnessarray.GetNcol(), (AAMutSelOmegaCodonSubMatrix *)0)) {
        Create();
    }

    ~AADiffSelCodonMatrixBidimArray() { Delete(); }

    //! return the number of rows (number of conditions)
    virtual int GetNrow() const override { return fitnessarray.GetNrow(); }
    //! return the number of columns (number of sites)
    virtual int GetNcol() const override { return fitnessarray.GetNcol(); }
    //! const access to the matrix for row (condition) i and column (site) j
    virtual const AAMutSelOmegaCodonSubMatrix &GetVal(int i, int j) const override {
        return *matrixarray[i][j];
    }
    //! non const access to the matrix for row (condition) i and column (site) j
    virtual AAMutSelOmegaCodonSubMatrix &operator()(int i, int j) override {
        return *matrixarray[i][j];
    }

    //! allocation and construction of all matrices
    void Create() {
        for (int i = 0; i < GetNrow(); i++) {
            for (int j = 0; j < GetNcol(); j++) {
                matrixarray[i][j] = new AAMutSelOmegaCodonSubMatrix(
                    &codonstatespace, &nucmatrix, fitnessarray.GetVal(i, j), 1.0, 1.0);
            }
        }
    }

    //! destruction of all matrices
    void Delete() {
        for (int i = 0; i < GetNrow(); i++) {
            for (int j = 0; j < GetNcol(); j++) { delete matrixarray[i][j]; }
        }
    }

    //! signal corruption of the parameters (matrices should recompute themselves)
    void Corrupt() {
        for (int i = 0; i < GetNrow(); i++) {
            for (int j = 0; j < GetNcol(); j++) { matrixarray[i][j]->CorruptMatrix(); }
        }
    }

    //! signal corruption for column (site) j
    void CorruptColumn(int j) {
        for (int i = 0; i < GetNrow(); i++) { matrixarray[i][j]->CorruptMatrix(); }
    }

    //! signal corruption for column (site) j, and only for those rows
    //! (conditions) that are flagged
    void CorruptColumn(int j, const std::vector<int> &flag) {
        for (int i = 0; i < GetNrow(); i++) {
            if (flag[i]) { matrixarray[i][j]->CorruptMatrix(); }
        }
    }

  private:
    const BidimSelector<std::vector<double>> &fitnessarray;
    const CodonStateSpace &codonstatespace;
    const SubMatrix &nucmatrix;
    std::vector<std::vector<AAMutSelOmegaCodonSubMatrix *>> matrixarray;
};

#endif
