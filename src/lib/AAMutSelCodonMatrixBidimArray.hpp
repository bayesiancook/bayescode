#pragma once

#include "AAMutSelOmegaCodonSubMatrix.hpp"
#include "BidimArray.hpp"

/**
 * \brief A BidimArray of AAMutSelOmegaCodonSubMatrix
 *
 */

class AAMutSelCodonMatrixBidimArray : public BidimArray<SubMatrix>,
                                      public BidimArray<AAMutSelOmegaCodonSubMatrix> {
  public:
    //! constructor parameterized by an array of fitness profiles, a codon
    //! state space and a single nucleotide matrix.
    AAMutSelCodonMatrixBidimArray(const CodonStateSpace *incodonstatespace,
        const SubMatrix *innucmatrix, const Selector<std::vector<double>> *infitnessarray,
        const Selector<double> *indelta_omega_array, double inomega_shift)
        : codonstatespace(incodonstatespace),
          nucmatrix(innucmatrix),
          fitnessarray(infitnessarray),
          delta_omega_array(indelta_omega_array),
          omega_shift(inomega_shift),
          matrixarray(infitnessarray->GetSize(),
              std::vector<AAMutSelOmegaCodonSubMatrix *>(
                  indelta_omega_array->GetSize(), (AAMutSelOmegaCodonSubMatrix *)0)) {
        Create();
    }

    ~AAMutSelCodonMatrixBidimArray() { Delete(); }

    //! return the number of rows (number of fitness Ncat)
    virtual int GetNrow() const override { return fitnessarray->GetSize(); }
    //! return the number of columns (number of omega Ncat)
    virtual int GetNcol() const override { return delta_omega_array->GetSize(); }

    //! const access to the matrix for row (fitness cat) i and column (omega cat) j
    virtual const AAMutSelOmegaCodonSubMatrix &GetVal(int i, int j) const override {
        return *matrixarray[i][j];
    }
    //! non const access to the matrix for row (fitness cat) i and column (omega cat) j
    virtual AAMutSelOmegaCodonSubMatrix &operator()(int i, int j) override {
        return *matrixarray[i][j];
    }

    //! allocation and construction of all matrices
    void Create() {
        for (int i = 0; i < GetNrow(); i++) {
            for (int j = 0; j < GetNcol(); j++) {
                matrixarray[i][j] = new AAMutSelOmegaCodonSubMatrix(codonstatespace, nucmatrix,
                    fitnessarray->GetVal(i), omega_shift + delta_omega_array->GetVal(j), 1.0);
            }
        }
    }

    //! destruction of all matrices
    void Delete() {
        for (int i = 0; i < GetNrow(); i++) {
            for (int j = 0; j < GetNcol(); j++) { delete matrixarray[i][j]; }
        }
    }

    //! update all matrices
    void UpdateCodonMatrices() {
        for (int i = 0; i < GetNrow(); i++) {
            UpdateCodonMatrices(i);
        }
    }

    //! update all matrices for component i of the fitness mixture
    void UpdateCodonMatrices(int i) {
        for (int j = 0; j < GetNcol(); j++) {
            matrixarray[i][j]->SetOmega(omega_shift + delta_omega_array->GetVal(j));
            matrixarray[i][j]->CorruptMatrix();
        }
    }

    //! update only those matrices for which occupancy[i] != 0
    void UpdateCodonMatrices(const Selector<int> &occupancy) {
        for (int i = 0; i < GetNrow(); i++) {
            if (!occupancy.GetVal(i)) {
                UpdateCodonMatrices(i);
            }
        }
    }

  private:
    const CodonStateSpace *codonstatespace;
    const SubMatrix *nucmatrix;

    const Selector<std::vector<double>> *fitnessarray;
    const Selector<double> *delta_omega_array;
    double omega_shift;
    std::vector<std::vector<AAMutSelOmegaCodonSubMatrix *>> matrixarray;
};
