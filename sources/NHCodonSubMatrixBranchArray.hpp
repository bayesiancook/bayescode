
#pragma once

#include "AAMutSelOmegaCodonSubMatrix.hpp"
#include "BranchArray.hpp"
#include "CodonSubMatrix.hpp"
#include "SubMatrix.hpp"

/**
 * \brief A BranchArray of MGOmegaCodonSubMatrix
 *
 * The constructor takes as arguments a single nucleotide matrix and a
 * BranchSelector<double> which returns the value of omega for each site. It
 * then constructs an array of MGOmegaCodonSubMatrix of same size as the omega
 * array, and such that the matrix on branch j takes the jth value of the omega
 * array as its dN/dS.
 *
 * This array derives from BranchArray and not SimpleBranchArray, because it is
 * implemented as a vector<MGOmegaCodonSubMatrix*>.
 */

class NHMGOmegaCodonSubMatrixBranchArray : public BranchArray<SubMatrix>,
                                         public BranchArray<MGOmegaCodonSubMatrix> {
  public:
    //! constructor parameterized by a codon state space, a single nucleotide
    //! matrix and an array (in fact, a BranchSelector) of omega's
    NHMGOmegaCodonSubMatrixBranchArray(const CodonStateSpace *incodonstatespace,
                                     const BranchSelector<SubMatrix> *innucmatrixarray,
                                     const BranchSelector<double> *inomegaarray)
        : codonstatespace(incodonstatespace),
          nucmatrixarray(innucmatrixarray),
          omegaarray(inomegaarray),
          matrixarray(inomegaarray->GetNbranch()) {
        Create();
    }

    ~NHMGOmegaCodonSubMatrixBranchArray() { Delete(); }

    //! return array size
    const Tree &GetTree() const { return omegaarray->GetTree(); }
    //! const access to matrix i
    int GetNbranch() const { return omegaarray->GetNbranch(); }
    const MGOmegaCodonSubMatrix &GetVal(int i) const { return *matrixarray[i]; }
    //! non-const access to matrix i
    MGOmegaCodonSubMatrix &operator[](int i) { return *matrixarray[i]; }

    //! update all matrices
    void UpdateCodonMatrices() {
        for (int i = 0; i < GetNbranch(); i++) {
            (*this)[i].SetOmega(omegaarray->GetVal(i));
            (*this)[i].CorruptMatrix();
        }
    }

  private:
    void Create() {
        for (int i = 0; i < GetNbranch(); i++) {
            matrixarray[i] =
                new MGOmegaCodonSubMatrix(codonstatespace, &nucmatrixarray->GetVal(i), omegaarray->GetVal(i));
        }
    }

    void Delete() {
        for (int i = 0; i < GetNbranch(); i++) {
            delete matrixarray[i];
        }
    }

    const CodonStateSpace *codonstatespace;
    const BranchSelector<SubMatrix> *nucmatrixarray;
    const BranchSelector<double> *omegaarray;
    vector<MGOmegaCodonSubMatrix *> matrixarray;
};

