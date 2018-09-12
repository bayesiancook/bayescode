
#ifndef AAMUTSELMATRIXARRAY_H
#define AAMUTSELMATRIXARRAY_H

#include "AAMutSelSubMatrix.hpp"
#include "Array.hpp"

/**
 * \brief An Array of AAMutSelSubMatrix
 */

class AAMutSelSubMatrixArray : public Array<SubMatrix>, public Array<AAMutSelSubMatrix> {
  public:
    AAMutSelSubMatrixArray(const CodonStateSpace *incodonstatespace, const SubMatrix *innucmatrix,
        double inxi, const Selector<vector<double>> *inprofilearray, int innormalise)
        : codonstatespace(incodonstatespace),
          nucmatrix(innucmatrix),
          xi(inxi),
          profilearray(inprofilearray),
          normalise(innormalise),
          matrixarray(inprofilearray->GetSize()) {
        Create();
    }

    ~AAMutSelSubMatrixArray() { Delete(); }

    //! return array size
    int GetSize() const { return profilearray->GetSize(); }
    //! const access to matrix i
    const AAMutSelSubMatrix &GetVal(int i) const { return *matrixarray[i]; }
    //! non-const access to matrix i
    AAMutSelSubMatrix &operator[](int i) { return *matrixarray[i]; }

    void SetXi(double inxi) { xi = inxi; }

    //! update all matrices
    void UpdateMatrices() {
        for (int i = 0; i < GetSize(); i++) {
            (*this)[i].SetXi(xi);
            (*this)[i].CorruptMatrix();
        }
    }

    void UpdateMatrix(int i) { (*this)[i].CorruptMatrix(); }

  private:
    void Create() {
        for (int i = 0; i < GetSize(); i++) {
            matrixarray[i] = new AAMutSelSubMatrix(
                codonstatespace, nucmatrix, xi, profilearray->GetVal(i), 1.0, normalise);
        }
    }

    void Delete() {
        for (int i = 0; i < GetSize(); i++) { delete matrixarray[i]; }
    }

    const CodonStateSpace *codonstatespace;
    const SubMatrix *nucmatrix;
    double xi;
    const Selector<vector<double>> *profilearray;
    int normalise;
    vector<AAMutSelSubMatrix *> matrixarray;
};

#endif
