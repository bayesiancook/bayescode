
#ifndef GTRMATRIXARRAY_H
#define GTRMATRIXARRAY_H

#include "Array.hpp"
#include "GTRSubMatrix.hpp"

/**
 * \brief An Array of GTRSubMatrix
 */

class GTRSubMatrixArray : public Array<SubMatrix>, public Array<GTRSubMatrix> {
  public:
    GTRSubMatrixArray(const vector<double> &inrelrate,
        const Selector<vector<double>> *inprofilearray, int innormalise)
        : relrate(inrelrate),
          profilearray(inprofilearray),
          nstate(inprofilearray->GetVal(0).size()),
          normalise(innormalise),
          matrixarray(inprofilearray->GetSize()) {
        Create();
    }

    ~GTRSubMatrixArray() { Delete(); }

    //! return array size
    int GetSize() const { return profilearray->GetSize(); }
    //! const access to matrix i
    const GTRSubMatrix &GetVal(int i) const { return *matrixarray[i]; }
    //! non-const access to matrix i
    GTRSubMatrix &operator[](int i) { return *matrixarray[i]; }

    //! update all matrices
    void UpdateMatrices() {
        for (int i = 0; i < GetSize(); i++) {
            (*this)[i].CopyStationary(profilearray->GetVal(i));
            (*this)[i].CorruptMatrix();
        }
    }

    void UpdateMatrix(int i) {
        (*this)[i].CopyStationary(profilearray->GetVal(i));
        (*this)[i].CorruptMatrix();
    }

  private:
    void Create() {
        for (int i = 0; i < GetSize(); i++) {
            matrixarray[i] = new GTRSubMatrix(nstate, relrate, profilearray->GetVal(i), normalise);
        }
    }

    void Delete() {
        for (int i = 0; i < GetSize(); i++) { delete matrixarray[i]; }
    }

    const vector<double> &relrate;
    const Selector<vector<double>> *profilearray;
    int nstate;
    int normalise;
    vector<GTRSubMatrix *> matrixarray;
};

#endif
