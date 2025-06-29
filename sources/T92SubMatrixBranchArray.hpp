
#pragma once

#include "T92SubMatrix.hpp"
#include "BranchArray.hpp"

class T92SubMatrixBranchArray : public BranchArray<SubMatrix>,
                                         public BranchArray<T92SubMatrix> {
  public:
    //! constructor parameterized by a codon state space, a single nucleotide
    //! matrix and an array (in fact, a BranchSelector) of omega's
    T92SubMatrixBranchArray(const BranchSelector<double> *intstv, const BranchSelector<double>* ingc, bool innormalize)
        : tstv(intstv), gc(ingc), normalize(innormalize),
          matrixarray(intstv->GetNbranch()) {
        Create();
    }

    ~T92SubMatrixBranchArray() { Delete(); }

    //! return array size
    const Tree &GetTree() const { return tstv->GetTree(); }
    //! const access to matrix i
    int GetNbranch() const { return tstv->GetNbranch(); }
    const T92SubMatrix &GetVal(int i) const { return *matrixarray[i]; }
    //! non-const access to matrix i
    T92SubMatrix &operator[](int i) { return *matrixarray[i]; }

    //! update all matrices
    void UpdateMatrices() {
        for (int i = 0; i < GetNbranch(); i++) {
            (*this)[i].SetKappa(tstv->GetVal(i));
            (*this)[i].SetGC(gc->GetVal(i));
            (*this)[i].CorruptMatrix();
        }
    }

  private:
    void Create() {
        for (int i = 0; i < GetNbranch(); i++) {
            matrixarray[i] =
                new T92SubMatrix(tstv->GetVal(i), gc->GetVal(i), normalize);
        }
    }

    void Delete() {
        for (int i = 0; i < GetNbranch(); i++) {
            delete matrixarray[i];
        }
    }

    const BranchSelector<double> *tstv;
    const BranchSelector<double> *gc;
    bool normalize;
    vector<T92SubMatrix *> matrixarray;
};

