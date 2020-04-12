#pragma once

#include "CodonSubMatrix.hpp"

class PolyMutNucCodonSubMatrix : public virtual CodonSubMatrix {
  public:
    PolyMutNucCodonSubMatrix(const CodonStateSpace *instatespace, const SubMatrix *inNucMatrix, double inu,
                      bool innormalise)
        : SubMatrix(instatespace->GetNstate(), innormalise),
          CodonSubMatrix(instatespace, innormalise), u(inu) {
        SetNucMatrix(inNucMatrix);
    }

    const SubMatrix *GetNucMatrix() const { return NucMatrix; }

    double GetU() const {
        return u;
    }

    void SetU(double inu)   {
        u = inu;
        CorruptMatrix();
    }

  protected:

    void SetNucMatrix(const SubMatrix *inmatrix) {
        NucMatrix = inmatrix;
        if (NucMatrix->GetNstate() != Nnuc) {
            std::cerr << "error in CodonSubMatrix: underyling mutation process "
                         "should be a 4x4 matrix\n";
            throw;
        }
    }

    const SubMatrix *NucMatrix;
    double u;
};
