#ifndef T92SUBMATRIX_H
#define T92SUBMATRIX_H

#include "BiologicalSequences.hpp"  //FIXME only used for Naa (const int)
#include "SubMatrix.hpp"

class T92SubMatrix : public virtual SubMatrix {
  public:
    //! constructor parameterized by an array of relative rates (size
    //! Nstate*(Nstate-1)/2) and an array of equilibrium frequencies (size Nstate)
    T92SubMatrix(double inkappa, double ingc, bool innormalise = false)
        : SubMatrix(4, innormalise), kappa(inkappa), gc(ingc) {}
    ~T92SubMatrix() override = default;

    void SetKappa(double inkappa) {
        kappa = inkappa;
        CorruptMatrix();
    }

    void SetGC(double ingc) {
        gc = ingc;
        CorruptMatrix();
    }

  protected:
    void ComputeArray(int i) const override;
    void ComputeStationary() const override;

    double kappa;
    double gc;
};

#endif
