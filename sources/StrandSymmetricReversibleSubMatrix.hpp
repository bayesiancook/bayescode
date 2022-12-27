#ifndef StrandSymmetricReversibleSUBMATRIX_H
#define StrandSymmetricReversibleSUBMATRIX_H

#include "BiologicalSequences.hpp"  //FIXME only used for Naa (const int)
#include "SubMatrix.hpp"

/* general strand symmetric - reversible
rho_AC rho_AG rho_AT rho_CG
*/

class StrandSymmetricReversibleSubMatrix : public virtual SubMatrix {
  public:

    StrandSymmetricReversibleSubMatrix(const vector<double>& inrho, double ingc, bool innormalise = false)
        : SubMatrix(4, innormalise), rho(inrho), gc(ingc) {}
    ~StrandSymmetricReversibleSubMatrix() override = default;

    void SetGC(double ingc) {
        gc = ingc;
        CorruptMatrix();
    }

  protected:
    void ComputeArray(int i) const override;
    void ComputeStationary() const override;

    const vector<double>& rho;
    double gc;
};

#endif
