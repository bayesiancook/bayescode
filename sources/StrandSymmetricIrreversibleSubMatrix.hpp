
#ifndef StrandSymmetricIrreversibleSUBMATRIX_H
#define StrandSymmetricIrreversibleSUBMATRIX_H

#include "BiologicalSequences.hpp"  //FIXME only used for Naa (const int)
#include "SubMatrix.hpp"

/* general strand symmetric - irreversible
AC - TG
AG - TC
AT - TA
CA - GT
CG - GC
CT - GA
*/

class StrandSymmetricIrreversibleSubMatrix : public virtual SubMatrix {
  public:

    StrandSymmetricIrreversibleSubMatrix(const vector<double>& inrates, bool innormalise = false)
        : SubMatrix(4, innormalise), rates(inrates) {}
    ~StrandSymmetricIrreversibleSubMatrix() override = default;

  protected:
    void ComputeArray(int i) const override;
    void ComputeStationary() const override;

    const vector<double>& rates;
};

/*
template<class Rates>
auto make_strand_symmetric_matrix(Rates rates, bool normalize = false)  {
    return new StrandSymmetricIrreversibleSubMatrix<Rates>(rates, normalize);
}
*/

#endif
