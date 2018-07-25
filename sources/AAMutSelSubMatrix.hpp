#ifndef AASUBSELSUBMATRIX_H
#define AASUBSELSUBMATRIX_H

#include "BiologicalSequences.hpp"  //FIXME only used for Naa (const int)
#include "CodonStateSpace.hpp"
#include "SubMatrix.hpp"

class AAMutSelSubMatrix : public virtual SubMatrix {
  public:
    AAMutSelSubMatrix(const CodonStateSpace *incodonstatespace, const SubMatrix *innucmatrix,
                      double inxi, const vector<double> &inaa, double inNe,
                      bool innormalise = false)
        : SubMatrix(Naa, innormalise),
          codonstatespace(incodonstatespace),
          nucmatrix(innucmatrix),
          xi(inxi),
          aa(inaa),
          Ne(inNe) {}
    ~AAMutSelSubMatrix() override = default;

    //! const access (by reference) to amino-acid fitness vector
    const vector<double> &GetAAFitnessProfile() const { return aa; }

    //! \brief access by copy to fitness of a given amino-acid
    //!
    //! Note: to avoid numerical errors, this function returns aa[a] + 1e-8.
    double GetFitness(int a) const { return Ne * aa[a] + 1e-8; }

    void SetNe(double inNe) {
        Ne = inNe;
        CorruptMatrix();
    }

    void SetXi(double inxi) {
        xi = inxi;
        CorruptMatrix();
    }

  protected:
    void ComputeArray(int i) const override;
    void ComputeStationary() const override;

    const CodonStateSpace *codonstatespace;
    const SubMatrix *nucmatrix;
    double xi;
    const vector<double> &aa;
    double Ne;
};

#endif
