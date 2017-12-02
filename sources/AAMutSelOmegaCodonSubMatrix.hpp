
#ifndef AAMUTSELCODONMATRIX_H
#define AAMUTSELCODONMATRIX_H

#include "CodonSubMatrix.hpp"
#include "Random.hpp"

using namespace std;

class AAMutSelOmegaCodonSubMatrix : public virtual NucCodonSubMatrix, public virtual OmegaCodonSubMatrix {
  public:

    AAMutSelOmegaCodonSubMatrix(const CodonStateSpace *instatespace, const SubMatrix *inNucMatrix, const vector<double>& inaa, double inomega,
                          bool innormalise = false)
        : SubMatrix(instatespace->GetNstate(), innormalise),
          CodonSubMatrix(instatespace, innormalise),
          NucCodonSubMatrix(instatespace, inNucMatrix, innormalise),
          OmegaCodonSubMatrix(instatespace,inomega,innormalise),
          aa(inaa) {}

    const vector<double>& GetAAFitnessProfile() const {return aa;}
    double GetFitness(int a) const {return aa[a] + 1e-8;}

  protected:

    void ComputeArray(int i) const /*override*/;
    void ComputeStationary() const /*override*/;

    // data members

    const vector<double>& aa;
};

#endif


