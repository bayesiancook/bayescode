
#ifndef AAMUTSELCODONMATRIX_H
#define AAMUTSELCODONMATRIX_H

#include "CodonSubMatrix.hpp"
#include "Random.hpp"

using namespace std;

class AAMutSelOmegaCodonSubMatrix : public virtual NucCodonSubMatrix {
  public:

    AAMutSelOmegaCodonSubMatrix(const CodonStateSpace *instatespace, const SubMatrix *inNucMatrix, const vector<double>& inaa, double inomega,
                          bool innormalise = false)
        : SubMatrix(instatespace->GetNstate(), innormalise),
          CodonSubMatrix(instatespace, innormalise),
          NucCodonSubMatrix(instatespace, inNucMatrix, innormalise),
          aa(inaa),
          omega(inomega) {}

    double GetOmega() const { return omega + omegamin; }
    void SetOmega(double inomega) { omega = inomega; CorruptMatrix();}

    const vector<double>& GetAAFitnessProfile() const {return aa;}
    double GetFitness(int a) const {return aa[a];}

  protected:

    void ComputeArray(int i) const /*override*/;
    void ComputeStationary() const /*override*/;

    // data members

    const vector<double>& aa;
    double omega;
};

#endif


