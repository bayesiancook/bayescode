#ifndef AASUBSELSUBMATRIX_H
#define AASUBSELSUBMATRIX_H

#include "BiologicalSequences.hpp"  //FIXME only used for Naa (const int)
#include "SubMatrix.hpp"

class AASubSelSubMatrix : public virtual SubMatrix {
  public:
    //! constructor parameterized by an array of relative rates (size Nstate*(Nstate-1)/2) and an array of equilibrium frequencies (size Nstate)
    AASubSelSubMatrix(int inNstate, const std::vector<double>& rr, const std::vector<double>& stat, bool innormalise = false);
    ~AASubSelSubMatrix() override = default;

    //! return number of relative rates
    int GetNRelativeRate() const { return Nrr; }
    //! converter returning relrate(i,j) == relrate(j,i), for any pair of states i!=j
    double RelativeRate(int i, int j) const { return mRelativeRate[rrindex(i, j)]; }

    //! make a copy of the entries of the equilibrium frequency vector; should be done each time this vector has been modified
    void CopyStationary(const std::vector<double>& instat);

  protected:
    void ComputeArray(int i) const override;
    void ComputeStationary() const override {}

    double GetFitness(int i) const {return mStationary(i);}

    const std::vector<double>&  mRelativeRate;
    int Nrr;

  private:
    int rrindex(int i, int j) const {
        return (i < j) ? (2 * Nstate - i - 1) * i / 2 + j - i - 1
                       : (2 * Nstate - j - 1) * j / 2 + i - j - 1;
    }

};

#endif
