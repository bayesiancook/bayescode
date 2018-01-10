#ifndef GTRSUBMATRIX_H
#define GTRSUBMATRIX_H

#include "BiologicalSequences.hpp"  //FIXME only used for Naa (const int)
#include "SubMatrix.hpp"
#include "tinycompo.hpp"

/**
 * \brief A general time reversible substitution matrix
 *
 * Parameterized in terms of an array of relative rates (size Nstate*(Nstate-1)/2) and an array of equilibrium frequencies
 * (size Nstate)
 */

class GTRSubMatrix : public virtual SubMatrix, public tc::Component {
    const vector<double>* mRelativeRate;

  public:
    //! constructor parameterized by an array of relative rates (size Nstate*(Nstate-1)/2) and an array of equilibrium
    //! frequencies (size Nstate)
    GTRSubMatrix(int inNstate, bool innormalise = false) : SubMatrix(inNstate, innormalise), Nrr(Nstate * (Nstate - 1) / 2) {
        port("mRelativeRate", &GTRSubMatrix::mRelativeRate);
        port("CopyStationary", &GTRSubMatrix::CopyStationary);
    }
    ~GTRSubMatrix() = default;

    //! return number of relative rates
    int GetNRelativeRate() const { return Nrr; }
    //! converter returning relrate(i,j) == relrate(j,i), for any pair of states i!=j
    double RelativeRate(int i, int j) const { return (*mRelativeRate)[rrindex(i, j, GetNstate())]; }

    //! make a copy of the entries (not of the pointer) of the equilibrium frequency vector; should be done each time this
    //! vector has been modified
    void CopyStationary(const std::vector<double>* instat);

  protected:
    void ComputeArray(int i) const /*override*/;
    void ComputeStationary() const /*override*/ {}

    // const std::vector<double>& mRelativeRate;
    int Nrr;

  private:
    static int rrindex(int i, int j, int nstate) {
        return (i < j) ? (2 * nstate - i - 1) * i / 2 + j - i - 1 : (2 * nstate - j - 1) * j / 2 + i - j - 1;
    }
};

#endif
