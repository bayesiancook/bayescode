#include "GTRSubMatrix.hpp"
#include <iostream>
using namespace std;

// ---------------------------------------------------------------------------
//     GTRSubMatrix
// ---------------------------------------------------------------------------

GTRSubMatrix::GTRSubMatrix(const int inNstate, const std::vector<double> &rr,
    const std::vector<double> &stat, const bool innormalise)
    : SubMatrix(inNstate, innormalise), mRelativeRate(rr) {
    Nrr = Nstate * (Nstate - 1) / 2;
    CopyStationary(stat);
}

GTRSubMatrix::GTRSubMatrix(const GTRSubMatrix &other)
    : GTRSubMatrix(other.Nstate, other.mRelativeRate,
          std::vector<double>(
              other.mStationary.data(), other.mStationary.data() + other.mStationary.size()),
          other.normalise) {}

void GTRSubMatrix::CopyStationary(const std::vector<double> &instat) {
    for (int k = 0; k < Nstate; k++) { mStationary[k] = instat[k]; }
}

// ---------------------------------------------------------------------------
//     ComputeArray
// ---------------------------------------------------------------------------

void GTRSubMatrix::ComputeArray(int i) const {
    double total = 0;
    for (int j = 0; j < Nstate; j++) {
        if (i != j) {
            Q(i, j) = RelativeRate(i, j) * mStationary[j];
            total += Q(i, j);
        }
    }

    Q(i, i) = -total;
}
