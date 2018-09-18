#include "GTRSubMatrix.hpp"
#include <iostream>
using namespace std;

// ---------------------------------------------------------------------------
//     GTRSubMatrix
// ---------------------------------------------------------------------------

GTRSubMatrix::GTRSubMatrix(
    int inNstate, const std::vector<double> &rr, const std::vector<double> &stat, bool innormalise)
    : SubMatrix(inNstate, innormalise), mRelativeRate(rr) {
    Nrr = Nstate * (Nstate - 1) / 2;
    CopyStationary(stat);
}

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
