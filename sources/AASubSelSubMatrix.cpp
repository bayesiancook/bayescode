#include "AASubSelSubMatrix.hpp"
#include <iostream>
using namespace std;

// ---------------------------------------------------------------------------
//     AASubSelSubMatrix
// ---------------------------------------------------------------------------

AASubSelSubMatrix::AASubSelSubMatrix(int inNstate, const std::vector<double> &rr,
                                     const std::vector<double> &stat, bool innormalise)
    : SubMatrix(inNstate, innormalise), mRelativeRate(rr) {
    Nrr = Nstate * (Nstate - 1) / 2;
    CopyStationary(stat);
}

void AASubSelSubMatrix::CopyStationary(const std::vector<double> &instat) {
    for (int k = 0; k < Nstate; k++) {
        mStationary[k] = instat[k];
    }
}

// ---------------------------------------------------------------------------
//     ComputeArray
// ---------------------------------------------------------------------------

void AASubSelSubMatrix::ComputeArray(int i) const {
    double total = 0;
    for (int j = 0; j < Nstate; j++) {
        if (i != j) {
            Q(i, j) = RelativeRate(i, j);

            double S = log(GetFitness(j)) - log(GetFitness(i));
            if ((fabs(S)) < 1e-30) {
                Q(i, j) *= 1 + S / 2;
            } else if (S > 50) {
                Q(i, j) *= S;
            } else if (S < -50) {
                Q(i, j) = 0;
            } else {
                Q(i, j) *= S / (1.0 - exp(-S));
            }

            total += Q(i, j);
        }
    }

    Q(i, i) = -total;
}
