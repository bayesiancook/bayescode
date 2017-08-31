
#include "AAMutSelOmegaCodonSubMatrix.hpp"

void AAMutSelOmegaCodonSubMatrix::ComputeStationary() const {

    // compute stationary probabilities
    double total = 0;
    for (int i = 0; i < Nstate; i++) {
        mStationary[i] = NucMatrix->Stationary(GetCodonPosition(0, i)) *
                         NucMatrix->Stationary(GetCodonPosition(1, i)) *
                         NucMatrix->Stationary(GetCodonPosition(2, i)) *
                         GetFitness(GetCodonStateSpace()->Translation(i));
        total += mStationary[i];
    }

    // renormalize stationary probabilities
    // double min = 1;
    for (int i = 0; i < Nstate; i++) {
        mStationary[i] /= total;
    }
}


void AAMutSelOmegaCodonSubMatrix::ComputeArray(int i) const {

    double total = 0;
    for (int j = 0; j < Nstate; j++) {
        if (i != j) {
            int pos = GetDifferingPosition(i, j);
            if ((pos != -1) && (pos != 3)) {
                int a = GetCodonPosition(pos, i);
                int b = GetCodonPosition(pos, j);

                Q(i,j) = (*NucMatrix)(a, b);

                double deltaS = 0;
                if (!Synonymous(i, j)) {
                    deltaS = log(GetFitness(GetCodonStateSpace()->Translation(j))) - log(GetFitness(GetCodonStateSpace()->Translation(i)));
		}
                if ((fabs(deltaS)) < 1e-30) {
                    Q(i,j) *= 1 + deltaS / 2;
                } else if (deltaS > 50) {
                    Q(i,j) *= deltaS;
                } else if (deltaS < -50) {
                    Q(i,j) = 0;
                }
                if (deltaS != 0) {
                    Q(i,j) *= deltaS / (1.0 - exp(-deltaS));
                }
            } else {
                Q(i,j) = 0;
            }
            total += Q(i,j);

            if (isinf(Q(i,j))) {
                cerr << "Q matrix infinite: " << Q(i,j) << '\n';
                exit(1);
            }

            if (Q(i,j) < 0) {
                cerr << "Q matrix negative: " << Q(i,j) << '\n';
                exit(1);
            }
        }
    }

    Q(i,i) = -total;

    if (total < 0) {
        cerr << "negative rate away\n";
        exit(1);
    }
}
