
#include "AAMutSelSubMatrix.hpp"

void AAMutSelSubMatrix::ComputeStationary() const {
    for (int a = 0; a < Nstate; a++) { mStationary[a] = 0; }

    for (int i = 0; i < codonstatespace->GetNstate(); i++) {
        int a = codonstatespace->Translation(i);
        mStationary[a] += nucmatrix->Stationary(codonstatespace->GetCodonPosition(0, i)) *
                          nucmatrix->Stationary(codonstatespace->GetCodonPosition(1, i)) *
                          nucmatrix->Stationary(codonstatespace->GetCodonPosition(2, i)) *
                          GetFitness(a);
    }

    double total = 0;
    for (int a = 0; a < Nstate; a++) { total += mStationary[a]; }
    for (int a = 0; a < Nstate; a++) { mStationary[a] /= total; }
}

void AAMutSelSubMatrix::ComputeArray(int a) const {
    for (int b = 0; b < Nstate; b++) { Q(a, b) = 0; }

    double totweight = 0;
    for (int i = 0; i < codonstatespace->GetNstate(); i++) {
        if (codonstatespace->Translation(i) == a) {
            double weight = nucmatrix->Stationary(codonstatespace->GetCodonPosition(0, i)) *
                            nucmatrix->Stationary(codonstatespace->GetCodonPosition(1, i)) *
                            nucmatrix->Stationary(codonstatespace->GetCodonPosition(2, i)) *
                            GetFitness(a);
            totweight += weight;

            for (int j = 0; j < codonstatespace->GetNstate(); j++) {
                int b = codonstatespace->Translation(j);
                if (b != a) {
                    double mut = 1.0;
                    int ndiff = 0;
                    for (int pos = 0; pos < 3; pos++) {
                        int n1 = codonstatespace->GetCodonPosition(pos, i);
                        int n2 = codonstatespace->GetCodonPosition(pos, j);

                        if (n1 != n2) {
                            mut *= (*nucmatrix)(n1, n2);
                            if (ndiff) { mut *= xi; }
                            ndiff++;
                        }
                    }

                    if (mut > 0) {
                        double S = log(GetFitness(b)) - log(GetFitness(a));
                        double pfix = 1.0;

                        if ((fabs(S)) < 1e-30) {
                            pfix = 1 + S / 2;
                        } else if (S > 50) {
                            pfix = S;
                        } else if (S < -50) {
                            pfix = 0;
                        } else {
                            pfix = S / (1.0 - exp(-S));
                        }

                        Q(a, b) += weight * mut * pfix;
                    }
                }
            }
        }
    }

    double total = 0;
    for (int b = 0; b < Nstate; b++) {
        Q(a, b) /= totweight;
        total += Q(a, b);
    }

    Q(a, a) = -total;
}
