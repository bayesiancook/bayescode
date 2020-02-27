
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

    if (isinf(total))   {
        cerr << "in codon submatrix compute stationary: inf\n";
        exit(1);
    }

    if (isnan(total))   {
        cerr << "in codon submatrix compute stationary: nan\n";
        exit(1);
    }

    if (! total)    {
        cerr << "tot stat is 0\n";
        exit(1);
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

                Q(i, j) = (*NucMatrix)(a, b);

                double deltaS = 0;
                if (!Synonymous(i, j)) {
                    deltaS = log(GetFitness(GetCodonStateSpace()->Translation(j))) -
                             log(GetFitness(GetCodonStateSpace()->Translation(i)));
                }
                if ((fabs(deltaS)) < 1e-30) {
                    Q(i, j) *= 1 + deltaS / 2;
                } else if (deltaS > 50) {
                    Q(i, j) *= deltaS;
                } else if (deltaS < -50) {
                    Q(i, j) = 0;
                } else {
                    Q(i, j) *= deltaS / (1.0 - exp(-deltaS));
                }
                if (!Synonymous(i, j)) {
                    Q(i, j) *= GetOmega();
                }
            } else {
                Q(i, j) = 0;
            }
            total += Q(i, j);

            /*
            if (std::isinf(Q(i, j))) {
                cerr << "Q matrix infinite: " << Q(i, j) << '\n';
                exit(1);
            }

            if (Q(i, j) < 0) {
                cerr << "Q matrix negative: " << Q(i, j) << '\n';
                exit(1);
            }
            */
        }
    }

    if (isinf(total))   {
        cerr << "in codon submatrix compute Q: inf\n";
        exit(1);
    }

    if (isnan(total))   {
        cerr << "in codon submatrix compute Q: nan\n";
        exit(1);
    }

    if (total < 0) {
        cerr << "negative rate away\n";
        exit(1);
    }

    Q(i, i) = -total;
}

double AAMutSelOmegaCodonSubMatrix::GetPredictedDNDS() const {

	// UpdateMatrix();
    double totom = 0;
    double totweight = 0;
    for (int i = 0; i < Nstate; i++) {

		double weight = 0;
		double om = 0;

        for (int j = 0; j < Nstate; j++) {
            if (i != j) {
                if (!Synonymous(i, j)) {
                    int pos = GetDifferingPosition(i, j);
                    if ((pos != -1) && (pos != 3)) {
                        int a = GetCodonPosition(pos, i);
                        int b = GetCodonPosition(pos, j);

                        double nucrate = (*NucMatrix)(a, b);

                        double deltaS = log(GetFitness(GetCodonStateSpace()->Translation(j))) - log(GetFitness(GetCodonStateSpace()->Translation(i)));
                        double pfix = 1.0;
                        if ((fabs(deltaS)) < 1e-30) {
                            pfix = 1 + deltaS / 2;
                        } else if (deltaS > 50) {
                            pfix = deltaS;
                        } else if (deltaS < -50) {
                            pfix = 0;
                        } else {
                            pfix = deltaS / (1.0 - exp(-deltaS));
                        }

						om += nucrate*pfix;
						weight += nucrate;
                    }
                }
            }
        }

		totom += Stationary(i) * om;
		totweight += Stationary(i) * weight;
    }

    double om = totom / totweight;

    if (isinf(om))   {
        cerr << "in codon submatrix compute dN/dS: inf\n";
        cerr << totom << '\t' << totweight << '\n';
        exit(1);
    }

    if (isnan(om))   {
        cerr << "in codon submatrix compute dN/dS: nan\n";
        cerr << totom << '\t' << totweight << '\n';
        exit(1);
    }

    return om;
}

