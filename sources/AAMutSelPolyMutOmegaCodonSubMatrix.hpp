#pragma once

#include "PolyMutNucCodonSubMatrix.hpp"
#include "Random.hpp"

using namespace std;

/**
 * \brief A mutation-selection codon substitution process.
 *
 * This codon substitution process describes the evolution of a coding position,
 * under a constant fitness landscape over the 20 amino-acids. The model is
 * parameterized by a nucleotide matrix, specifying the mutation process, and a
 * vector of 20 scaled fitness parameters (summing to 1), for the 20
 * amino-acids. The process also takes a real parameter, omega, which acts as a
 * multiplier in front of all non-synonymous substitutions. The standard
 * mutation-selection process is obtained by setting omega=1. Letting omega be
 * different from 1 was explored in Rodrigue and Lartillot, 2017, MBE (detecting
 * deviations from expected non-syn rate under the standard mut-sel model).
 */

class AAMutSelPolyMutOmegaCodonSubMatrix : public virtual PolyMutNucCodonSubMatrix,
                                    public virtual OmegaCodonSubMatrix {
  public:
    //! constructor, parameterized by a codon state space (genetic code), a
    //! nucleotide mutation matrix, a 20-vector of amino-acid fitnesss, and a
    //! positive real parameter omega (=1 in the standard model).
    AAMutSelPolyMutOmegaCodonSubMatrix(const CodonStateSpace *instatespace, const SubMatrix *inNucMatrix,
                                const vector<double> &inaa, double inu, double inomega, double inNe,
                                bool innormalise = false)
        : SubMatrix(instatespace->GetNstate(), innormalise),
          CodonSubMatrix(instatespace, innormalise),
          PolyMutNucCodonSubMatrix(instatespace, inNucMatrix, inu, innormalise),
          OmegaCodonSubMatrix(instatespace, inomega, innormalise),
          aa(inaa),
          Ne(inNe) {}

    //! const access (by reference) to amino-acid fitness vector
    const vector<double> &GetAAFitnessProfile() const { return aa; }

    //! \brief access by copy to fitness of a given amino-acid
    //!
    //! Note: to avoid numerical errors, this function returns aa[a] + 1e-8.
    double GetFitness(int a) const { return Ne * aa[a] + 1e-8; }

  protected:
    void SetNe(double inNe) {
        Ne = inNe;
        CorruptMatrix();
    }

    void ComputeStationary() const override {
        // compute stationary probabilities
        double total = 0;
        for (int i = 0; i < Nstate; i++) {
            mStationary[i] = NucMatrix->Stationary(GetCodonPosition(0, i)) *
                             NucMatrix->Stationary(GetCodonPosition(1, i)) *
                             NucMatrix->Stationary(GetCodonPosition(2, i)) *
                             GetFitness(GetCodonStateSpace()->Translation(i));
            total += mStationary[i];
        }

        if (std::isinf(total))   {
            cerr << "in codon submatrix compute stationary: inf\n";
            exit(1);
        }

        if (std::isnan(total))   {
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

    void ComputeArray(int i) const override {
        double total = 0;
        for (int j = 0; j < Nstate; j++) {
            if (i != j) {
                int npos = 0;
                double nucrate = 1.0;
                for (int pos=0; pos<3; pos++) {
                    if (GetCodonStateSpace()->GetCodonPosition(pos, i) != GetCodonStateSpace()->GetCodonPosition(pos, j))   {
                        if (npos)   {
                            nucrate *= u;
                        }
                        npos++;

                        int a = GetCodonPosition(pos, i);
                        int b = GetCodonPosition(pos, j);
                        nucrate *= (*NucMatrix)(a, b);
                    }
                }
                Q(i,j) = nucrate;

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
                total += Q(i, j);
            }
        }

        if (std::isinf(total))   {
            cerr << "in codon submatrix compute Q: inf\n";
            exit(1);
        }

        if (std::isnan(total))   {
            cerr << "in codon submatrix compute Q: nan\n";
            exit(1);
        }

        if (total < 0) {
            cerr << "negative rate away\n";
            cerr << total << '\n';
            exit(1);
        }

        Q(i, i) = -total;
    }

  public:

    double GetPredictedDNDS() const {

        // UpdateMatrix();
        double totom = 0;
        double totweight = 0;
        for (int i = 0; i < Nstate; i++) {

            double weight = 0;
            double om = 0;

            for (int j = 0; j < Nstate; j++) {
                if (i != j) {
                    if (!Synonymous(i, j)) {

                        int npos = 0;
                        double nucrate = 1.0;
                        for (int pos=0; pos<3; pos++) {
                            if (GetCodonStateSpace()->GetCodonPosition(pos, i) != GetCodonStateSpace()->GetCodonPosition(pos, j))   {
                                if (npos)   {
                                    nucrate *= u;
                                }
                                npos++;

                                int a = GetCodonPosition(pos, i);
                                int b = GetCodonPosition(pos, j);
                                nucrate *= (*NucMatrix)(a, b);
                            }
                        }

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

            totom += Stationary(i) * om;
            totweight += Stationary(i) * weight;
        }

        double om = totom / totweight;

        if (std::isinf(om))   {
            cerr << "in codon submatrix compute dN/dS: inf\n";
            cerr << totom << '\t' << totweight << '\n';
            exit(1);
        }

        if (std::isnan(om))   {
            cerr << "in codon submatrix compute dN/dS: nan\n";
            cerr << totom << '\t' << totweight << '\n';
            exit(1);
        }

        return om;
    }


  private:

    // data members

    const vector<double> &aa;
    double Ne;
};

