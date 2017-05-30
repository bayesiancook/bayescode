#ifndef MSCODONSUBMATRIX_H
#define MSCODONSUBMATRIX_H

#include <iostream>
#include "CodonSubMatrix.hpp"

// square root
class MGSRFitnessCodonUsageSubMatrix : public MGCodonSubMatrix {
  public:
    MGSRFitnessCodonUsageSubMatrix(const CodonStateSpace *instatespace, const SubMatrix *inNucMatrix,
                                   const double *infitness, const double *incodonusageselection,
                                   bool innormalise = false)
        : SubMatrix(instatespace->GetNstate(), innormalise),
          CodonSubMatrix(instatespace, innormalise),
          MGCodonSubMatrix(instatespace, inNucMatrix, innormalise),
          fitness(infitness),
          codonusageselection(incodonusageselection) {}

    double GetFitness(int aastate) const { return fitness[aastate] + 1e-6; }

    double GetCodonUsageSelection(int codonstate) const { return codonusageselection[codonstate]; }

  protected:
    // look at how ComputeArray and ComputeStationary are implemented in
    // CodonSubMatrix.cpp

    void ComputeArray(int i) const override;
    void ComputeStationary() const override;
    // data members
    const double *fitness;
    const double *codonusageselection;
};

// mutation selection
class MGMSFitnessCodonUsageSubMatrix : public MGCodonSubMatrix {
  public:
    MGMSFitnessCodonUsageSubMatrix(const CodonStateSpace *instatespace, const SubMatrix *inNucMatrix,
                                   const double *infitness, const double *incodonusageselection,
                                   bool innormalise = false)
        : SubMatrix(instatespace->GetNstate(), innormalise),
          CodonSubMatrix(instatespace, innormalise),
          MGCodonSubMatrix(instatespace, inNucMatrix, innormalise),
          fitness(infitness),
          codonusageselection(incodonusageselection) {}

    double GetFitness(int aastate) const { return fitness[aastate] + 1e-6; }

    double GetCodonUsageSelection(int codonstate) const { return codonusageselection[codonstate]; }

  protected:
    // look at how ComputeArray and ComputeStationary are implemented in
    // CodonSubMatrix.cpp

    void ComputeArray(int i) const override;
    void ComputeStationary() const override;

/*
    void ToStream(std::ostream & /*os*/) override {  // FIXME unused parameter
        std::cerr << "nucmatrix : \n";
        NucMatrix->CheckReversibility();
        std::cerr << '\n';
        NucMatrix->CorruptMatrix();
        NucMatrix->UpdateMatrix();
        NucMatrix->ToStream(std::cerr);
        std::cerr << '\n';
        NucMatrix->CheckReversibility();
        std::cerr << '\n';

        for (int i = 0; i < Naa; i++) {
            std::cerr << fitness[i] << '\t';
        }
        std::cerr << '\n';
        std::cerr << '\n';
        double total = 0;
        for (int i = 0; i < GetNstate(); i++) {
            total += NucMatrix->Stationary(GetCodonPosition(0, i)) *
                     NucMatrix->Stationary(GetCodonPosition(1, i)) *
                     NucMatrix->Stationary(GetCodonPosition(2, i)) *
                     GetFitness(GetCodonStateSpace()->Translation(i)) * GetCodonUsageSelection(i);
        }
        for (int i = 0; i < GetNstate(); i++) {
            std::cerr << i << '\t' << NucMatrix->Stationary(GetCodonPosition(0, i));
            std::cerr << '\t' << NucMatrix->Stationary(GetCodonPosition(1, i));
            std::cerr << '\t' << NucMatrix->Stationary(GetCodonPosition(2, i));
            std::cerr << '\t' << GetFitness(GetCodonStateSpace()->Translation(i));
            std::cerr << '\t' << GetCodonUsageSelection(i);
            std::cerr << '\t'
                      << NucMatrix->Stationary(GetCodonPosition(0, i)) *
                             NucMatrix->Stationary(GetCodonPosition(1, i)) *
                             NucMatrix->Stationary(GetCodonPosition(2, i)) *
                             GetFitness(GetCodonStateSpace()->Translation(i)) *
                             GetCodonUsageSelection(i) / total;
            std::cerr << '\n';
        }
        std::cerr << '\n';
        for (int i = 0; i < GetNstate(); i++) {
            for (int j = 0; j < GetNstate(); j++) {
                if (i != j) {
                    int pos = GetDifferingPosition(i, j);
                    if ((pos != -1) && (pos != 3)) {
                        int a = GetCodonPosition(pos, i);
                        int b = GetCodonPosition(pos, j);

                        double mut = (*NucMatrix)(a, b);

                        double deltaS;
                        if (!Synonymous(i, j)) {
                            deltaS =
                                (log(GetFitness(GetCodonStateSpace()->Translation(j))) -
                                 log(GetFitness(GetCodonStateSpace()->Translation(i)))) +
                                (log(GetCodonUsageSelection(j)) - log(GetCodonUsageSelection(i)));
                        } else {
                            deltaS =
                                log(GetCodonUsageSelection(j)) - log(GetCodonUsageSelection(i));
                        }
                        double fix = 1;
                        if (deltaS != 0) {
                            fix = deltaS / (1.0 - exp(-deltaS));
                        }
                        std::cerr << i << '\t' << j << '\t' << mut << '\t' << deltaS << '\t' << fix
                                  << '\t' << mut * fix << '\t' << Q[i][j] << '\n';
                    }
                }
            }
        }
    }
*/

    // data members
    const double *fitness;
    const double *codonusageselection;
};

#endif
