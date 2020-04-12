#pragma once

#include "BiologicalSequences.hpp"  //FIXME only used for Naa (const int)
#include "CodonSubMatrix.hpp"

class AAProjectedCodonSubMatrix : public virtual SubMatrix {
  public:
    //! constructor parameterized by an array of relative rates (size
    //! Nstate*(Nstate-1)/2) and an array of equilibrium frequencies (size Nstate)
    AAProjectedCodonSubMatrix(const CodonSubMatrix* infrom) : SubMatrix(20,false), from(infrom) {

    }
    ~AAProjectedCodonSubMatrix() override = default;

    void GetRelRateArray(vector<double>& rr)    {
        int i = 0;
        double tot = 0;
        for (int a=0; a<GetNstate(); a++)   {
            for (int b=a+1; b<GetNstate(); b++) {
                double rhoab = Q(a,b) / Stationary(b);
                double rhoba = Q(b,a) / Stationary(a);
                double diff = abs(rhoab - rhoba);
                if (diff > 1e-6)    {
                    cerr << "error in compute rel rates: non reversible\n";
                    exit(1);
                }
                rr[i] = rhoab;
                tot += rhoab;
                i++;
            }
        }
        for (int i=0; i<Naarr; i++)   {
            rr[i] /= tot;
            rr[i] *= 190;
        }
    }

    void GetStationaryArray(vector<double>& stat)   {
        for (int a=0; a<GetNstate(); a++)   {
            stat[a] = Stationary(a);
        }
    }

    void ComputeFullArray() {
        ComputeStationary();
        for (int a=0; a<GetNstate(); a++)   {
            for (int b=0; b<GetNstate(); b++)   {
                Q(a,b) = 0;
            }
        }
        for (int c=0; c<from->GetNstate(); c++) {
            int a = from->GetCodonStateSpace()->Translation(c);
            for (int d=0; d<from->GetNstate(); d++) {
                int b = from->GetCodonStateSpace()->Translation(d);
                if (a != b) {
                    Q(a,b) += from->Stationary(c) * (*from)(c,d);
                }
            }
        }
        for (int a=0; a<GetNstate(); a++)   {
            double stat = Stationary(a);
            double tot = 0;
            for (int b=0; b<GetNstate(); b++)   {
                Q(a,b) /= stat;
                tot += Q(a,b);
            }
            Q(a,a) = -tot;
        }
    }

  protected:

    void ComputeStationary() const override {
        for (int a=0; a<GetNstate(); a++) {
            mStationary[a] = 0;
        }
        for (int c=0; c<from->GetNstate(); c++) {
            int a = from->GetCodonStateSpace()->Translation(c);
            mStationary[a] += from->Stationary(c);
        }
    }

    void ComputeArray(int a) const override {
        for (int b=0; b<GetNstate(); b++)   {
            Q(a,b) = 0;
        }
        for (int c=0; c<from->GetNstate(); c++) {
            int a2 = from->GetCodonStateSpace()->Translation(c);
            if (a == a2)    {
                for (int d=0; d<from->GetNstate(); d++) {
                    int b = from->GetCodonStateSpace()->Translation(d);
                    if (a != b) {
                        Q(a,b) += from->Stationary(c) * (*from)(c,d);
                    }
                }
            }
        }
        double stat = Stationary(a);
        double tot = 0;
        for (int b=0; b<GetNstate(); b++)   {
            Q(a,b) /= stat;
            tot += Q(a,b);
        }
        Q(a,a) = -tot;
    }

  private:
    int rrindex(int i, int j) const {
        return (i < j) ? (2 * Nstate - i - 1) * i / 2 + j - i - 1
                       : (2 * Nstate - j - 1) * j / 2 + i - j - 1;
    }

    const CodonSubMatrix* from;
};

