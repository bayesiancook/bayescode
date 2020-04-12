#pragma once

#include "Array.hpp"
#include "BidimArray.hpp"
#include "CodonSubMatrix.hpp"

class MeanCodonSubMatrixFromBidimArray : public virtual CodonSubMatrix    {

    public:
    MeanCodonSubMatrixFromBidimArray(const CodonStateSpace& instatespace, const BidimSelector<SubMatrix>& inarray, const vector<double>& inrowweights, const vector<double>& incolweights) : 
        SubMatrix(instatespace.GetNstate(), false),
        CodonSubMatrix(&instatespace, false),
        array(inarray),
        rowweights(inrowweights),
        colweights(incolweights)    {
            if (array.GetNrow() != int(rowweights.size())) {
                cerr << "error row weight vector has incorrect dimension\n";
                exit(1);
            }
            if (array.GetNcol() != int(colweights.size())) {
                cerr << "error col weight vector has incorrect dimension\n";
                exit(1);
            }
            // ComputeFullArray();
    }

    void ComputeFullArray() {
        ComputeStationary();
        for (int a=0; a<GetNstate(); a++)   {
            for (int b=0; b<GetNstate(); b++)   {
                Q(a,b) = 0;
            }
        }
        for (int i=0; i<array.GetNrow(); i++) {
            for (int j=0; j<array.GetNcol(); j++) {
                for (int a=0; a<GetNstate(); a++)   {
                    for (int b=0; b<GetNstate(); b++)   {
                        Q(a,b) += rowweights[i] * colweights[j] * array.GetVal(i,j).Stationary(a) * array.GetVal(i,j)(a,b);
                    }
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
            if (fabs(tot) > 1e-6)   {
                cerr << "error in mean codon sub matrix: row does not sum to 0\n";
                exit(1);
            }
        }
    }

    protected:

    void ComputeStationary() const override   {
        for (int a=0; a<GetNstate(); a++) {
            mStationary[a] = 0;
        }

        for (int i=0; i<array.GetNrow(); i++) {
            for (int j=0; j<array.GetNcol(); j++) {
                for (int a=0; a<GetNstate(); a++) {
                    mStationary[a] += rowweights[i] * colweights[j] * array.GetVal(i,j).Stationary(a);
                }
            }
        }

        // check
        double tot = 0;
        for (int a=0; a<GetNstate(); a++) {
            tot += mStationary[a];
        }
        if (fabs(tot-1) > 1e-6) {
            cerr << "error in mean codon sub matrix: mean total mass is not 1\n";
            exit(1);
        }
    }

    void ComputeArray(int a) const override {
        for (int b=0; b<GetNstate(); b++)   {
            Q(a,b) = 0;
        }
        for (int i=0; i<array.GetNrow(); i++) {
            for (int j=0; j<array.GetNcol(); j++) {
                for (int b=0; b<GetNstate(); b++)   {
                    Q(a,b) += rowweights[i] * colweights[j] * array.GetVal(i,j).Stationary(a) * array.GetVal(i,j)(a,b);
                }
            }
        }

        double tot = 0;
        double stat = Stationary(a);
        for (int b=0; b<GetNstate(); b++)   {
            Q(a,b) /= stat;
            tot += Q(a,b);
        }
        if (fabs(tot) > 1e-6)   {
            cerr << "error in mean codon sub matrix: row does not sum to 0\n";
            exit(1);
        }
    }

    const BidimSelector<SubMatrix>& array;
    const vector<double>& rowweights;
    const vector<double>& colweights;

};

