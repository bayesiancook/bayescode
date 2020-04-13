#pragma once

#include "AAMutSelPolyMutOmegaCodonSubMatrix.hpp"
#include "BidimArray.hpp"

class AAMutSelPolyMutCodonMatrixBidimArray : public BidimArray<SubMatrix>,
                                       public BidimArray<AAMutSelPolyMutOmegaCodonSubMatrix> {
  public:

    //! constructor parameterized by a bidim array of fitness profiles, a codon
    //! state space and a single nucleotide matrix.
    AAMutSelPolyMutCodonMatrixBidimArray(const BidimSelector<vector<double>> &infitnessarray,
                                   const CodonStateSpace &incodonstatespace,
                                   const SubMatrix &innucmatrix,
                                   double inu = 0,
                                   double inomega = 1.0)
        : fitnessarray(infitnessarray),
          codonstatespace(incodonstatespace),
          nucmatrix(innucmatrix),
          u(inu),
          omega(inomega),
          matrixarray(infitnessarray.GetNrow(),
                      vector<AAMutSelPolyMutOmegaCodonSubMatrix *>(infitnessarray.GetNcol(),
                                                            (AAMutSelPolyMutOmegaCodonSubMatrix *)0)) {
        Create();
    }

    ~AAMutSelPolyMutCodonMatrixBidimArray() { Delete(); }

    //! return the number of rows (number of conditions)
    virtual int GetNrow() const override { return fitnessarray.GetNrow(); }
    //! return the number of columns (number of sites)
    virtual int GetNcol() const override { return fitnessarray.GetNcol(); }
    //! const access to the matrix for row (condition) i and column (site) j
    virtual const AAMutSelPolyMutOmegaCodonSubMatrix &GetVal(int i, int j) const override {
        return *matrixarray[i][j];
    }
    //! non const access to the matrix for row (condition) i and column (site) j
    virtual AAMutSelPolyMutOmegaCodonSubMatrix &operator()(int i, int j) override {
        return *matrixarray[i][j];
    }

    //! allocation and construction of all matrices
    void Create() {
        for (int i = 0; i < GetNrow(); i++) {
            for (int j = 0; j < GetNcol(); j++) {
                matrixarray[i][j] = new AAMutSelPolyMutOmegaCodonSubMatrix(
                    &codonstatespace, &nucmatrix, fitnessarray.GetVal(i, j), u, omega, 1.0);
            }
        }
    }

    //! destruction of all matrices
    void Delete() {
        for (int i = 0; i < GetNrow(); i++) {
            for (int j = 0; j < GetNcol(); j++) {
                delete matrixarray[i][j];
            }
        }
    }

    void SetU(double inu)   {
        u = inu;
    }

    void SetOmega(double inomega)   {
        omega = inomega;
    }

    //! signal corruption of the parameters (matrices should recompute themselves)
    void Corrupt() {
        for (int i = 0; i < GetNrow(); i++) {
            for (int j = 0; j < GetNcol(); j++) {
                matrixarray[i][j]->SetU(u);
                matrixarray[i][j]->SetOmega(omega);
                matrixarray[i][j]->CorruptMatrix();
            }
        }
    }

    //! signal corruption for row i
    void CorruptRow(int i) {
        for (int j = 0; j < GetNcol(); j++) {
            matrixarray[i][j]->SetU(u);
            matrixarray[i][j]->SetOmega(omega);
            matrixarray[i][j]->CorruptMatrix();
        }
    }

    //! signal corruption for column (site) j
    void CorruptColumn(int j) {
        for (int i = 0; i < GetNrow(); i++) {
            matrixarray[i][j]->SetU(u);
            matrixarray[i][j]->SetOmega(omega);
            matrixarray[i][j]->CorruptMatrix();
        }
    }

    //! signal corruption for column (site) j, and only for those rows
    //! (conditions) that are flagged
    void CorruptColumn(int j, const vector<int> &flag) {
        for (int i = 0; i < GetNrow(); i++) {
            if (flag[i]) {
                matrixarray[i][j]->SetU(u);
                matrixarray[i][j]->SetOmega(omega);
                matrixarray[i][j]->CorruptMatrix();
            }
        }
    }

    // assuming uniform weights
    void GetMeanAAFrequencies(vector<vector<double>>& compaafreq) const {

        if (GetNrow() != Naa)   {
            cerr << "error in get mean aa freqs\n";
            exit(1);
        }

        for (int a=0; a<Naa; a++)   {
            for (int b=0; b<Naa; b++)   {
                compaafreq[a][b] = 0;
            }
        }
        for (int i = 0; i < GetNrow(); i++) {
            for (int j=0; j<GetNcol(); j++) {
                const AAMutSelPolyMutOmegaCodonSubMatrix& mat = GetVal(i,j);
                for (int k=0; k<mat.GetNstate(); k++)   {
                    compaafreq[i][mat.GetCodonStateSpace()->Translation(k)] += mat.Stationary(k);
                }
            }
        }
        for (int i = 0; i < GetNrow(); i++) {
            for (int a=0; a<Naa; a++)   {
                compaafreq[i][a] /= GetNcol();
            }
        }
    }

  private:
    const BidimSelector<vector<double>> &fitnessarray;
    const CodonStateSpace &codonstatespace;
    const SubMatrix &nucmatrix;
    double u;
    double omega;
    vector<vector<AAMutSelPolyMutOmegaCodonSubMatrix *>> matrixarray;
};

