
#ifndef AADIFFSELCODONMATRIXARRAY_H
#define AADIFFSELCODONMATRIXARRAY_H

#include "Array.hpp"
#include "BidimArray.hpp"
#include "AAMutSelOmegaCodonSubMatrix.hpp"
#include "IIDDirichlet.hpp"

/**
 * \brief A BidimArray of AAMutSelOmegaCodonSubMatrix (used in DiffSelModel and DiffSelSparseModel)
 *
 * The constructor takes as arguments a single nucleotide matrix and a bidim array of amino-acid fitness profiles (across sites and conditions).
 * It then constructs an array of AAMutSelOmegaCodonSubMatrix of same size as the array of fitness profiles.
 */

class AAMutSelNeCodonMatrixBidimArray : public BidimArray<SubMatrix>, public BidimArray<AAMutSelOmegaCodonSubMatrix>  {

    public:

    // two possible ways to do it:
    // AAMutSelNeCodonMatrixBidimArray(const Selector<vector<double> >& infitnessarray, const Selector<double>& inNe, const CodonStateSpace& incodonstatespace, const SubMatrix& innucmatrix);
    // AAMutSelNeCodonMatrixBidimArray(const Selector<vector<double> >& infitnessarray, const vector<double>& inNe, const CodonStateSpace& incodonstatespace, const SubMatrix& innucmatrix);

    //! constructor parameterized by a bidim array of fitness profiles, a codon state space and a single nucleotide matrix.
    AAMutSelNeCodonMatrixBidimArray(const MultiDirichlet& infitnessarray, const vector<double>& inNe, const CodonStateSpace& incodonstatespace, const SubMatrix& innucmatrix) :
        fitnessarray(infitnessarray), Ne(inNe), codonstatespace(incodonstatespace), nucmatrix(innucmatrix),
        matrixarray(Ne.size(),vector<AAMutSelOmegaCodonSubMatrix*>(infitnessarray.GetSize(),(AAMutSelOmegaCodonSubMatrix*)0)) {
        Create();
    }

    ~AAMutSelNeCodonMatrixBidimArray()   {
        Delete();
    }

    //! return the number of rows (number of conditions)
	virtual int GetNrow() const override {return Ne.size();}
    //! return the number of columns (number of sites)
    virtual int GetNcol() const override {return fitnessarray.GetSize();}
    //! const access to the matrix for row (condition) i and column (site) j
	virtual const AAMutSelOmegaCodonSubMatrix& GetVal(int i, int j) const override {return *matrixarray[i][j];}
    //! non const access to the matrix for row (condition) i and column (site) j
	virtual AAMutSelOmegaCodonSubMatrix& operator()(int i, int j) override {return *matrixarray[i][j];}

  virtual void SetOmega(double omega) {
    for (int i=0; i<this->GetNrow(); i++)  {
        for (int j=0; j<this->GetNcol(); j++)   {
            (*this)(i,j).SetOmega(omega);
        }
    }
  }

  virtual void SetNe(vector<double> &Ne) {
    for (int i=0; i<this->GetNrow(); i++)  {
        for (int j=0; j<this->GetNcol(); j++)   {
            (*this)(i,j).SetNe(Ne[i]);
        }
    }
  }

    //! allocation and construction of all matrices
    void Create()   {
        for (int i=0; i<GetNrow(); i++)  {
          std::cout << "Ne[i]: "<<Ne[i] <<std::endl;
            for (int j=0; j<GetNcol(); j++)   {
                matrixarray[i][j] = new AAMutSelOmegaCodonSubMatrix(&codonstatespace,&nucmatrix,fitnessarray.GetVal(j),1.0,Ne[i]);
            }
        }
    }

    //! destruction of all matrices
    void Delete()   {
        for (int i=0; i<GetNrow(); i++)  {
            for (int j=0; j<GetNcol(); j++)   {
                delete matrixarray[i][j];
            }
        }
    }

    //! signal corruption of the parameters (matrices should recompute themselves)
    void Corrupt()    {
        for (int i=0; i<GetNrow(); i++)  {
            for (int j=0; j<GetNcol(); j++)   {
                matrixarray[i][j]->CorruptMatrix();
            }
        }
    }

    //! signal corruption for column (site) j
    void CorruptColumn(int j)    {
        for (int i=0; i<GetNrow(); i++)  {
            matrixarray[i][j]->CorruptMatrix();
        }
    }

    //! signal corruption for column (site) j, and only for those rows (conditions) that are flagged
    void CorruptColumn(int j, const vector<int>& flag)  {
        for (int i=0; i<GetNrow(); i++) {
            if (flag[i])    {
                matrixarray[i][j]->CorruptMatrix();
            }
        }
    }

    private:

    const Selector<vector<double> >& fitnessarray;
    const vector<double> Ne;
    const CodonStateSpace& codonstatespace;
    const SubMatrix& nucmatrix;
    vector<vector<AAMutSelOmegaCodonSubMatrix*> > matrixarray;
};

#endif
