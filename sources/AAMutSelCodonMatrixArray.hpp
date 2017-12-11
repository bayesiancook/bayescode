
#ifndef AAMUTSELCODONMATRIXARRAY_H
#define AAMUTSELCODONMATRIXARRAY_H

#include "BidimArray.hpp"
#include "AAMutSelOmegaCodonSubMatrix.hpp"

class AAMutSelCodonMatrixArray : public BidimArray<SubMatrix>, public BidimArray<AAMutSelOmegaCodonSubMatrix>  {

    public:

    AAMutSelCodonMatrixArray(const BidimSelector<vector<double> >& infitnessarray, const CodonStateSpace& incodonstatespace, const SubMatrix& innucmatrix) :
        fitnessarray(infitnessarray), codonstatespace(incodonstatespace), nucmatrix(innucmatrix),
        matrixarray(infitnessarray.GetNrow(),vector<AAMutSelOmegaCodonSubMatrix*>(infitnessarray.GetNcol(),0)) {
        Create();
    }

    ~AAMutSelCodonMatrixArray()   {
        Delete();
    }

	virtual int GetNrow() const override {return fitnessarray.GetNrow();}
    virtual int GetNcol() const override {return fitnessarray.GetNcol();}
	virtual const AAMutSelOmegaCodonSubMatrix& GetVal(int i, int j) const override {return *matrixarray[i][j];}
	virtual AAMutSelOmegaCodonSubMatrix& operator()(int i, int j) override {return *matrixarray[i][j];}

    void Create()   {
        for (int i=0; i<GetNrow(); i++)  {
            for (int j=0; j<GetNcol(); j++)   {
                matrixarray[i][j] = new AAMutSelOmegaCodonSubMatrix(&codonstatespace,&nucmatrix,fitnessarray.GetVal(i,j),1.0,false);
            }
        }
    }

    void Delete()   {
        for (int i=0; i<GetNrow(); i++)  {
            for (int j=0; j<GetNcol(); j++)   {
                delete matrixarray[i][j];
            }
        }
    }

    void Corrupt()    {
        for (int i=0; i<GetNrow(); i++)  {
            for (int j=0; j<GetNcol(); j++)   {
                matrixarray[i][j]->CorruptMatrix();
            }
        }
    }

    void CorruptColumn(int j)    {
        for (int i=0; i<GetNrow(); i++)  {
            matrixarray[i][j]->CorruptMatrix();
        }
    }

    void CorruptColumn(int j, const vector<int>& flag)  {
        for (int i=0; i<GetNrow(); i++) {
            if (flag[i])    {
                matrixarray[i][j]->CorruptMatrix();
            }
        }
    }

    private:

    const BidimSelector<vector<double> >& fitnessarray;
    const CodonStateSpace& codonstatespace;
    const SubMatrix& nucmatrix;
    vector<vector<AAMutSelOmegaCodonSubMatrix*> > matrixarray;
};

#endif
