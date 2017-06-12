
#ifndef CODONMATRIXARRAY_H
#define CODONMATRIXARRAY_H

#include "GTRSubMatrix.hpp"
#include "CodonSubMatrix.hpp"
#include "Array.hpp"

class MGOmegaHeterogeneousCodonSubMatrixArray : public Array<SubMatrix>, public Array<MGOmegacodonSubMatrix>	{

	public:
	MGOmegaHeterogeneousCodonSubMatrixArray(const CodonStateSpace* incodonstatespace, const GTRSubMatrix* innucmatrix, const Array<double>* inomegaarray) : codonstatespace(incodonstatespace), nucmatrix(innucmatrix), omegaarray(inomegaarray), matrixarray(inomegaarray->GetSize())   {
        Create();
	}

	~MGOmegaHeterogeneousCodonSubMatrixArray()	{
        Delete();
	}
		
    int GetSize() const {return omegaarray->GetSize();}
    const MGOmegaCodonSubMatrix& GetVal(int i) const {return *matrixarray[i];}
    MGOmegaCodonSubMatrix& operator[](int i) {return *matrixarray[i];}

	const GTRSubMatrix& GetNucMatrix() const {return *nucmatrix;}

    void UpdateCodonMatrices()  {
		for (int i=0; i<GetSize(); i++)	{
            (*this)[i]->SetOmega(omegaarray->GetVal(i));
            (*this)[i]->CorruptMatrix();
		}
    }

	private:

    void Create()   {
		for (int i=0; i<GetSize(); i++)	{
			matrixarray[i] = new MGOmegaCodonSubMatrix(codonstatespace,nucmatrix,omegaarray->GetVal(i));
		}
    }

    void Delete()   {
		for (int i=0; i<GetSize(); i++)	{
			delete matrixarray[i];
		}
    }
        

	const CodonStateSpace* codonstatespace;
	const GTRSubMatrix* nucmatrix;
	const Array<double>* omegaarray;
    vector<MGOmegaCodonSubMatrix*> matrixarray;
};

#endif
