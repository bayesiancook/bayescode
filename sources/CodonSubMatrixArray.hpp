
#ifndef CODONMATRIXARRAY_H
#define CODONMATRIXARRAY_H

#include "GTRSubMatrix.hpp"
#include "CodonSubMatrix.hpp"
#include "AAMutSelOmegaCodonSubMatrix.hpp"
#include "Array.hpp"

class MGOmegaCodonSubMatrixArray : public Array<SubMatrix>, public Array<MGOmegaCodonSubMatrix>	{

	public:
	MGOmegaCodonSubMatrixArray(const CodonStateSpace* incodonstatespace, const GTRSubMatrix* innucmatrix, const Array<double>* inomegaarray) : codonstatespace(incodonstatespace), nucmatrix(innucmatrix), omegaarray(inomegaarray), matrixarray(inomegaarray->GetSize())   {
        Create();
	}

	~MGOmegaCodonSubMatrixArray()	{
        Delete();
	}
		
    int GetSize() const {return omegaarray->GetSize();}
    const MGOmegaCodonSubMatrix& GetVal(int i) const {return *matrixarray[i];}
    MGOmegaCodonSubMatrix& operator[](int i) {return *matrixarray[i];}

	const GTRSubMatrix& GetNucMatrix() const {return *nucmatrix;}

    void UpdateCodonMatrices()  {
		for (int i=0; i<GetSize(); i++)	{
            (*this)[i].SetOmega(omegaarray->GetVal(i));
            (*this)[i].CorruptMatrix();
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

class AAMutSelOmegaCodonSubMatrixArray : public Array<SubMatrix>, public Array<AAMutSelOmegaCodonSubMatrix>	{

	public:
	AAMutSelOmegaCodonSubMatrixArray(const CodonStateSpace* incodonstatespace, const GTRSubMatrix* innucmatrix, const Array<vector<double> >* inaafitnessarray, double inomega) : codonstatespace(incodonstatespace), nucmatrix(innucmatrix), aafitnessarray(inaafitnessarray), omega(inomega), matrixarray(inaafitnessarray->GetSize())   {
	// AAMutSelOmegaCodonSubMatrixArray(const CodonStateSpace* incodonstatespace, const GTRSubMatrix* innucmatrix, const Array<vector<double> >* inaafitnessarray, const Array<double>* inomegaarray) : codonstatespace(incodonstatespace), nucmatrix(innucmatrix), aafitnessarray(inaafitnessarray), omegaarray(inomegaarray), matrixarray(inomegaarray->GetSize())   {
        Create();
	}

	~AAMutSelOmegaCodonSubMatrixArray()	{
        Delete();
	}
		
    void SetOmega(double inomega)   {
        omega = inomega;
    }

    int GetSize() const {return aafitnessarray->GetSize();}
    const AAMutSelOmegaCodonSubMatrix& GetVal(int i) const {return *matrixarray[i];}
    AAMutSelOmegaCodonSubMatrix& operator[](int i) {return *matrixarray[i];}

	const GTRSubMatrix& GetNucMatrix() const {return *nucmatrix;}

    void UpdateCodonMatrices()  {
		for (int i=0; i<GetSize(); i++)	{
            (*this)[i].SetOmega(omega);
            // (*this)[i].SetOmega(omegaarray->GetVal(i));
            (*this)[i].CorruptMatrix();
		}
    }

	private:

    void Create()   {
		for (int i=0; i<GetSize(); i++)	{
			matrixarray[i] = new AAMutSelOmegaCodonSubMatrix(codonstatespace,nucmatrix,aafitnessarray->GetVal(i),omega);
			// matrixarray[i] = new AAMutSelOmegaCodonSubMatrix(codonstatespace,nucmatrix,aafitnessarray->GetVal(i),omegaarray->GetVal(i));
		}
    }

    void Delete()   {
		for (int i=0; i<GetSize(); i++)	{
			delete matrixarray[i];
		}
    }
        

	const CodonStateSpace* codonstatespace;
	const GTRSubMatrix* nucmatrix;
    const Array<vector<double> >* aafitnessarray;
	// const Array<double>* omegaarray;
    double omega;
    vector<AAMutSelOmegaCodonSubMatrix*> matrixarray;
};

#endif
