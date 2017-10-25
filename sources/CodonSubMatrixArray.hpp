
#ifndef CODONMATRIXARRAY_H
#define CODONMATRIXARRAY_H

#include "GTRSubMatrix.hpp"
#include "CodonSubMatrix.hpp"
#include "AAMutSelOmegaCodonSubMatrix.hpp"
#include "Array.hpp"

class MGOmegaCodonSubMatrixArray : public Array<SubMatrix>, public Array<MGOmegaCodonSubMatrix>	{

	public:
	MGOmegaCodonSubMatrixArray(const CodonStateSpace* incodonstatespace, const GTRSubMatrix* innucmatrix, const ConstArray<double>* inomegaarray) : codonstatespace(incodonstatespace), nucmatrix(innucmatrix), omegaarray(inomegaarray), matrixarray(inomegaarray->GetSize())   {
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

    void UpdateCodonMatrices(const vector<int>& occupancy)  {
        if (((int) occupancy.size()) != GetSize())  {
            cerr << "error in UpdateCodonMatrices: occupancy vector size does not match array size\n";
            exit(1);
        }
		for (int i=0; i<GetSize(); i++)	{
            if (occupancy[i])   {
                (*this)[i].SetOmega(omegaarray->GetVal(i));
                (*this)[i].CorruptMatrix();
            }
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
	const ConstArray<double>* omegaarray;
    vector<MGOmegaCodonSubMatrix*> matrixarray;
};

class AAMutSelOmegaCodonSubMatrixArray : public Array<SubMatrix>, public Array<AAMutSelOmegaCodonSubMatrix>	{

	public:
	AAMutSelOmegaCodonSubMatrixArray(const CodonStateSpace* incodonstatespace, const GTRSubMatrix* innucmatrix, const ConstArray<vector<double> >* inaafitnessarray, double inomega) : codonstatespace(incodonstatespace), nucmatrix(innucmatrix), aafitnessarray(inaafitnessarray), omega(inomega), omegaarray(0), matrixarray(inaafitnessarray->GetSize())   {
        Create();
	}

	AAMutSelOmegaCodonSubMatrixArray(const CodonStateSpace* incodonstatespace, const GTRSubMatrix* innucmatrix, const ConstArray<vector<double> >* inaafitnessarray, const ConstArray<double>* inomegaarray) : codonstatespace(incodonstatespace), nucmatrix(innucmatrix), aafitnessarray(inaafitnessarray), omegaarray(inomegaarray), matrixarray(inomegaarray->GetSize())   {
        if (aafitnessarray->GetSize() != omegaarray->GetSize()) {
            cerr << "error in constructor of AAMutSelOmegaCodonSubMatrixArray: arrays of aafitness and omega values should be of same size\n";
            exit(1);
        }
        Create();
	}

	~AAMutSelOmegaCodonSubMatrixArray()	{
        Delete();
	}
		
    void SetOmega(double inomega)   {
        if (omegaarray) {
            cerr << "error in AAMutSelOmegaCodonSubMatrixArray::SetOmega\n";
            exit(1);
        }
        omega = inomega;
    }

    int GetSize() const {return aafitnessarray->GetSize();}
    const AAMutSelOmegaCodonSubMatrix& GetVal(int i) const {return *matrixarray[i];}
    AAMutSelOmegaCodonSubMatrix& operator[](int i) {return *matrixarray[i];}

	const GTRSubMatrix& GetNucMatrix() const {return *nucmatrix;}

    void UpdateCodonMatrices()  {
        if (omegaarray) {
            for (int i=0; i<GetSize(); i++)	{
                (*this)[i].SetOmega(omegaarray->GetVal(i));
                (*this)[i].CorruptMatrix();
            }
        }
        else    {
            for (int i=0; i<GetSize(); i++)	{
                (*this)[i].SetOmega(omega);
                (*this)[i].CorruptMatrix();
            }
        }
    }

    void UpdateCodonMatrices(const vector<int>& occupancy)  {
        if (omegaarray) {
            for (int i=0; i<GetSize(); i++)	{
                if (! occupancy[i]) {
                    (*this)[i].SetOmega(omegaarray->GetVal(i));
                    (*this)[i].CorruptMatrix();
                }
            }
        }
        else    {
            for (int i=0; i<GetSize(); i++)	{
                if (! occupancy[i]) {
                    (*this)[i].SetOmega(omega);
                    (*this)[i].CorruptMatrix();
                }
            }
        }
    }

	private:

    void Create()   {
		for (int i=0; i<GetSize(); i++)	{
            if (omegaarray) {
                matrixarray[i] = new AAMutSelOmegaCodonSubMatrix(codonstatespace,nucmatrix,aafitnessarray->GetVal(i),omegaarray->GetVal(i));
            }
            else    {
                matrixarray[i] = new AAMutSelOmegaCodonSubMatrix(codonstatespace,nucmatrix,aafitnessarray->GetVal(i),omega);
            }
		}
    }

    void Delete()   {
		for (int i=0; i<GetSize(); i++)	{
			delete matrixarray[i];
		}
    }
        

	const CodonStateSpace* codonstatespace;
	const GTRSubMatrix* nucmatrix;
    const ConstArray<vector<double> >* aafitnessarray;
    double omega;
	const ConstArray<double>* omegaarray;
    vector<AAMutSelOmegaCodonSubMatrix*> matrixarray;
};

#endif
