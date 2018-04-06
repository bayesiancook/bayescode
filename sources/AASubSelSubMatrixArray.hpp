
#ifndef AASUBSELMATRIXARRAY_H
#define AASUBSELMATRIXARRAY_H

#include "AASubSelSubMatrix.hpp"
#include "Array.hpp"

/**
 * \brief An Array of AASubSelSubMatrix
 */

class AASubSelSubMatrixArray : public Array<SubMatrix>, public Array<AASubSelSubMatrix>	{

	public:
	AASubSelSubMatrixArray(const vector<double>& inrelrate, const Selector<vector<double> >* inprofilearray, int innormalise) : relrate(inrelrate), profilearray(inprofilearray), nstate(inprofilearray->GetVal(0).size()), normalise(innormalise),  matrixarray(inprofilearray->GetSize()) {
        Create();
	}

	~AASubSelSubMatrixArray()	{
        Delete();
	}
		
    //! return array size
    int GetSize() const {return profilearray->GetSize();}
    //! const access to matrix i
    const AASubSelSubMatrix& GetVal(int i) const {return *matrixarray[i];}
    //! non-const access to matrix i
    AASubSelSubMatrix& operator[](int i) {return *matrixarray[i];}

    //! update all matrices
    void UpdateMatrices()  {
        for (int i=0; i<GetSize(); i++)	{
            (*this)[i].CopyStationary(profilearray->GetVal(i));
            (*this)[i].CorruptMatrix();
        }
    }

    void UpdateMatrix(int i)    {
        (*this)[i].CopyStationary(profilearray->GetVal(i));
        (*this)[i].CorruptMatrix();
    }

	private:

    void Create()   {
		for (int i=0; i<GetSize(); i++)	{
            matrixarray[i] = new AASubSelSubMatrix(nstate,relrate,profilearray->GetVal(i),normalise);
		}
    }

    void Delete()   {
		for (int i=0; i<GetSize(); i++)	{
			delete matrixarray[i];
		}
    }
        

    const vector<double>& relrate;
    const Selector<vector<double> >* profilearray;
    int nstate;
    int normalise;
    vector<AASubSelSubMatrix*> matrixarray;
};

#endif
