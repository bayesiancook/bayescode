
#ifndef BRANCHARRAY_H
#define BRANCHARRAY_H

#include "Tree.hpp"
#include <vector>
#include "MPIBuffer.hpp"

template<class T> class BranchSelector	{

	public:
	virtual ~BranchSelector() {}

	int GetNbranch() const {return GetTree().GetNbranch();}

	virtual const Tree& GetTree() const = 0;
	virtual const T& GetVal(int index) const = 0;

    unsigned int GetMPISize() const {return this->GetNbranch() * MPISize(this->GetVal(0));}

    void MPIPut(MPIBuffer& buffer) const {
        for (int i=0; i<this->GetNbranch(); i++)  {
            buffer << this->GetVal(i);
        }
    }

    void ToStream(ostream& os) const    {
        for (int i=0; i<this->GetNbranch(); i++)  {
            os << this->GetVal(i) << '\t';
        }
        os << '\n';
    }
};

template<class T> class BranchArray : public BranchSelector<T>	{

	public:
	virtual ~BranchArray() {}

    void Copy(const BranchSelector<T>& from)  {

        if (this->GetNbranch() != from.GetNbranch())    {
            cerr << "error: branch arrays do not have same size\n";
            exit(1);
        }
        for (int i=0; i<this->GetNbranch(); i++)  {
            (*this)[i] = from.GetVal(i);
        }
    }

	virtual T& operator[](int index) = 0;

    void MPIGet(const MPIBuffer& buffer)    {
        for (int i=0; i<this->GetNbranch(); i++)  {
            buffer >> (*this)[i];
        }
    }

    void FromStream(istream& is)    {
        for (int i=0; i<this->GetNbranch(); i++)  {
            is >> (*this)[i];
        }
    }
};

template<class T> ostream& operator<<(ostream& os, const BranchSelector<T>& array)  {
    array.ToStream(os);
    return os;
}

template<class T> istream& operator>>(istream& is, BranchArray<T>& array) {
    array.FromStream(is);
    return is;
}

template<class T> class BranchHomogeneousSelector : public BranchSelector<T>	{

	public:
	BranchHomogeneousSelector(const Tree* intree, const T& invalue) : tree(intree), value(invalue) {}
	~BranchHomogeneousSelector() {}

	const Tree& GetTree() const /*override*/ {return tree;}
	const T& GetVal(int index) const /*override*/ {return value;}

	private:
	const Tree& tree;
	const T& value;
};

template<class T> class SimpleBranchArray : public BranchArray<T>	{

	public:
	SimpleBranchArray(const Tree& intree) : tree(intree), array(intree.GetNbranch()) {}
	SimpleBranchArray(const Tree& intree, const T& initval) : tree(intree), array(intree.GetNbranch(),initval) {}
	virtual ~SimpleBranchArray() {}

	const Tree& GetTree() const /*override*/ {return tree;}
	T& operator[](int index) /*override*/ {return array[index];}
	const T& GetVal(int index) const /*override*/ {return array[index];}

	private:
	const Tree& tree;
	vector<T> array;
};

#endif
