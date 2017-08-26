
#ifndef BRANCHARRAY_H
#define BRANCHARRAY_H

#include "Tree.hpp"
#include <vector>

template<class T> class ConstBranchArray	{

	public:
	virtual ~ConstBranchArray() {}

	int GetNbranch() const {return GetTree().GetNbranch();}

	virtual const Tree& GetTree() const = 0;
	virtual const T& GetVal(int index) const = 0;

};

template<class T> class BranchArray : public ConstBranchArray<T>	{

	public:
	virtual ~BranchArray() {}

	virtual T& operator[](int index) = 0;
};

template<class T> class HomogeneousBranchArray : public ConstBranchArray<T>	{

	public:
	HomogeneousBranchArray(const Tree* intree, const T& invalue) : tree(intree), value(invalue) {}
	~HomogeneousBranchArray() {}

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
