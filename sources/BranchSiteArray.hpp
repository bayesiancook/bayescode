
#ifndef BRANCHSITEARRAY_H
#define BRANCHSITEARRAY_H

#include "Array.hpp"
#include "BranchArray.hpp"

template<class T> class ConstBranchSiteArray	{

	public:
	virtual ~ConstBranchSiteArray() {}

	const Tree* GetTree() const  = 0;
	int GetNbranch() const {return tree->GetNbranch();}

	int GetSize() const  = 0;
	virtual const T& GetVal(int branch, int site) const = 0;
};

template<class T> class BranchSiteArray : public ConstBranchSiteArray<T> {

	public:
	virtual ~BranchSiteArray() {}

	virtual T& operator()(int branch, int site) = 0;
};

template<class T> class HomogeneousBranchSiteArray : public virtual ConstBranchSiteArray<T> {

	public:
	HomogeneousBranchSiteArray(const Tree* intree, int insize, const T& invalue) : tree(intree), size(insize), value(invalue) {}
	~HomogeneousBranchSiteArray() {}

    const Tree* GetTree() const override {return tree;}
    int GetSize() const override {return size;}
	const T& GetVal(int branch, int site) const override {return value;}

	private:
    const Tree* tree;
    int size;
	const T& value;
};

template<class T> class BranchHomogeneousSiteHeterogeneousArray : public virtual ConstBranchSiteArray<T>	{

	public:
	BranchHomogeneousSiteHeterogeneousArray(const Tree* intree, const ConstArray<T>* inarray) : tree(intree), array(inarray) {}
	~BranchHomogeneousSiteHeterogeneousArray() {}

    const Tree* GetTree() const override {return tree;}
    int GetSize() const override {return array->GetSize();}
	const T& GetVal(int branch, int site) const override {return array->GetVal(site);}

	private:
    const Tree* tree;
	const ConstArray<T>* array;
};

template<class T> class BranchHeterogeneousSiteHomogeneousArray : public virtual ConstBranchSiteArray<T> {

	public:
	BranchHeterogeneousSiteHomogeneousArray(const ConstBranchArray<T>* inbrancharray, int insize) : brancharray(inbrancharray), size(insize) {}
	~BranchHeterogeneousSiteHomogeneousArray() {}

    const Tree* GetTree() const override {return brancharray->GetTree();}
    int GetSize() const override {return size;}
	const T& GetVal(int branch, int site) const override {return brancharray->GetVal(branch);}

	private:
	const ConstBranchArray<T>* brancharray;
    int size;
};

#endif
