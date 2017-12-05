
#ifndef BRANCHSITEARRAY_H
#define BRANCHSITEARRAY_H

#include "Array.hpp"
#include "BranchArray.hpp"

template<class T> class BranchSiteSelector	{

	public:
	virtual ~BranchSiteSelector() {}

	virtual const Tree& GetTree() const  = 0;
	int GetNbranch() const {return GetTree().GetNbranch();}

	virtual int GetSize() const  = 0;
	virtual const T& GetVal(int branch, int site) const = 0;
};

template<class T> class BranchHomogeneousSiteHomogeneousSelector : public BranchSiteSelector<T> {

	public:
	BranchHomogeneousSiteHomogeneousSelector(const Tree& intree, int insize, const T& invalue) : tree(intree), size(insize), value(invalue) {}
	~BranchHomogeneousSiteHomogeneousSelector() {}

    const Tree& GetTree() const /*override*/ {return tree;}
    int GetSize() const /*override*/ {return size;}
	const T& GetVal(int branch, int site) const /*override*/ {return value;}

	private:
    const Tree& tree;
    int size;
	const T& value;
};

template<class T> class BranchHomogeneousSiteHeterogeneousSelector : public BranchSiteSelector<T>	{

	public:
	BranchHomogeneousSiteHeterogeneousSelector(const Tree& intree, const Selector<T>& inarray) : tree(intree), array(inarray) {}
	~BranchHomogeneousSiteHeterogeneousSelector() {}

    const Tree& GetTree() const /*override*/ {return tree;}
    int GetSize() const /*override*/ {return array.GetSize();}
	const T& GetVal(int branch, int site) const /*override*/ {return array.GetVal(site);}

	private:
    const Tree& tree;
	const Selector<T>& array;
};

template<class T> class BranchHeterogeneousSiteHomogeneousSelector : public BranchSiteSelector<T> {

	public:
	BranchHeterogeneousSiteHomogeneousSelector(const BranchSelector<T>& inbrancharray, int insize) : brancharray(inbrancharray), size(insize) {}
	~BranchHeterogeneousSiteHomogeneousSelector() {}

    const Tree& GetTree() const /*override*/ {return brancharray.GetTree();}
    int GetSize() const /*override*/ {return size;}
	const T& GetVal(int branch, int site) const /*override*/ {return brancharray.GetVal(branch);}

	private:
	const BranchSelector<T>& brancharray;
    int size;
};

#endif
