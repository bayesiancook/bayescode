
#ifndef MATRIXSELECTOR_H
#define MATRIXSELECTOR_H

#include "Array.hpp"
#include "BidimArray.hpp"
#include "BranchAllocationSystem.hpp"
#include "SubMatrix.hpp"

class SubMatrixSelector : public BranchSiteSelector<SubMatrix> {
  public:
    SubMatrixSelector(const BidimSelector<SubMatrix>& inmatrixbidimarray, const BranchAllocationSystem& inbranchalloc)
        : matrixbidimarray(inmatrixbidimarray), branchalloc(inbranchalloc) {}

    virtual const Tree& GetTree() const override { return branchalloc.GetTree(); }
    virtual int GetSize() const override { return matrixbidimarray.GetNcol(); }

    virtual const SubMatrix& GetVal(int branch, int site) const override {
        return matrixbidimarray.GetVal(branchalloc.GetBranchAlloc(branch), site);
    }

  private:
    const BidimSelector<SubMatrix>& matrixbidimarray;
    const BranchAllocationSystem& branchalloc;
};

class RootSubMatrixSelector : public Selector<SubMatrix> {
  public:
    RootSubMatrixSelector(const BidimSelector<SubMatrix>& inmatrixbidimarray) : matrixbidimarray(inmatrixbidimarray) {}

    ~RootSubMatrixSelector() {}

    virtual int GetSize() const override { return matrixbidimarray.GetNcol(); }
    virtual const SubMatrix& GetVal(int site) const override { return matrixbidimarray.GetVal(0, site); }

  private:
    const BidimSelector<SubMatrix>& matrixbidimarray;
};

#endif
