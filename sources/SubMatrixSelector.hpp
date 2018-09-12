
#ifndef MATRIXSELECTOR_H
#define MATRIXSELECTOR_H

#include "Array.hpp"
#include "BidimArray.hpp"
#include "BranchAllocationSystem.hpp"
#include "SubMatrix.hpp"

/**
 * \brief A BranchSiteSelector for DiffSelModel and DiffSelSparseModel
 *
 * In DiffSelModel and DiffSelSparseModel, the substitution process varies
 * across sites and across the k=0..K-1 conditions. In addition, each branch is
 * allocated to one of the K conditions (see BranchAllocationSystem).
 * Accordingly, SubMatrixSelector returns, for a given branch and a given site,
 * a pointer to the substitution matrix corresponding to that site, and under
 * the condition being active for that branch.
 */

class SubMatrixSelector : public BranchSiteSelector<SubMatrix> {
  public:
    //! constructor, parameterized by bidim array of substitution matrices and
    //! branch allocation system
    SubMatrixSelector(const BidimSelector<SubMatrix> &inmatrixbidimarray,
        const BranchAllocationSystem &inbranchalloc)
        : matrixbidimarray(inmatrixbidimarray), branchalloc(inbranchalloc) {}

    virtual const Tree &GetTree() const override { return branchalloc.GetTree(); }
    virtual int GetSize() const override { return matrixbidimarray.GetNcol(); }

    virtual const SubMatrix &GetVal(int branch, int site) const override {
        return matrixbidimarray.GetVal(branchalloc.GetBranchAlloc(branch), site);
    }

  private:
    const BidimSelector<SubMatrix> &matrixbidimarray;
    const BranchAllocationSystem &branchalloc;
};

/**
 * \brief A Selector of substitution matrices at the root, for DiffSelModel and
 * DiffSelSparseModel
 *
 * In DiffSelModel and DiffSelSparseModel, the root of the tree is assumed to be
 * under condition 0. Accorddingly, RootSubMatrixSelector returns, for a given
 * site, a pointer to the substitution matrix corresponding to that site, and
 * under condition 0.
 */

class RootSubMatrixSelector : public Selector<SubMatrix> {
  public:
    //! constructor, parameterized by bidim array of substitution matrices
    RootSubMatrixSelector(const BidimSelector<SubMatrix> &inmatrixbidimarray)
        : matrixbidimarray(inmatrixbidimarray) {}

    ~RootSubMatrixSelector() {}

    virtual int GetSize() const override { return matrixbidimarray.GetNcol(); }
    virtual const SubMatrix &GetVal(int site) const override {
        return matrixbidimarray.GetVal(0, site);
    }

  private:
    const BidimSelector<SubMatrix> &matrixbidimarray;
};

#endif
