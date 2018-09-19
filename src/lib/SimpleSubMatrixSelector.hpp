
#ifndef MATRIXSELECTOR_H
#define MATRIXSELECTOR_H

#include "Array.hpp"
#include "CodonSubMatrix.hpp"

/**
 * \brief A BranchSelector for DiffSelModel and DiffSelSparseModel
 *
 * In DiffSelModel and DiffSelSparseModel, the substitution process varies
 * across sites and across the k=0..K-1 conditions. In addition, each branch is
 * allocated to one of the K conditions
 * Accordingly, SimpleMGOmegaCodonSubMatrixSelector returns, for a given branch
 * and a given site, a pointer to the substitution matrix corresponding to that
 * site, and under the condition being active for that branch.
 */

class SimpleSubMatrixSelector : public BranchSelector<SubMatrix> {
  public:
    //! constructor, parameterized by bidim array of substitution matrices and
    //! branch allocation system
    SimpleSubMatrixSelector(
        const Selector<SubMatrix> &inmatrixarray, const BranchSelector<int> &inbranchalloc)
        : matrixarray(inmatrixarray), branchalloc(inbranchalloc) {}

    virtual const Tree &GetTree() const override { return branchalloc.GetTree(); }

    virtual const SubMatrix &GetVal(int branch) const override {
        return matrixarray.GetVal(branchalloc.GetVal(branch));
    }

  private:
    const Selector<SubMatrix> &matrixarray;
    const BranchSelector<int> &branchalloc;
};

class SimpleMGOmegaCodonSubMatrixSelector : public BranchSelector<MGOmegaCodonSubMatrix> {
  public:
    //! constructor, parameterized by bidim array of substitution matrices and
    //! branch allocation system
    SimpleMGOmegaCodonSubMatrixSelector(const Selector<MGOmegaCodonSubMatrix> &inmatrixarray,
        const BranchSelector<int> &inbranchalloc)
        : matrixarray(inmatrixarray), branchalloc(inbranchalloc) {}

    virtual const Tree &GetTree() const override { return branchalloc.GetTree(); }

    virtual const MGOmegaCodonSubMatrix &GetVal(int branch) const override {
        return matrixarray.GetVal(branchalloc.GetVal(branch));
    }

  private:
    const Selector<MGOmegaCodonSubMatrix> &matrixarray;
    const BranchSelector<int> &branchalloc;
};

/**
 * \brief A Selector of substitution matrices at the root, for DiffSelModel and
 * DiffSelSparseModel
 *
 * In DiffSelModel and DiffSelSparseModel, the root of the tree is assumed to be
 * under condition 0. Accorddingly, RootSimpleMGOmegaCodonSubMatrixSelector
 * returns, for a given site, a pointer to the substitution matrix corresponding
 * to that site, and under condition 0.
 */

#endif
