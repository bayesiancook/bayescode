#ifndef BRANCHALLOC_H
#define BRANCHALLOC_H

#include "Tree.hpp"

/**
 * \brief An object that associates an integer (in 0..K-1) to each branch of a
 * phylogenetic tree
 *
 * should perhaps derive it from BranchArray
 *
 */

class BranchAllocationSystem {
  public:
    //! \brief Constructor with tree and number of conditions
    //!
    //! names of branches of the tree are assumed to encode integers, specifying
    //! the allocation for each branch
    BranchAllocationSystem(const Tree &intree, int inNcond)
        : tree(intree),
          Ncond(inNcond),
          branchalloc(intree.GetNbranch(), 0),
          Nbranch(intree.GetNbranch()) {
        MakeBranchAllocations();
    }

    //! return allocation status of branch j
    int GetBranchAlloc(int j) const { return branchalloc[j]; }

    //! return a const ref to underlying tree
    const Tree &GetTree() const { return tree; }

    //! return number of branches
    int GetNbranch() const { return Nbranch; }

    //! return array of branch allocations as a simple vector<int>
    const vector<int> &GetAllocVector() const { return branchalloc; }

  private:
    //! read out branch names (recursively) and fill-in allocation map
    void MakeBranchAllocations() {

        if (Ncond == Nbranch)   {
            for (int j = 0; j < Nbranch; j++) {
                branchalloc[j] = j;
            }
        }

        else    {
            // default pre-initialization
            for (int j = 0; j < Nbranch; j++) {
                branchalloc[j] = -1;
            }

            RecursiveMakeBranchAllocations(tree.GetRoot());

            // check that all branches have been correctly initialized
            for (int j = 0; j < Nbranch; j++) {
                if ((branchalloc[j] < 0) || (branchalloc[j] >= Ncond)) {
                    std::cerr << "error in make branch allocation\n";
                    cerr << j << '\t' << branchalloc[j] << '\n';
                    exit(1);
                }
            }
        }
    }

    //! recursive helper function for MakeBranchAllocations
    void RecursiveMakeBranchAllocations(const Link *from) {
        if (!from->isRoot()) {
            int k = atoi(from->GetBranch()->GetName().c_str());
            if (k >= Ncond) {
                k = Ncond - 1;
            }
            if (k < 0) {
                std::cerr << "error : allocation out of bound\n";
                std::cerr << "k" << '\t' << "Ncond" << '\n';
                exit(1);
            }
            branchalloc[from->GetBranch()->GetIndex()] = k;
        }
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            RecursiveMakeBranchAllocations(link->Out());
        }
    }

  private:
    const Tree &tree;
    int Ncond;
    vector<int> branchalloc;
    int Nbranch;
};

#endif
