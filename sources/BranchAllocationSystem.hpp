
#ifndef BRANCHALLOC_H
#define BRANCHALLOC_H

class BranchAllocationSystem {
  public:
    BranchAllocationSystem(const Tree& intree, int inNcond)
        : tree(intree), Ncond(inNcond), branchalloc(intree.GetNbranch(), 0), Nbranch(intree.GetNbranch()) {
        MakeBranchAllocations();
    }

    int GetBranchAlloc(int j) const { return branchalloc[j]; }

    const Tree& GetTree() const { return tree; }

    int GetNbranch() const { return Nbranch; }
    const vector<int>& GetAllocVector() const { return branchalloc; }

    void MakeBranchAllocations() {
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

    void RecursiveMakeBranchAllocations(const Link* from) {
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
        for (const Link* link = from->Next(); link != from; link = link->Next()) {
            RecursiveMakeBranchAllocations(link->Out());
        }
    }

  private:
    const Tree& tree;
    int Ncond;
    vector<int> branchalloc;
    int Nbranch;
};

#endif
