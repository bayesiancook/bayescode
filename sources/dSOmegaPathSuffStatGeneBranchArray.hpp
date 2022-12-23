
#include "dSOmegaPathSuffStat.hpp"

class dSOmegaPathSuffStatGeneBranchArray : public Array<dSOmegaPathSuffStatBranchArray> {
  public:
    //! constructor, parameterized by underlying tree and size (number of genes)
    dSOmegaPathSuffStatGeneBranchArray(const Tree& intree, int insize)
        : tree(intree), array(insize, (dSOmegaPathSuffStatBranchArray *)0) {
        for (int i = 0; i < GetSize(); i++) {
            array[i] = new dSOmegaPathSuffStatBranchArray(tree);
        }
    }

    ~dSOmegaPathSuffStatGeneBranchArray() {
        for (int i = 0; i < GetSize(); i++) {
            delete[] array[i];
        }
    }

    int GetSize() const override { return array.size(); }

    const dSOmegaPathSuffStatBranchArray &GetVal(int i) const override { return *array[i]; }

    dSOmegaPathSuffStatBranchArray &operator[](int i) override { return *array[i]; }

    //! clear all suff stats
    void Clear() {
        for (int i = 0; i < GetSize(); i++) {
            array[i]->Clear();
        }
    }

  private:
    const Tree& tree;
    vector<dSOmegaPathSuffStatBranchArray *> array;
};
