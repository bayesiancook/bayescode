
#ifndef BRANCHARRAY_H
#define BRANCHARRAY_H

#include <vector>
#include "MPIBuffer.hpp"
#include "Tree.hpp"

/**
 * \brief An interface for an abstract tree-structured array (indexed by branches) of constant references over objects of
 * type T
 *
 * This abstract class is meant as a general read-only interface
 * returning a reference over a T object for any branch of a given phylogenetic tree (through the GetVal(int index) method).
 * Branches are indexed by integers (ranging from 0 to GetNbranch()-1).
 * This class, and its various implementations/specializations, follows the same logic as the classes deriving from Array,
 * except for the fact that the Branch selectors are meant to define distributions across branches.
 *
 */

template <class T>
class BranchSelector {
  public:
    virtual ~BranchSelector() {}

    //! return the number of branches of the underlying tree
    int GetNbranch() const { return GetTree().GetNbranch(); }

    //! return a const reference to the underlying tree
    virtual const Tree& GetTree() const = 0;
    //! const access to array element by index
    virtual const T& GetVal(int index) const = 0;

    //! return size of entire array, when put into an MPI buffer
    unsigned int GetMPISize() const { return this->GetNbranch() * MPISize(this->GetVal(0)); }

    //! write array into MPI buffer
    void MPIPut(MPIBuffer& buffer) const {
        for (int i = 0; i < this->GetNbranch(); i++) {
            buffer << this->GetVal(i);
        }
    }

    //! write array into generic output stream
    void ToStream(ostream& os) const {
        for (int i = 0; i < this->GetNbranch(); i++) {
            os << this->GetVal(i) << '\t';
        }
        os << '\n';
    }
};

/**
 * \brief An interface for an abstract tree-structured array (indexed by branches) of references over objects of type T
 *
 * Deriving from BranchSelector, BranchArray gives const access through GetVal(int index).
 * In addition, it provides a non-const access to through the [] operator.
 * In practice, this class is (should be?) used only for implementing arrays of all-distinct instances of objects of type T.
 * See also Array (and Selector) for a similar class hierarchy.
 */

template <class T>
class BranchArray : public BranchSelector<T> {
  public:
    virtual ~BranchArray() {}

    //! element-by-element copy (arrays should be of same size)
    void Copy(const BranchSelector<T>& from) {
        if (this->GetNbranch() != from.GetNbranch()) {
            cerr << "error: branch arrays do not have same size\n";
            exit(1);
        }
        for (int i = 0; i < this->GetNbranch(); i++) {
            (*this)[i] = from.GetVal(i);
        }
    }

    //! non-const access to array element by index
    virtual T& operator[](int index) = 0;

    //! get array from MPI buffer
    void MPIGet(const MPIBuffer& buffer) {
        for (int i = 0; i < this->GetNbranch(); i++) {
            buffer >> (*this)[i];
        }
    }

    //! get array from generic input stream
    void FromStream(istream& is) {
        for (int i = 0; i < this->GetNbranch(); i++) {
            is >> (*this)[i];
        }
    }
};

/**
 * \brief output stream operator for a BranchSelector<T>
 *
 * in practice, calls the ToStream method of BranchSelector<T>
 */

template <class T>
ostream& operator<<(ostream& os, const BranchSelector<T>& array) {
    array.ToStream(os);
    return os;
}

/**
 * \brief input stream operator for a BranchArray<T>
 *
 * in practice, calls the FromStream method of BranchArray<T>
 */

template <class T>
istream& operator>>(istream& is, BranchArray<T>& array) {
    array.FromStream(is);
    return is;
}

/**
 * \brief A BranchSelector<T> that returns a reference to the same value T for any branch
 *
 * Useful for implementing simple models assuming, e.g. the same substitution matrix for all branches.
 */

template <class T>
class BranchHomogeneousSelector : public BranchSelector<T> {
  public:
    //! \brief Constructor, taking as its arguments the tree and the value to be returned for any branch
    BranchHomogeneousSelector(const Tree* intree, const T& invalue) : tree(intree), value(invalue) {}
    ~BranchHomogeneousSelector() {}

    const Tree& GetTree() const /*override*/ { return tree; }

    //! return a reference to the same value (i.e. value) for any index
    const T& GetVal(int index) const /*override*/ { return value; }

  private:
    const Tree& tree;
    const T& value;
};

/**
 * \brief The 'standard' implementation of an BranchArray<T>, simply as a std::vector<T>
 *
 */

template <class T>
class SimpleBranchArray : public BranchArray<T> {
  public:
    //! Constructor (with only the tree given as argument)
    SimpleBranchArray(const Tree* intree) : tree(intree), array(intree->GetNbranch()) {}

    //! Constructor with tree and initializer value
    SimpleBranchArray(const Tree* intree, const T& initval) : tree(intree), array(intree->GetNbranch(), initval) {}
    virtual ~SimpleBranchArray() {}

    const Tree& GetTree() const /*override*/ { return *tree; }
    T& operator[](int index) /*override*/ { return array[index]; }
    const T& GetVal(int index) const /*override*/ { return array[index]; }

  private:
    const Tree* tree;
    vector<T> array;
};

#endif
