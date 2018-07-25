
#ifndef BIDIMARRAY_H
#define BIDIMARRAY_H

#include <vector>
#include "MPIBuffer.hpp"

/**
 * \brief An interface for an abstract indexed bidim-array of constant
 * references over objects of type T
 *
 * This abstract class is meant as a general read-only interface
 * returning a reference over a T object for any pair of indices (through the
 * GetVal(int i, int j) method). As an example, DiffSelSparseModel implements a
 * bidim array of site- and condition-specific fitness vectors over the 20 amino
 * acids (fitness vector)
 */

template <class T>
class BidimSelector {
  public:
    virtual ~BidimSelector() {}

    //! return number of rows
    virtual int GetNrow() const = 0;
    //! return number of columns
    virtual int GetNcol() const = 0;
    //! const access to element array
    virtual const T &GetVal(int i, int j) const = 0;

    //! return size of entire array, when put into an MPI buffer
    unsigned int GetMPISize() const {
        return this->GetNrow() * this->GetNcol() * MPISize(this->GetVal(0, 0));
    }

    //! write array into MPI buffer
    void MPIPut(MPIBuffer &buffer) const {
        for (int i = 0; i < this->GetNrow(); i++) {
            for (int j = 0; j < this->GetNcol(); j++) {
                buffer << this->GetVal(i, j);
            }
        }
    }

    //! write array into generic output stream
    void ToStream(ostream &os) const {
        for (int i = 0; i < this->GetNrow(); i++) {
            for (int j = 0; j < this->GetNcol(); j++) {
                os << this->GetVal(i, j) << '\t';
            }
        }
    }
};

/**
 * \brief An interface for an abstract indexed array of references over objects
 * of type T
 *
 * Deriving from BidimSelector, BidimArray gives const access through GetVal(int
 * i, int j). In addition, it provides a non-const access to through the ()
 * operator. In practice, this class is (should be?) used only for implementing
 * arrays of all-distinct instances of objects of type T. The SimpleBidimArray
 * class template provides a more explicit implementation (by encapsulating a
 * vector<vector<T> >, see SimpleBidimArray). However, in some cases, it can be
 * useful to implement a BidimArray<T> as a vector<vector<T*> > instead, hence
 * the specific declaration of this intermediate class.
 */

template <class T>
class BidimArray : public BidimSelector<T> {
  public:
    virtual ~BidimArray() {}

    //! non-const access to array element
    virtual T &operator()(int i, int j) = 0;

    //! get array from MPI buffer
    void MPIGet(const MPIBuffer &buffer) {
        for (int i = 0; i < this->GetNrow(); i++) {
            for (int j = 0; j < this->GetNcol(); j++) {
                buffer >> (*this)(i, j);
            }
        }
    }

    //! get array from generic input stream
    void FromStream(istream &is) {
        for (int i = 0; i < this->GetNrow(); i++) {
            for (int j = 0; j < this->GetNcol(); j++) {
                is >> (*this)(i, j);
            }
        }
    }
};

/**
 * \brief output stream operator for a BidimSelector<T>
 *
 * in practice, calls the ToStream method of BidimSelector<T>
 */

template <class T>
ostream &operator<<(ostream &os, const BidimSelector<T> &array) {
    array.ToStream(os);
    return os;
}

/**
 * \brief input stream operator for a BidimArray<T>
 *
 * in practice, calls the FromStream method of BidimArray<T>
 */

template <class T>
istream &operator>>(istream &is, BidimArray<T> &array) {
    array.FromStream(is);
    return is;
}

/**
 * \brief A BidimSelector<T> that returns a reference to the same value T for
 * any pair of indices (i,j)
 */

template <class T>
class BidimHomogeneousSelector : public BidimSelector<T> {
  public:
    //! \brief Constructor, taking as its arguments the number of rows and columns
    //! of the array and the value to be returned for any pair of indices
    BidimHomogeneousSelector(int innrow, int inncol, const T &invalue)
        : nrow(innrow), ncol(inncol), value(invalue) {}
    ~BidimHomogeneousSelector() {}

    int GetNrow() const override { return nrow; }
    int GetNcol() const override { return ncol; }

    //! return a reference to the same value (i.e. value) for any pair of indices
    const T &GetVal(int i, int j) const override { return value; }

  private:
    int nrow;
    int ncol;
    const T &value;
};

/**
 * \brief The 'standard' implementation of a BidimArray<T>, simply as a
 * vector<vector<T> >
 *
 */

template <class T>
class SimpleBidimArray : public BidimArray<T> {
  public:
    //! Constructor with array dimensions and initializer value
    SimpleBidimArray(int innrow, int inncol, const T &initval)
        : nrow(innrow), ncol(inncol), array(innrow, vector<T>(inncol, initval)) {}
    virtual ~SimpleBidimArray() {}

    int GetNrow() const override { return nrow; }
    int GetNcol() const override { return ncol; }

    T &operator()(int i, int j) override { return array[i][j]; }
    const T &GetVal(int i, int j) const override { return array[i][j]; }

    //! return a const ref to the vector<T> corresponding to row i
    const vector<T> &GetSubArray(int i) const { return array[i]; }

  protected:
    int nrow;
    int ncol;
    vector<vector<T>> array;
};

#endif
