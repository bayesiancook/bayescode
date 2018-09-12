#ifndef PERM_H
#define PERM_H

#include "Array.hpp"

/**
 * \brief A permutation of 1..N
 */

class Permutation : public SimpleArray<int> {
  public:
    //! constructor (parameterized by N)
    Permutation(int size) : SimpleArray<int>(size) { Reset(); }

    ~Permutation() {}

    //! set equal to identity permuation
    void Reset() {
        for (int i = 0; i < GetSize(); i++) { (*this)[i] = i; }
    }

    //! return size when put into an MPI buffer
    unsigned int GetMPISize() const { return GetSize(); }

    //! put current value of count and beta into an MPI buffer
    void MPIPut(MPIBuffer &buffer) const {
        for (int i = 0; i < GetSize(); i++) { buffer << GetVal(i); }
    }

    //! get value from MPI buffer
    void MPIGet(const MPIBuffer &buffer) {
        for (int i = 0; i < GetSize(); i++) { buffer >> (*this)[i]; }
    }
};

#endif
