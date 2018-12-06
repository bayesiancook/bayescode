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
};

#endif
