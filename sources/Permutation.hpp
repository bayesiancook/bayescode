#ifndef PERM_H
#define PERM_H

#include "Array.hpp"

class Permutation : public SimpleArray<int> {

    public:

    Permutation(int size) : SimpleArray<int>(size) {
        Reset();
    }

    ~Permutation() {}

    void Reset()    {
        for (int i=0; i<GetSize(); i++) {
            (*this)[i] = i;
        }
    }

    unsigned int GetMPISize() const {return GetSize();}

    void MPIPut(MPIBuffer& buffer) const    {
        for (int i=0; i<GetSize(); i++) {
            buffer << GetVal(i);
        }
    }

    void MPIGet(const MPIBuffer& buffer)    {
        for (int i=0; i<GetSize(); i++) {
            buffer >> (*this)[i];
        }
    }
};

#endif
