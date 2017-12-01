
#ifndef OCCSUFFSTAT_H
#define OCCSUFFSTAT_H

#include "SuffStat.hpp"
#include "Array.hpp"

class OccupancySuffStat : public SimpleArray<int>, public SuffStat   {

    public:

    OccupancySuffStat(int insize) : SimpleArray<int>(insize,0) {}
    ~OccupancySuffStat() {}

    void Clear()    {
        for (int i=0; i<GetSize(); i++) {
            (*this)[i] = 0;
        }
    }

    void Add(const OccupancySuffStat& from) {
        if (from.GetSize() != GetSize())    {
            cerr << "error in OccupancySuffStat::Add: non matching array size\n";
            exit(1);
        }
        for (int i=0; i<GetSize(); i++) {
            (*this)[i] += from.GetVal(i);
        }
    }

    OccupancySuffStat& operator+=(const OccupancySuffStat& from)    {
        Add(from);
        return *this;
    }

    void Increment(int i)   {
        (*this)[i]++;
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

    void Add(const MPIBuffer& buffer)   {
        for (int i=0; i<GetSize(); i++) {
            int temp;
            buffer >> temp;
            (*this)[i] += temp;
        }
    }
};

#endif

