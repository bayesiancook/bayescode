
#ifndef OCCSUFFSTAT_H
#define OCCSUFFSTAT_H

#include "Array.hpp"
#include "SuffStat.hpp"

class OccupancySuffStat : public SimpleArray<int>, public SuffStat {
  public:
    OccupancySuffStat(int insize) : SimpleArray<int>(insize, 0) {}
    ~OccupancySuffStat() {}

    int GetNcluster() const {
        int n = 0;
        for (int i = 0; i < GetSize(); i++) {
            if (GetVal(i)) {
                n++;
            }
        }
        return n;
    }

    void Clear() {
        for (int i = 0; i < GetSize(); i++) {
            (*this)[i] = 0;
        }
    }

    void Add(const OccupancySuffStat& from) {
        if (from.GetSize() != GetSize()) {
            cerr << "error in OccupancySuffStat::Add: non matching array size\n";
            exit(1);
        }
        for (int i = 0; i < GetSize(); i++) {
            (*this)[i] += from.GetVal(i);
        }
    }

    OccupancySuffStat& operator+=(const OccupancySuffStat& from) {
        Add(from);
        return *this;
    }

    void Increment(int i) { (*this)[i]++; }

    void AddSuffStat(const Selector<int>& alloc) {
        for (int i = 0; i < alloc.GetSize(); i++) {
            (*this)[alloc.GetVal(i)]++;
        }
    }

    unsigned int GetMPISize() const { return GetSize(); }

    void MPIPut(MPIBuffer& buffer) const {
        for (int i = 0; i < GetSize(); i++) {
            buffer << GetVal(i);
        }
    }

    void MPIGet(const MPIBuffer& buffer) {
        for (int i = 0; i < GetSize(); i++) {
            buffer >> (*this)[i];
        }
    }

    void Add(const MPIBuffer& buffer) {
        for (int i = 0; i < GetSize(); i++) {
            int temp;
            buffer >> temp;
            (*this)[i] += temp;
        }
    }
};

#endif
