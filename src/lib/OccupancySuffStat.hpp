
#ifndef OCCSUFFSTAT_H
#define OCCSUFFSTAT_H

#include "Array.hpp"
#include "SuffStat.hpp"

/**
 * \brief A sufficient statistic for a multinomial mixture allocation vector
 * (see MultinomialAllocationVector)
 *
 * Given a mixture of K components, and a vector of N integers,
 * specifying the allocation of N items to the component of the mixtures,
 * OccupancySuffStat simply stores the counts (numbers of sites allocated to
 * each component), which is a sufficient statistic for the underlying component
 * weights.
 *
 */

class OccupancySuffStat : public SimpleArray<int>, public SuffStat {
  public:
    //! \brief constructor (parameterized by mixture size)
    OccupancySuffStat(int insize) : SimpleArray<int>(insize, 0) {}
    ~OccupancySuffStat() {}

    //! \brief return number of non-empty components
    int GetNcluster() const {
        int n = 0;
        for (int i = 0; i < GetSize(); i++) {
            if (GetVal(i)) { n++; }
        }
        return n;
    }

    //! reset count vector
    void Clear() {
        for (int i = 0; i < GetSize(); i++) { (*this)[i] = 0; }
    }

    //! implement additive behavior of OccupancySuffStat
    void Add(const OccupancySuffStat &from) {
        if (from.GetSize() != GetSize()) {
            cerr << "error in OccupancySuffStat::Add: non matching array size\n";
            exit(1);
        }
        for (int i = 0; i < GetSize(); i++) { (*this)[i] += from.GetVal(i); }
    }

    //! implement additive behavior of OccupancySuffStat
    OccupancySuffStat &operator+=(const OccupancySuffStat &from) {
        Add(from);
        return *this;
    }

    //! increment count for component i
    void Increment(int i) { (*this)[i]++; }

    //! add suff stat based on an allocation vector
    void AddSuffStat(const Selector<int> &alloc) {
        for (int i = 0; i < alloc.GetSize(); i++) { (*this)[alloc.GetVal(i)]++; }
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

    //! get an OccupancySuffStat from MPI buffer and then add it to this object
    void Add(const MPIBuffer &buffer) {
        for (int i = 0; i < GetSize(); i++) {
            int temp;
            buffer >> temp;
            (*this)[i] += temp;
        }
    }
};

#endif
