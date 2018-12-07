
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
            std::cerr << "error in OccupancySuffStat::Add: non matching array size\n";
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
};

#endif
