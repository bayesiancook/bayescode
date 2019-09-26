#pragma once

#include "lib/CodonSuffStat.hpp"

template <class T>
class SuffstatInterface {
    virtual T _get() = 0;

  public:
    virtual void gather() = 0;

    T get() {
#ifndef NDEBUG
        auto tmp = _get();
        gather();
        assert(_get() == tmp);
#endif
        return _get();
    }
    // no virtual destructor because not meant to be used for owning pointers
};

struct omega_suffstat_t {
    int count;
    double beta;
    bool operator==(const omega_suffstat_t& other) const {
        return count == other.count && beta == other.beta;
    }
};

class OmegaSSW : public SuffstatInterface<omega_suffstat_t> {  // SSW = suff stat wrapper
    const OmegaCodonSubMatrix& _codon_submatrix;
    const PathSuffStat& _path_suffstat;
    OmegaPathSuffStat _ss;

    omega_suffstat_t _get() final { return {_ss.GetCount(), _ss.GetBeta()}; }

  public:
    OmegaSSW(const OmegaCodonSubMatrix& codon_submatrix, const PathSuffStat& pathsuffstat)
        : _codon_submatrix(codon_submatrix), _path_suffstat(pathsuffstat) {}

    void gather() final {
        _ss.Clear();
        _ss.AddSuffStat(_codon_submatrix, _path_suffstat);
    }
};