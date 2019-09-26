#pragma once

#include "lib/CodonSuffStat.hpp"

template <class T>
class SuffstatInterface {
    virtual T _get() = 0;

  protected:
    ~SuffstatInterface() = default;

  public:
    // protected non-virtual destructor, as this interface is not
    // meant to be sued with owning pointers
    virtual void gather() = 0;

    T get() {
#ifndef NDEBUG
        auto tmp = _get();
        gather();
        assert(_get() == tmp);
#endif
        return _get();
    }
};

class PathSSW final : public SuffstatInterface<PathSuffStat&> {
    PathSuffStat _ss;
    PhyloProcess& _phyloprocess;

    PathSuffStat& _get() final { return _ss; }

  public:
    PathSSW(PhyloProcess& phyloprocess) : _phyloprocess(phyloprocess) {}

    void gather() final {
        _ss.Clear();
        _ss.AddSuffStat(_phyloprocess);
    }
};

struct omega_suffstat_t {
    int count;
    double beta;
    bool operator==(const omega_suffstat_t& other) const {
        return count == other.count && beta == other.beta;
    }
};

class OmegaSSW final : public SuffstatInterface<omega_suffstat_t> {  // SSW = suff stat wrapper
    const OmegaCodonSubMatrix& _codon_submatrix;
    SuffstatInterface<PathSuffStat&>& _path_suffstat;
    OmegaPathSuffStat _ss;

    omega_suffstat_t _get() final { return {_ss.GetCount(), _ss.GetBeta()}; }

  public:
    OmegaSSW(
        const OmegaCodonSubMatrix& codon_submatrix, SuffstatInterface<PathSuffStat&>& pathsuffstat)
        : _codon_submatrix(codon_submatrix), _path_suffstat(pathsuffstat) {}

    void gather() final {
        _ss.Clear();
        _ss.AddSuffStat(_codon_submatrix, _path_suffstat.get());
    }
};