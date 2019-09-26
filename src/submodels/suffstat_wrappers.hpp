#pragma once

#include "Proxy.hpp"
#include "lib/CodonSuffStat.hpp"

class PathSSW final : public Proxy<PathSuffStat&> {
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

class OmegaSSW final : public Proxy<omega_suffstat_t> {  // SSW = suff stat wrapper
    const OmegaCodonSubMatrix& _codon_submatrix;
    Proxy<PathSuffStat&>& _path_suffstat;
    OmegaPathSuffStat _ss;

    omega_suffstat_t _get() final { return {_ss.GetCount(), _ss.GetBeta()}; }

  public:
    OmegaSSW(const OmegaCodonSubMatrix& codon_submatrix, Proxy<PathSuffStat&>& pathsuffstat)
        : _codon_submatrix(codon_submatrix), _path_suffstat(pathsuffstat) {}

    void gather() final {
        _ss.Clear();
        _ss.AddSuffStat(_codon_submatrix, _path_suffstat.get());
    }
};