#pragma once

#include "Proxy.hpp"
#include "lib/CodonSuffStat.hpp"
#include "lib/GTRSubMatrix.hpp"

// =================================================================================================
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

// =================================================================================================
class NucMatrixProxy : public Proxy<GTRSubMatrix&> {
    GTRSubMatrix& _mat;
    const std::vector<double>& _eq_freqs;

    GTRSubMatrix& _get() final { return _mat; }

  public:
    NucMatrixProxy(GTRSubMatrix& mat, const std::vector<double>& eq_freqs)
        : _mat(mat), _eq_freqs(eq_freqs) {}

    void gather() final {
        _mat.CopyStationary(_eq_freqs);
        _mat.CorruptMatrix();
    }
};

// =================================================================================================
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

// =================================================================================================
struct brancharray_poisson_suffstat_t {
    int count;
    double beta;
    bool operator==(const brancharray_poisson_suffstat_t& other) const {
        return count == other.count && beta == other.beta;
    }
};

class BranchArrayPoissonSSW final : public Proxy<brancharray_poisson_suffstat_t, int> {
    PoissonSuffStatBranchArray _ss;
    PhyloProcess& _phyloprocess;

    brancharray_poisson_suffstat_t _get(int i) final {
        auto& local_ss = _ss.GetVal(i);
        return {local_ss.GetCount(), local_ss.GetBeta()};
    }

  public:
    BranchArrayPoissonSSW(const Tree& tree, PhyloProcess& phyloprocess)
        : _ss(tree), _phyloprocess(phyloprocess) {}

    void gather() final {
        _ss.Clear();
        _ss.AddLengthPathSuffStat(_phyloprocess);
    }
};