#pragma once

#include "lib/CodonSuffStat.hpp"

template <class T = double>
class SuffstatWrapper {
    std::function<void()> _gather;
    std::function<T()> _get_value;

  public:
    template <class Gather, class GetValue>
    SuffstatWrapper(Gather gather, GetValue get_value) : _gather(gather), _get_value(get_value) {}

    double get_value() {
#ifndef NDEBUG
        double tmp = _get_value();
        _gather();
        assert(_get_value() == tmp);
#endif
        return _get_value();
    }

    void gather() { _gather(); }
};


class OmegaSSW {  // SSW = suff stat wrapper
    const OmegaCodonSubMatrix& _codon_submatrix;
    const PathSuffStat& _path_suffstat;

    OmegaPathSuffStat _ss;
    SuffstatWrapper<int> _ss_count;
    SuffstatWrapper<double> _ss_beta;

  public:
    void gather() {
        _ss.Clear();
        _ss.AddSuffStat(_codon_submatrix, _path_suffstat);
    }

    OmegaSSW(const OmegaCodonSubMatrix& codon_submatrix, const PathSuffStat& pathsuffstat)
        : _codon_submatrix(codon_submatrix),
          _path_suffstat(pathsuffstat),
          _ss_count([this]() { gather(); }, [this]() { return _ss.GetCount(); }),
          _ss_beta([this]() { gather(); }, [this]() { return _ss.GetBeta(); }) {}

    OmegaSSW(const OmegaSSW& other)
        : _codon_submatrix(other._codon_submatrix),
          _path_suffstat(other._path_suffstat),
          _ss(other._ss),
          _ss_count([this]() { gather(); }, [this]() { return _ss.GetCount(); }),
          _ss_beta([this]() { gather(); }, [this]() { return _ss.GetBeta(); }) {}

    auto& count_ssw() { return _ss_count; }
    auto& beta_ssw() { return _ss_beta; }
};