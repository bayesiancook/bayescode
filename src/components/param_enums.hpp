#pragma once

#include <ostream>

// global == "fixed"
// global but estimated == "shared"
// gene specific, with hyperparameters estimated across genes == "shrunk"
// gene-specific, with fixed hyperparameters == "independent"
enum param_mode_t { independent, shrunk, shared, fixed };

bool resampled(param_mode_t p) { return p == independent || p == shrunk; }

std::ostream &operator<<(std::ostream &os, const param_mode_t &c) {
    if (c == independent) {
        os << "independent";
    } else if (c == shrunk) {
        os << "shrunk";
    } else if (c == shared) {
        os << "shared";
    } else if (c == param_mode_t::fixed) {
        os << "fixed";
    }
    return os;
}

template <class ValueType, class HyperValueType>
class MultiGeneParameter {
    // invariant: ensures value and hyper_value are present iff needed
    
    param_mode_t _mode;
    ValueType _value;             // assumes it's default-initializable
    HyperValueType _hyper_value;  // assumes it's default-initializable
  public:
    MultiGeneParameter(param_mode_t mode) : _mode(mode) {
        assert(mode == shared or mode == shrunk);
    }
    MultiGeneParameter(param_mode_t mode, ValueType value) : _mode(mode), _value(value) {
        assert(mode == fixed);
    }
    MultiGeneParameter(param_mode_t mode, HyperValueType hyper_value)
        : _mode(mode), _hyper_value(hyper_value) {
        assert(mode == shrunk);
    }
    param_mode_t mode() const { return mode; }
    ValueType value() const {
        assert(mode == fixed);
        return value;
    }
    HyperValueType hyper_value() const {
        assert(mode == shrunk);
        return hyper_value;
    }
};

enum mask_mode_t { gene_spec_mask_fixed_hyper, gene_spec_mask_est_hyper, shared_mask, no_mask };

std::ostream &operator<<(std::ostream &os, const mask_mode_t &c) {
    if (c == gene_spec_mask_fixed_hyper) {
        os << "gene_spec_mask_fixed_hyper";
    } else if (c == gene_spec_mask_est_hyper) {
        os << "gene_spec_mask_est_hyper";
    } else if (c == shared_mask) {
        os << "shared_mask";
    } else if (c == no_mask) {
        os << "no_mask";
    }
    return os;
}

bool gene_specific_mask_mode(mask_mode_t p) {
    return p == gene_spec_mask_fixed_hyper || p == gene_spec_mask_est_hyper;
}