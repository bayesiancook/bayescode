#pragma once

#include <ostream>

//! - mode == 3: global == "fixed"
//! - mode == 2: global but estimated == "shared"
//! - mode == 1: gene specific, with hyperparameters estimated across genes == "shrunk"
//! - mode == 0: gene-specific, with fixed hyperparameters == "independent"
enum param_mode_t { independent, shrunk, shared, fixed };

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

bool resampled(param_mode_t p) { return p == independent || p == shrunk; }

//! Used in a multigene context.
//! - mode == 3: no mask "no_mask"
//! - mode == 2: parameter (maskprob) shared across genes "shared_mask"
//! - mode == 1: gene-specific parameter (maskprob), hyperparameters estimated
//! across genes "gene_spec_mask_est_hyper"
//! - mode == 0: gene-specific parameter (maskprob) with fixed hyperparameters
//! "gene_spec_mask_fixed_hyper"
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