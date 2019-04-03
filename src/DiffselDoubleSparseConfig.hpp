#pragma once

#include "components/param_enums.hpp"
#include "global/logging.hpp"
using std::istream;
using std::ostream;
using std::string;

struct DiffselDoubleSparseConfig {
    string datafile{""};
    string treefile{""};

    int Ncond{0};
    int Nlevel{0};
    int codonmodel{0};

    double pihypermean{0.};
    double shiftprobmean{0.};
    double shiftprobinvconc{0.};
    double maskepsilon{0.};
    double fitnessshape{0.};

    // branch lengths
    param_mode_t blmode;  // branch lengths fixed or sampled
    double bl_multiplier{0.};

    // fitness-related parameters
    param_mode_t fitnessshapemode;   // estimation method for fitness hyperparameter
    double fixed_shape{0.};          // (shape of multi-gamma distribution)
    param_mode_t fitnesscentermode;  // estimation method for fitness hyperparameter (center of
                                     // multi-gamma distribution)

    param_mode_t nucmode;  // mutation matrix parameters fixed or sampled


    mask_mode_t maskmode;  // estimation method for site profile masks. Used in a multigene context.

    // epsilon
    int maskepsilonmode;
    double epsilon{0.};

    // sparse model variants
    bool withtoggle;
    bool site_wise;

    void init() {
        if (shape > 0) {
            fitnessshapemode = param_mode_t::fixed;
            fitnessshape = shape;
        } else {
            fitnessshapemode = independent;
            fitnessshape = 20.0;
        }

        if (epsilon == 1) {
            maskepsilon = 1;
            maskmode = no_mask;
            maskepsilonmode = 3;
        } else if (epsilon >= 0) {
            maskepsilon = epsilon;
            maskepsilonmode = 3;
            maskmode = gene_spec_mask_fixed_hyper;
        } else {
            maskepsilonmode = 0;
            maskmode = gene_spec_mask_fixed_hyper;
            maskepsilon = 0.01;
        }
    }
};

ostream& operator<<(ostream& os, DiffselDoubleSparseConfig& config) {
    os << fmt::format(
        "datafile: {}\n\ttreefile: {}\n\tfitnesscentermode: {}\n\twithtoggle: {}\n\tpihypermean: "
        "{}\n\tshiftprobmean: {}\n\tshiftprobinvconc: {}\n\tcodonmodel: {}\n\tblmode: "
        "{}\n\tnucmode: {}\n\tfitnessshapemode: {}\n\tfitnessshape: {}\n\tmaskepsilon: "
        "{}\n\tmaskmode: {}\n\tmaskepsilonmode: {}\n\tNcond: {}\n\tNlevel: {}\n\tsite_wise: {}",
        config.datafile, config.treefile, config.fitnesscentermode, config.withtoggle,
        config.pihypermean, config.shiftprobmean, config.shiftprobinvconc, config.codonmodel,
        config.blmode, config.nucmode, config.fitnessshapemode, config.fitnessshape,
        config.maskepsilon, config.maskmode, config.maskepsilonmode, config.Ncond, config.Nlevel,
        config.site_wise);
    return os;
}