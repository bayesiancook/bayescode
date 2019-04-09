#pragma once

#include "model/param_classes.hpp"

struct DiffselDoubleSparseConfig {
    string datafile{""};
    string treefile{""};

    // number of diff sel categories
    int Ncond{0};

    // number of levels of the model
    // with 2 levels, structure of the model is as follows:
    // baseline (condition 0)
    // baseline  || fitness1 (for condition 1)
    // baseline || fitness1  || fitnessk  (for condition k=2..Ncond)
    int Nlevel{0};

    MultiGeneParameter<double, hyper_mean_invshape> branch_lengths;
    MultiGeneParameter<double, hyper_rate> fitness_shape;  // fitness is a multi-gamma
    MultiGeneParameter<per_aa<double>, no_hyper> fitness_center;
    MultiGeneParameter<double, hyper_mean_invconc> mask_prob;
    MultiGeneParameter<double, no_hyper> epsilon;  // epsilon is uniform

    param_mode_t nucmode;  // mutation matrix parameters fixed or sampled

    // sparse model variants
    bool withtoggle;
    bool site_wise;

    // pi-related
    double pihypermean{0.};

    void check() const {
        assert(Ncond > 0);
        assert(Nlevel > 0);
        assert(branch_lengths.check());
        assert(fitness_shape.check());
        assert(fitness_center.check());
        assert(mask_prob.check());
        assert(epsilon.check());
    }
};

ostream& operator<<(ostream& os, DiffselDoubleSparseConfig& config) {
    // os << fmt::format(
    //     "datafile: {}\n\ttreefile: {}\n\tfitnesscentermode: {}\n\twithtoggle: {}\n\tpihypermean:
    //     "
    //     "{}\n\tshiftprobmean: {}\n\tshiftprobinvconc: {}\n\tcodonmodel: {}\n\tblmode: "
    //     "{}\n\tnucmode: {}\n\tfitnessshapemode: {}\n\tfitnessshape: {}\n\tmaskepsilon: "
    //     "{}\n\tmaskmode: {}\n\tmaskepsilonmode: {}\n\tNcond: {}\n\tNlevel: {}\n\tsite_wise: {}",
    //     config.datafile, config.treefile, config.fitnesscentermode, config.withtoggle,
    //     config.pihypermean, config.shiftprobmean, config.shiftprobinvconc, config.codonmodel,
    //     config.blmode, config.nucmode, config.fitnessshapemode, config.fitnessshape,
    //     config.maskepsilon, config.maskmode, config.maskepsilonmode, config.Ncond, config.Nlevel,
    //     config.site_wise);
    return os;
}