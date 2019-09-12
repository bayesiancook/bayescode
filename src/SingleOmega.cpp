#include <cmath>
#include "components/ChainCheckpoint.hpp"
#include "components/ChainDriver.hpp"
#include "components/ConsoleLogger.hpp"
#include "components/InferenceAppArgParse.hpp"
#include "components/MoveScheduler.hpp"
#include "components/StandardTracer.hpp"
#include "components/restart_check.hpp"
#include "data_preparation.hpp"
#include "lib/CodonSubMatrix.hpp"
#include "lib/CodonSuffStat.hpp"
#include "lib/PoissonSuffStat.hpp"
#include "submodels/branch_array.hpp"
#include "submodels/global_omega.hpp"
#include "submodels/nuc_rates.hpp"


using namespace std;


class LegacyArrayProxy : public BranchSelector<double> {
    std::vector<double>& data_ref;
    const Tree& tree_ref;

  public:
    LegacyArrayProxy(std::vector<double>& data_ref, const Tree& tree_ref)
        : data_ref(data_ref), tree_ref(tree_ref) {}

    virtual const Tree& GetTree() const override { return tree_ref; }
    virtual const double& GetVal(int index) const override { return data_ref[index]; }
};


int main(int argc, char* argv[]) {
    // parsing command-line arguments
    ChainCmdLine cmd{argc, argv, "SingleOmega", ' ', "0.1"};
    InferenceAppArgParse args(cmd);
    cmd.parse();

    // input data
    auto data = prepare_data(args.alignment.getValue(), args.treefile.getValue());

    // random generator
    auto gen = make_generator();

    // model
    auto global_omega = globom::make_fixed(1.0, 1.0, gen);
    auto branch_lengths = make_branchlength_array(data.parser, 0.1, 1.0);
    auto nuc_rates = make_nucleotide_rate(
        normalize({1, 1, 1, 1, 1, 1}), 1. / 6, normalize({1, 1, 1, 1}), 1. / 4, gen);
    MGOmegaCodonSubMatrix codon_sub_matrix(
        dynamic_cast<const CodonStateSpace*>(data.alignment.GetStateSpace()),
        &get<nuc_matrix>(nuc_rates), get<omega, value>(global_omega));
    LegacyArrayProxy branch_adapter(get<bl_array, value>(branch_lengths), *data.tree);
    PhyloProcess phyloprocess(
        data.tree.get(), &data.alignment, &branch_adapter, 0, &codon_sub_matrix);
    phyloprocess.Unfold();

    // suff stats
    PoissonSuffStatBranchArray bl_suffstats{*data.tree};
    PathSuffStat path_suffstats;

    // move schedule
    auto touch_matrices = [&global_omega, &codon_sub_matrix, &nuc_rates]() {
        get<nuc_matrix>(nuc_rates).CopyStationary(get<eq_freq, value>(nuc_rates));
        get<nuc_matrix>(nuc_rates).CorruptMatrix();
        codon_sub_matrix.SetOmega(get<omega, value>(global_omega));
        codon_sub_matrix.CorruptMatrix();
    };

    auto scheduler = make_move_scheduler(
        [&gen, &global_omega, &phyloprocess, &touch_matrices, &path_suffstats]() {  //
            // move phyloprocess
            touch_matrices();
            phyloprocess.Move(1.0);

            // move omega
            for (int rep = 0; rep < 10; rep++) {
                // move omega
                path_suffstats.Clear();
                path_suffstats.AddSuffStat(phyloprocess);
                globom::move(global_omega, []() { return 0.; }, gen);  //@fixme with real logprob

                // move nuc rates
                touch_matrices();
                // CollectNucPathSuffStat();

                // Move::Profile(nucrelrate, 0.1, 1, 3, &SingleOmegaModel::NucRatesLogProb,
                //     &SingleOmegaModel::TouchNucMatrix, this);
                // Move::Profile(nucrelrate, 0.03, 3, 3, &SingleOmegaModel::NucRatesLogProb,
                //     &SingleOmegaModel::TouchNucMatrix, this);
                // Move::Profile(nucrelrate, 0.01, 3, 3, &SingleOmegaModel::NucRatesLogProb,
                //     &SingleOmegaModel::TouchNucMatrix, this);

                // Move::Profile(nucstat, 0.1, 1, 3, &SingleOmegaModel::NucRatesLogProb,
                //     &SingleOmegaModel::TouchNucMatrix, this);
                // Move::Profile(nucstat, 0.01, 1, 3, &SingleOmegaModel::NucRatesLogProb,
                //     &SingleOmegaModel::TouchNucMatrix, this);

                touch_matrices();
            }


        });

    // initializing components
    ChainDriver chain_driver{cmd.chain_name(), args.every.getValue(), args.until.getValue()};

    ConsoleLogger console_logger;
    // ChainCheckpoint chain_checkpoint(cmd.chain_name() + ".param", *chain_driver, *model);
    // StandardTracer trace(*model, cmd.chain_name());

    // registering components to chain driver
    // chain_driver->add(*model);
    chain_driver.add(scheduler);
    chain_driver.add(console_logger);
    // chain_driver->add(chain_checkpoint);
    // chain_driver->add(trace);

    // launching chain!
    chain_driver.go();
}
