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
    PoissonSuffStatBranchArray bl_suffstats{*data.tree};
    auto nuc_rates = make_nucleotide_rate(
        normalize({1, 1, 1, 1, 1, 1}), 1. / 6, normalize({1, 1, 1, 1}), 1. / 4, gen);
    MGOmegaCodonSubMatrix codon_sub_matrix(
        dynamic_cast<const CodonStateSpace*>(data.alignment.GetStateSpace()),
        &get<nuc_matrix>(nuc_rates), get<omega, value>(global_omega));
    LegacyArrayProxy branch_adapter(get<bl_array, value>(branch_lengths), *data.tree);
    PhyloProcess phyloprocess(
        data.tree.get(), &data.alignment, &branch_adapter, 0, &codon_sub_matrix);
    phyloprocess.Unfold();

    // initializing components
    ChainDriver chain_driver{cmd.chain_name(), args.every.getValue(), args.until.getValue()};
    auto scheduler = make_move_scheduler([&gen, &global_omega]() {  //
        globom::move(global_omega, []() { return 0.; }, gen);
    });

    // auto model =
    //     make_singleomega_model(globom, args.alignment.getValue(), args.treefile.getValue());
    // model->Update();
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
