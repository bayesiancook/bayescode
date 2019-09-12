#include <cmath>
#include <fstream>
#include "CodonSequenceAlignment.hpp"
#include "components/ChainCheckpoint.hpp"
#include "components/ChainDriver.hpp"
#include "components/ConsoleLogger.hpp"
#include "components/InferenceAppArgParse.hpp"
#include "components/MoveScheduler.hpp"
#include "components/StandardTracer.hpp"
#include "components/restart_check.hpp"
#include "submodels/branch_array.hpp"
#include "submodels/global_omega.hpp"

using namespace std;

int main(int argc, char* argv[]) {
    // parsing command-line arguments
    ChainCmdLine cmd{argc, argv, "SingleOmega", ' ', "0.1"};
    InferenceAppArgParse args(cmd);
    cmd.parse();

    // parsing tree
    std::ifstream tree_stream{args.treefile.getValue()};
    NHXParser parser{tree_stream};
    auto tree = make_from_parser(parser);
    assert(tree->nb_nodes() > 0);
    DEBUG("Parsed tree with {} nodes.", tree->nb_nodes());

    // sequence alignment
    FileSequenceAlignment nuc_align(args.alignment.getValue());
    assert(nuc_align.GetNtaxa() > 0 && nuc_align.GetNsite() > 0);
    DEBUG("Parsed alignment with {} sequences of length {}. Example taxon name: {}.",
        nuc_align.GetNtaxa(), nuc_align.GetNsite(), nuc_align.GetTaxonSet()->GetTaxon(0));

    CodonSequenceAlignment alignment(&nuc_align);
    assert(alignment.GetNtaxa() > 0 && alignment.GetNsite() > 0);
    DEBUG("Converted alignment to codons (new length: {}).", alignment.GetNsite());
    
    const TaxonSet taxon_set = *alignment.GetTaxonSet();
    DEBUG("Got a taxon set of length {}. Example taxon name: {}.", taxon_set.GetNtaxa(),
        taxon_set.GetTaxon(0));

    // random generator
    auto gen = make_generator();

    // model
    auto global_omega = globom::make_fixed(1.0, 1.0, gen);
    auto branch_lengths = make_branchlength_array(parser, 0.1, 1.0);

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
