#include <cmath>
#include <fstream>
#include "SingleOmegaModel.hpp"
#include "bayes_toolbox/src/operations/logprob.hpp"
#include "bayes_toolbox/src/structure/ValueView.hpp"
#include "components/ChainCheckpoint.hpp"
#include "components/ChainDriver.hpp"
#include "components/ConsoleLogger.hpp"
#include "components/InferenceAppArgParse.hpp"
#include "components/StandardTracer.hpp"
#include "components/restart_check.hpp"
#include "global_omega.hpp"
#include "branch_array.hpp"

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

    // random generator
    auto gen = make_generator();

    // model
    auto global_omega = make_fixed_globom(1.0, 1.0, gen);
    auto branch_lengths = make_branchlength_array(parser, 0.1, 1.0);

    // initializing components
    ChainDriver chain_driver{cmd.chain_name(), args.every.getValue(), args.until.getValue()};
    // auto model =
    //     make_singleomega_model(globom, args.alignment.getValue(), args.treefile.getValue());
    // model->Update();
    ConsoleLogger console_logger;
    // ChainCheckpoint chain_checkpoint(cmd.chain_name() + ".param", *chain_driver, *model);
    // StandardTracer trace(*model, cmd.chain_name());

    // registering components to chain driver
    // chain_driver->add(*model);
    chain_driver.add(console_logger);
    // chain_driver->add(chain_checkpoint);
    // chain_driver->add(trace);

    // launching chain!
    chain_driver.go();
}
