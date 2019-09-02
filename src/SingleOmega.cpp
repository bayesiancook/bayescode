#include <cmath>
#include <fstream>
#include "SingleOmegaModel.hpp"
#include "components/ChainCheckpoint.hpp"
#include "components/ChainDriver.hpp"
#include "components/ConsoleLogger.hpp"
#include "components/InferenceAppArgParse.hpp"
#include "components/StandardTracer.hpp"
#include "components/restart_check.hpp"

using namespace std;

int main(int argc, char *argv[]) {
    // parsing command-line arguments
    ChainCmdLine cmd{argc, argv, "SingleOmega", ' ', "0.1"};
    InferenceAppArgParse args(cmd);
    cmd.parse();

    // initializing components
    auto chain_driver =
        make_unique<ChainDriver>(cmd.chain_name(), args.every.getValue(), args.until.getValue());
    auto model = make_unique<SingleOmegaModel>(args.alignment.getValue(), args.treefile.getValue());
    model->Update();
    ConsoleLogger console_logger;
    ChainCheckpoint chain_checkpoint(cmd.chain_name() + ".param", *chain_driver, *model);
    StandardTracer trace(*model, cmd.chain_name());

    // registering components to chain driver
    chain_driver->add(*model);
    chain_driver->add(console_logger);
    chain_driver->add(chain_checkpoint);
    chain_driver->add(trace);

    // launching chain!
    chain_driver->go();
}
