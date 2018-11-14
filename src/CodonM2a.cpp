#include <cmath>
#include <fstream>
#include "CodonM2aModel.hpp"
#include "components/ChainCheckpoint.hpp"
#include "components/ChainDriver.hpp"
#include "components/ConsoleLogger.hpp"
#include "components/InferenceAppArgParse.hpp"
#include "components/StandardTracer.hpp"
#include "components/restart_check.hpp"

using namespace std;

int main(int argc, char *argv[]) {
    ChainCmdLine cmd{argc, argv, "CodonM2a", ' ', "0.1"};

    ChainDriver *chain_driver = nullptr;
    CodonM2aModel *model = nullptr;

    if (cmd.resume_from_checkpoint()) {
        std::ifstream is = cmd.checkpoint_file();
        chain_driver = new ChainDriver(is);
        model = new CodonM2aModel(is);
        check_restart(*model, cmd.chain_name() + ".trace");
    } else {
        InferenceAppArgParse args(cmd);

        TCLAP::ValueArg<double> pi{"p", "pi",
            "prior prob of being under positive selection (default: 0.1)", false, 0.1, "double",
            cmd.get()};

        cmd.parse();
        chain_driver =
            new ChainDriver(cmd.chain_name(), args.every.getValue(), args.until.getValue());
        model =
            new CodonM2aModel(args.alignment.getValue(), args.treefile.getValue(), pi.getValue());
    }

    ConsoleLogger console_logger;
    ChainCheckpoint chain_checkpoint(cmd.chain_name() + ".param", *chain_driver, *model);
    StandardTracer trace(*model, cmd.chain_name());
    chain_driver->add(*model);
    chain_driver->add(console_logger);
    chain_driver->add(chain_checkpoint);
    chain_driver->add(trace);
    chain_driver->go();

    delete chain_driver;
    delete model;
}
