#include <cmath>
#include <fstream>
#include "BaseArgParse.hpp"
#include "Chain.hpp"
#include "ChainCheckpoint.hpp"
#include "ChainDriver.hpp"
#include "ConsoleLogger.hpp"
#include "SingleOmegaModel.hpp"

using namespace std;

class SingleOmegaArgParse : public BaseArgParse {
  public:
    SingleOmegaArgParse(ChainCmdLine &cmd) : BaseArgParse(cmd) {}
    ValueArg<string> alignment{"a",      "alignment", "Alignment file (PHYLIP)", true, "",
                               "string", cmd};
    ValueArg<string> treefile{"t", "tree", "Tree file (NHX)", true, "", "string", cmd};
    ValueArg<int> every{"e",   "every", "Number of iterations between two traces", false, 1,
                        "int", cmd};
    ValueArg<int> until{"u",   "until", "Maximum number of (saved) iterations (-1 means unlimited)",
                        false, -1,      "int",
                        cmd};
    SwitchArg force{"f", "force", "Overwrite existing output files", cmd};
};

class SingleOmegaTrace : public ChainComponent {
    Tracer model_tracer;
    Tracer stats_tracer;
    std::string chain_name;

public:
    SingleOmegaTrace(SingleOmegaModel& m, std::string chain_name) : model_tracer(m,&SingleOmegaModel::declare_model),
                                                                    stats_tracer(m,&SingleOmegaModel::declare_stats),
                                                                    chain_name(chain_name) {}
    void start() {
        std::ofstream model_os{chain_name + ".chain", std::ios_base::trunc};
        model_tracer.write_header(model_os);
        std::ofstream stats_os{chain_name + ".trace", std::ios_base::trunc};
        stats_tracer.write_header(stats_os);
    }
    void savepoint(int) {
        std::ofstream model_os{chain_name + ".chain", std::ios_base::app};
        model_tracer.write_line(model_os);
        std::ofstream trace_os{chain_name + ".trace", std::ios_base::app};
        stats_tracer.write_line(trace_os);
    }
};

int main(int argc, char *argv[]) {
    ChainCmdLine cmd{argc, argv, "SingleOmega", ' ', "0.1"};

    ChainDriver *chain_driver = nullptr;
    SingleOmegaModel *model = nullptr;

    if (cmd.resume_from_checkpoint()) {
        std::ifstream is = cmd.checkpoint_file();
        chain_driver = new ChainDriver(is);
        model = new SingleOmegaModel(is);
    } else {
        SingleOmegaArgParse args(cmd);
        cmd.parse();
        chain_driver =
            new ChainDriver(cmd.chain_name(), args.every.getValue(), args.until.getValue());
        model = new SingleOmegaModel(args.alignment.getValue(), args.treefile.getValue());
    }

    ConsoleLogger console_logger;
    ChainCheckpoint chain_checkpoint(cmd.chain_name() + ".param", *chain_driver, *model);
    SingleOmegaTrace trace(*model, cmd.chain_name());
    chain_driver->add(*model);
    chain_driver->add(console_logger);
    chain_driver->add(chain_checkpoint);
    chain_driver->add(trace);
    chain_driver->go();

    delete chain_driver;
    delete model;
}
