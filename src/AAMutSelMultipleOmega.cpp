#include <cmath>
#include <fstream>
#include "AAMutSelMultipleOmegaModel.hpp"
#include "components/ChainCheckpoint.hpp"
#include "components/ChainDriver.hpp"
#include "components/ConsoleLogger.hpp"
#include "components/InferenceAppArgParse.hpp"
#include "components/StandardTracer.hpp"
#include "components/restart_check.hpp"

using namespace std;

class AAMutselArgParse : public BaseArgParse {
  public:
    AAMutselArgParse(ChainCmdLine &cmd) : BaseArgParse(cmd) {}

    SwitchArg flatfitness{
        "", "flatfitness", "Fitness landscape are flattened (and 'ncat' equals to 1).", cmd, false};
    ValueArg<int> ncat{"", "ncat",
        "Truncation of the first-level stick-breaking process (number of fitness profiles).", false,
        30, "int", cmd};
    ValueArg<int> basencat{"", "basencat", "Truncation of the second-level stick-breaking process.",
        false, 1, "int", cmd};
    ValueArg<std::string> fitness_profiles{"", "fitness_profiles",
        "Preferences profiles passed as input. If same size (number of sites) of the codon "
        "alignment, site allocations are considered fixed. If smaller than the alignment size, "
        "site allocations are computed and 'ncat' is given by the number of profiles in the file.",
        false, "Null", "string", cmd};

    ValueArg<double> omegashift{"", "omegashift",
        "Additive shift applied to omega (typically 1 for detecting adaptation, 0 for general case).", false,
        0.0, "double", cmd};
    ValueArg<std::string> deltaomegaarray{"", "deltaomegaarray",
        "Delta-omega array file passed as input. The values are considered fixed, hence "
        "'freeomega' is overridden to false and 'omegancat' equals to the number of values in the array",
        false, "Null", "string", cmd};
    ValueArg<int> omegancat{
        "", "omegancat", "Number of components of omega finite mixture.", false, 1, "int", cmd};
    SwitchArg freeomega{"", "freeomega",
        "Omega is allowed to vary. This parameter is used only if 'deltaomegaarray' is not passed as input.",
        cmd, false};


    //! - omegamode: omega fixed (3), shared across genes (2) or estimated with
    //! shrinkage across genes (1) or without shrinkage (0)
    int omegamode() {
        if (freeomega.getValue()) {
            return 0;
        } else {
            return 3;
        }
    }
};

int main(int argc, char *argv[]) {
    ChainCmdLine cmd{argc, argv, "AAMutSelMultipleOmega", ' ', "0.1"};

    ChainDriver *chain_driver = nullptr;
    AAMutSelMultipleOmegaModel *model = nullptr;

    if (cmd.resume_from_checkpoint()) {
        std::ifstream is = cmd.checkpoint_file();
        chain_driver = new ChainDriver(is);
        model = new AAMutSelMultipleOmegaModel(is);
        check_restart(*model, cmd.chain_name() + ".trace");
    } else {
        InferenceAppArgParse args(cmd);
        AAMutselArgParse aamutsel_args(cmd);
        cmd.parse();
        chain_driver =
            new ChainDriver(cmd.chain_name(), args.every.getValue(), args.until.getValue());
        model = new AAMutSelMultipleOmegaModel(args.alignment.getValue(), args.treefile.getValue(),
            aamutsel_args.fitness_profiles.getValue(), aamutsel_args.omegamode(),
            aamutsel_args.ncat.getValue(), aamutsel_args.basencat.getValue(),
            aamutsel_args.omegancat.getValue(), aamutsel_args.omegashift.getValue(),
            aamutsel_args.flatfitness.getValue(), aamutsel_args.deltaomegaarray.getValue());
    }

    ConsoleLogger console_logger;
    ChainCheckpoint chain_checkpoint(cmd.chain_name() + ".param", *chain_driver, *model);
    StandardTracer trace(*model, cmd.chain_name());
    chain_driver->add(*model);
    chain_driver->add(console_logger);
    chain_driver->add(chain_checkpoint);
    chain_driver->add(trace);
    chain_driver->go();
}
