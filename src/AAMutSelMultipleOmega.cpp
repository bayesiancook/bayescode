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

    ValueArg<int> ncat{"", "ncat",
        "Number of components for the amino-acid fitness profiles "
        "(truncation of the stick-breaking process).",
        false, 30, "int", cmd};
    ValueArg<std::string> profiles{"", "profiles",
        "File path the fitness profiles (tsv or csv), thus considered fixed. "
        "Each line must contains the fitness of each of the 20 amino-acid, thus summing to one. "
        "If same number of profiles as the codon alignment, site allocations are considered fixed. "
        "If smaller than the alignment size, site allocations are computed and `ncat` is given by "
        "the number of profiles in the file.",
        false, "Null", "string", cmd};
    SwitchArg flatfitness{"", "flatfitness",
        "Fitness profiles are flattened (and `ncat` equals to 1). "
        "This option is not compatible with the option `profiles`.",
        cmd, false};
    SwitchArg freeomega{"", "freeomega",
        "ω is allowed to vary (default ω is 1.0). "
        "Combined with the option `flatfitness`, we obtain the classical, ω-based codon model (Muse & Gaut). "
        "Without the option `flatfitness`, we obtain the mutation-selection codon model with a multiplicative factor (ω⁎).",
        cmd, false};
    ValueArg<int> omegancat{
        "", "omegancat", "Number of components for ω (finite mixture).", false, 1, "int", cmd};
    ValueArg<double> omegashift{"", "omegashift",
        "Additive shift applied to all ω (0.0 for the general case).", false, 0.0, "double", cmd};
    ValueArg<std::string> omegaarray{"", "omegaarray",
        "File path to ω values (one ω per line), thus considered fixed. "
        "`freeomega` is overridden to false and `omegancat` equals to the number of ω in the file.",
        false, "Null", "string", cmd};

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
    ChainCmdLine cmd{argc, argv, "AAMutSelMultipleOmega", ' ', "1.1.2"};

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
            aamutsel_args.profiles.getValue(), aamutsel_args.omegamode(),
            aamutsel_args.ncat.getValue(), 1, aamutsel_args.omegancat.getValue(),
            aamutsel_args.omegashift.getValue(), aamutsel_args.flatfitness.getValue(),
            aamutsel_args.omegaarray.getValue());
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
