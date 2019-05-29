#include <cmath>
#include <fstream>
#include "AAMutSelDM5Model.hpp"
#include "components/ChainCheckpoint.hpp"
#include "components/ChainDriver.hpp"
#include "components/ConsoleLogger.hpp"
#include "components/InferenceAppArgParse.hpp"
#include "components/StandardTracer.hpp"
#include "components/restart_check.hpp"

using namespace std;

class AAMutselDM5ArgParse : public BaseArgParse {
  public:
    AAMutselDM5ArgParse(ChainCmdLine &cmd) : BaseArgParse(cmd) {}

    ValueArg<int> ncat{
        "", "ncat", "truncation of the first-level stick-breaking process", false, 100, "int", cmd};
    ValueArg<int> basencat{"", "basencat", "truncation of the second-level stick-breaking process",
        false, 1, "int", cmd};
    ValueArg<double> omegashift{"", "omegashift",
        "the shift applied to omega (typically 1 for detecting adaptation, 0 for general case)",
        false, 1.0, "double", cmd};
    ValueArg<int> omegancat{
        "", "omegancat", "number of components of omega finite mixture", false, 1, "int", cmd};
    SwitchArg freeomega{"", "freeomega",
        "omega is allowed to vary with shrinkage (otherwise set to 1)", cmd, false};
    SwitchArg flatfitness{"", "flatfitness", "Fitness landscape are flattened", cmd, false};

    ValueArg<double> omega_weight_p0{"", "p0",
        "The initial value of the proportion of sites with delta_omega equal 0", false, 0.5,
        "double", cmd};
    SwitchArg omega_weight_fixed_p0{
        "", "fixp0", "The proportion of sites with delta_omega is fixed", cmd, false};
    ValueArg<double> omega_p0{
        "", "omega_p0", "the omega applied to the first category", false, 1.0, "double", cmd};

    ValueArg<double> hypermean_threshold{"", "hypermean_threshold",
        "The lower threshold for the mean of the Gamma distribution", false, 0.1, "double", cmd};
    ValueArg<double> hyperinvshape_threshold{"", "hyperinvshape_threshold",
        "The upper threshold for the inverse shape of the Gamma distribution", false, 1.0, "double",
        cmd};

    //! - omegamode: omega fixed (3), shared across genes (2) or estimated with
    //! shrinkage across genes (1) or without shrinkage (0)
    int omegamode() {
        if (freeomega.getValue()) {
            return 1;
        } else {
            return 3;
        }
    }

    //! Because we already have one category for the DeltaOmega=0
    int omegaNcat() { return omegancat.getValue() + 1; }
};

int main(int argc, char *argv[]) {
    ChainCmdLine cmd{argc, argv, "AAMutSelDM5", ' ', "0.1"};

    ChainDriver *chain_driver = nullptr;
    AAMutSelDM5Model *model = nullptr;

    if (cmd.resume_from_checkpoint()) {
        std::ifstream is = cmd.checkpoint_file();
        chain_driver = new ChainDriver(is);
        model = new AAMutSelDM5Model(is);
        check_restart(*model, cmd.chain_name() + ".trace");
    } else {
        InferenceAppArgParse args(cmd);
        AAMutselDM5ArgParse aamutseldm5_args(cmd);
        cmd.parse();
        chain_driver =
            new ChainDriver(cmd.chain_name(), args.every.getValue(), args.until.getValue());
        model = new AAMutSelDM5Model(args.alignment.getValue(), args.treefile.getValue(),
            aamutseldm5_args.omegamode(), aamutseldm5_args.ncat.getValue(),
            aamutseldm5_args.basencat.getValue(), aamutseldm5_args.omegaNcat(),
            aamutseldm5_args.omegashift.getValue(), aamutseldm5_args.flatfitness.getValue(),
            aamutseldm5_args.omega_weight_p0.getValue(),
            aamutseldm5_args.omega_weight_fixed_p0.getValue(), aamutseldm5_args.omega_p0.getValue(),
            aamutseldm5_args.hyperinvshape_threshold.getValue(),
            aamutseldm5_args.hypermean_threshold.getValue());
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
