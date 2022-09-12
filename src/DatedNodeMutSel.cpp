#include <cmath>
#include <fstream>
#include "DatedNodeMutSelModel.hpp"
#include "components/ChainCheckpoint.hpp"
#include "components/ChainDriver.hpp"
#include "components/ConsoleLogger.hpp"
#include "components/InferenceAppArgParse.hpp"
#include "components/StandardTracer.hpp"
#include "components/restart_check.hpp"

using namespace std;

class DatedNodeMutselArgParse : public BaseArgParse {
  public:
    explicit DatedNodeMutselArgParse(ChainCmdLine &cmd) : BaseArgParse(cmd) {}

    ValueArg<int> ncat{"", "ncat",
        "Number of components for the amino-acid fitness profiles "
        "(truncation of the stick-breaking process).",
        false, 30, "int", cmd};
    /*
    SwitchArg condition_aware{"b", "condition_aware",
        "One Ne per condition, if the tree doesn't have condition, then one Ne per branch
    (experimental).", cmd, false};
    */
    ValueArg<std::string> traitsfile{
        "", "traitsfile", "File path to the life-history trait (in log-space) in tsv format. "
        "The First column is `TaxonName` (taxon matching the name in the alignment) and the next columns are traits."
        "", false, "Null", "string", cmd};
    ValueArg<std::string> profiles{"", "profiles",
        "File path the fitness profiles (tsv or csv), thus considered fixed. "
        "Each line must contains the fitness of each of the 20 amino-acid, thus summing to one. "
        "If same number of profiles as the codon alignment, site allocations are considered fixed. "
        "If smaller than the alignment size, site allocations are computed and `ncat` is given by "
        "the number of profiles in the file.",
        false, "Null", "string", cmd};
    SwitchArg move_root_pop_size{"", "move_root_pop_size",
        "Move Ne at the root, for the equilibrium frequencies (experimental)", cmd, false};
    SwitchArg clamp_pop_sizes{
        "", "clamp_pop_sizes", "Clamp the branch population size (experimental).", cmd, false};
    SwitchArg clamp_nuc_matrix{
        "", "clamp_nuc_matrix", "Clamp the nucleotide matrix (experimental).", cmd, false};
    SwitchArg clamp_corr_matrix{
        "", "clamp_corr_matrix", "Clamp the correlation matrix (experimental).", cmd, false};
    SwitchArg polymorphism_aware{
        "p", "polymorphism_aware", "Use polymorphic data (experimental).", cmd, false};
    ValueArg<unsigned> precision{"", "precision",
        "The precision of Poisson-Random-Field computations (experimental).", false, 6, "unsigned",
        cmd};
    SwitchArg arithmetic{
        "d", "arithmetic", "Use arithmetic mean instead of geometric (experimental).", cmd, false};
    ValueArg<std::string> fossils{"", "fossils",
        "File path to the fossils calibration in tsv format "
        "with columns `NodeName`, `Age, `LowerBound` and `UpperBound`.",
        false, "Null", "string", cmd};
    ValueArg<int> prior_cov_df{
        "", "df", "Invert Wishart degree of freedom (experimental).", false, 0, "int", cmd};
    SwitchArg uniq_kappa{"", "uniq_kappa",
        "Unique kappa for the invert Wishart matrix prior, "
        "otherwise 1 for each dimension (experimental).",
        cmd, false};

    void check() {
        if (profiles.getValue() != "Null") {
            cout << "Preferences are clamped (option [--profiles <string>]), thus options [--ncat] "
                    "and [-basencat] are not used"
                 << endl;
        }
    }
};

int main(int argc, char *argv[]) {
    ChainCmdLine cmd{argc, argv, "DatedMutSel", ' ', "1.1.2"};

    ChainDriver *chain_driver = nullptr;
    unique_ptr<DatedNodeMutSelModel> model = nullptr;

    if (cmd.resume_from_checkpoint()) {
        std::ifstream is = cmd.checkpoint_file();
        chain_driver = new ChainDriver(is);
        is >> model;
        check_restart(*model, cmd.chain_name() + ".trace");
    } else {
        InferenceAppArgParse inference_args(cmd);
        DatedNodeMutselArgParse args(cmd);
        cmd.parse();
        args.check();
        chain_driver = new ChainDriver(
            cmd.chain_name(), inference_args.every.getValue(), inference_args.until.getValue());
        model = std::make_unique<DatedNodeMutSelModel>(inference_args.alignment.getValue(),
            inference_args.treefile.getValue(), args.traitsfile.getValue(),
            args.profiles.getValue(), args.ncat.getValue(), 1.0, false,
            args.polymorphism_aware.getValue(), args.precision.getValue(),
            args.arithmetic.getValue(), args.move_root_pop_size.getValue(),
            args.clamp_pop_sizes.getValue(), args.clamp_nuc_matrix.getValue(),
            args.clamp_corr_matrix.getValue(), args.fossils.getValue(),
            args.prior_cov_df.getValue(), args.uniq_kappa.getValue());
        model->Update();
    }
    model->ResampleSub(1.0);
    model->MoveParameters(10);
    ConsoleLogger console_logger;
    ChainCheckpoint chain_checkpoint(cmd.chain_name() + ".param", *chain_driver, *model);
    StandardTracer trace(*model, cmd.chain_name());
    chain_driver->add(*model);
    chain_driver->add(console_logger);
    chain_driver->add(chain_checkpoint);
    chain_driver->add(trace);
    chain_driver->go();
}
