#include <cmath>
#include <fstream>
#include "DiffSelDoublySparseModel.hpp"
#include "components/ChainCheckpoint.hpp"
#include "components/ChainDriver.hpp"
#include "components/ConsoleLogger.hpp"
#include "components/InferenceAppArgParse.hpp"
#include "components/StandardTracer.hpp"

using namespace std;


class DiffSelDoublySparseAppArgParse : public BaseArgParse {
  public:
    DiffSelDoublySparseAppArgParse(ChainCmdLine &cmd) : BaseArgParse(cmd) {}
    ValueArg<int> ncond{"", "ncond", "Number of conditions", false, 1, "int", cmd};
    ValueArg<int> nlevel{"", "nlevel", "Number of levels", false, 1, "int", cmd};
    ValueArg<double> fitnessshape{"", "shape", "Shape of the fitness distribution", false, 1,
        "double", cmd};  // We don't do the "free" option that was available in the previous version
    ValueArg<double> fitnesscentermode{
        "", "center", "Center mode of the fitness distribution", false, 3, "double", cmd};
    ValueArg<double> epsilon{
        "", "epsilon", "Epsilon value for rare amino acids", false, -1, "double", cmd};
    ValueArg<double> pihypermean{"", "pihypermean",
        "Probability that there is a differential effect in a condition", false, 0.1, "double",
        cmd};
    ValueArg<double> shiftprobmean{"", "shiftprobmean",
        "Mean of the Beta probability of profile change", false, 0.1, "double", cmd};
    ValueArg<double> shiftprobinvconc{"", "shiftprobinvconc",
        "Inverse concentration of the Beta probability of profile change (0 means very pointed "
        "probability)",
        false, 0.1, "double", cmd};
    SwitchArg sitewise{"", "sw", "Use model with site-wise convergence toggles", cmd};
    ValueArg<double> fixed_bl{"", "fixed-bl",
        "Use fixed branch lengths taken from tree file and multiplied by provided value", false,
        1.0, "double", cmd};
};


int main(int argc, char *argv[]) {
    string message =
        "A doubly-sparse version of the differential selection model. The model defines K "
        "conditions: background (k=0) and alternative (k=1..K-1); branches are a priori allocated "
        "to any one of the K conditions. The model is based on the following system of random "
        "variables:\n"
        " - G_kia > 0 : an array of site- and condition-specific pre-fitness parameters, for "
        "condition k=0..K-1, site i and amino-acid a;\n"
        " - m_ia = 0 or 1: an array of site-specific masks\n"
        " - d_kia = 0 or 1: an array of site and condition-specific toggles (for conditions "
        "k=1..K-1).\n"
        "The site-specific masks determine whether an amino-acid is or is not allowed at a given "
        "site; if allowed, then the toggles determine whether the fitness for the amino-acid at "
        "that site changes upon going from condition 0 to condition k; quantitatively, the fitness "
        "vector at site i under condition k F_kia is defined as: F_0ia = m_ia * G_0ia + (1 - m_ia) "
        "* epsilonF_kia = m_ia * d_kia * F_kia + m_ia * (1 - m_ia * d_kia) * F_0ia, for k=1..K-1 "
        "where epsilon << 1 is a residual fitness for background amino-acids priors and "
        "distributions: G_kia ~ Exponential(1.0)(m_ia)_{a=1..20} ~ IID_Bernoulli(maskprob), "
        "conditional on at least one m_ia being equal to 1d_kia ~ Bernoulli(shiftprob_k). maskprob "
        "and epsilon are estimated by default (uniform prior in both cases) shiftprob_k is "
        "estimated and is under a mixed prior:\n"
        " - with probability 1-pi_k, shiftprob_k = 0\n"
        " - with probability pi_k, shiftprob_k ~ Beta(shiftprobhypermean_k, "
        "shiftprobhyperinvconc_k).\n"
        "For each condition, the model has 3  key hyperparameters: pi_k, shiftprobhypermean_k and "
        "shiftprobhyperinvconc_k; in single gene analyses, these parameters are always the same "
        "across conditions and can be specified by the user; in multi-gene analyses, these "
        "parameters can be estimated across genes (separately for each condition, see "
        "multigenediffseldsparse). Statistical support for a differential selection effect in "
        "condition k is quantified:\n"
        " - for a given site i and a given amino acid a, by the posterior probability that d_kia "
        "is equal to 1;\n"
        " - for a given site i, by the posterior probability that d_kia == 1 for at least one "
        "amino-acid a;\n"
        " - for the gene: by the posterior probability that shiftprob_k > 0.";


    ChainCmdLine cmd{argc, argv, message, ' ', "0.1"};

    ChainDriver *chain_driver = nullptr;
    unique_ptr<DiffSelDoublySparseModel> model = nullptr;

    if (cmd.resume_from_checkpoint()) {
        ifstream is = cmd.checkpoint_file();
        chain_driver = new ChainDriver(is);
        FAIL("Resuming from checkpoint not implemented (yet) for diffsel");
        //  model = new DiffSelDoublySparseModel(is); TODO : IMPLEMENT
    } else {
        InferenceAppArgParse args(cmd);
        DiffSelDoublySparseAppArgParse ddargs(cmd);
        cmd.parse();
        chain_driver =
            new ChainDriver(cmd.chain_name(), args.every.getValue(), args.until.getValue());

        // diffsel params
        DiffselDoubleSparseConfig config;
        config.datafile = args.alignment.getValue();
        config.treefile = args.treefile.getValue();
        config.Ncond = ddargs.ncond.getValue();
        config.Nlevel = ddargs.nlevel.getValue();
        config.branch_lengths =
            ddargs.fixed_bl.isSet()
                ? MultiGeneParameter<double, hyper_mean_invshape>(
                      param_mode_t::fixed, ddargs.fixed_bl.getValue())
                : MultiGeneParameter<double, hyper_mean_invshape>(independent, {0.1, 1.0});

        model = make_unique<DiffSelDoublySparseModel>(config);
        model->Update();
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
}