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
    //     ValueArg<int> burnin{"b", "burnin", "Burnin iterations that will be discarded", false, 0,
    //     "int", cmd}; No more burnin: if the model has difficulty starting, we'll think about it
    //     again.
};


int main(int argc, char *argv[]) {
    string message =
        "\t\tDiffSelDoublySparse\n\tA doubly-sparse version of the differential selection "
        "model\nthe model defines K conditions: background (k=0) and alternative "
        "(k=1..K-1);\nbranches are a priori allocated to any one of the K conditions\nThe model is "
        "based on the following system of random variables:\n - G_kia > 0     : an array of site- "
        "and condition-specific pre-fitness parameters, for condition k=0..K-1, site i and "
        "amino-acid a;\n - m_ia =  0 or 1: an array of site-specific masks\n - d_kia = 0 or 1: an "
        "array of site and condition-specific toggles (for conditions k=1..K-1).\nthe "
        "site-specific masks determine whether an amino-acid is or is not allowed at a given "
        "site;\nif allowed, then the toggles determine whether the fitness for the amino-acid at "
        "that site changes upon going from condition 0 to condition k;\nquantitatively, the "
        "fitness vector at site i under condition k F_kia is defined as:\n  F_0ia = m_ia * G_0ia + "
        "(1 - m_ia) * epsilon\nF_kia = m_ia * d_kia * F_kia + m_ia * (1 - m_ia * d_kia) * F_0ia, "
        "for k=1..K-1\nwhere epsilon << 1 is a residual fitness for background amino-acids\npriors "
        "and distributions:\nG_kia ~ Exponential(1.0)\n(m_ia)_{a=1..20} ~ IID_Bernoulli(maskprob), "
        "conditional on at least one m_ia being equal to 1\nd_kia ~ "
        "Bernoulli(shiftprob_k)\n\nmaskprob and epsilon are estimated by default (uniform prior in "
        "both cases)\nshiftprob_k is estimated and is under a mixed prior:\n - with probability "
        "1-pi_k, shiftprob_k = 0\n - with probability pi_k, shiftprob_k ~ "
        "Beta(shiftprobhypermean_k, shiftprobhyperinvconc_k)\nfor each condition, the model has 3 "
        "thus three key hyperparameters: pi_k, shiftprobhypermean_k and "
        "shiftprobhyperinvconc_k;\nin single gene analyses, these parameters are always the same "
        "across conditions and can be specified by the user;\nin multi-gene analyses, these "
        "parameters can be estimated across genes (separately for each condition, see "
        "multigenediffseldsparse).\nStatistical support for a differential selection effect in "
        "condition k is quantified:\n - for a given site i and a given amino acid a, by the "
        "posterior probability that d_kia is equal to 1;\n - for a given site i, by the posterior "
        "probability that d_kia == 1 for at least one amino-acid a;\n - for the gene: by the "
        "posterior probability that shiftprob_k > 0;\nprogram options:\n\t-f: force overwrite of "
        "already existing chain\n\t-x <every> <until>: saving frequency and stopping time "
        "(default: every = 1, until = -1)\n model options:\n\t-ncond <ncond>:  specify number of "
        "conditions\n\t-pi <pi>: specify the value for pi_k (same for all k=1..K-1, default value "
        "is 0.1)\n\t-shiftprob <mean> <invconc>: specify the values for shiftprobhypermean_k and "
        "shiftprobhyperinvconc_k (same for all k=1..K-1, default: 0.1 amd 0.1\n";


    ChainCmdLine cmd{argc, argv, message, ' ', "0.1"};

    ChainDriver *chain_driver = nullptr;
    unique_ptr<DiffSelDoublySparseModel> model = nullptr;

    // Default values, as in the original version:
    int codonmodel = 1;

    if (cmd.resume_from_checkpoint()) {
        ifstream is = cmd.checkpoint_file();
        chain_driver = new ChainDriver(is);
        //  model = new DiffSelDoublySparseModel(is); TODO : IMPLEMENT
    } else {
        InferenceAppArgParse args(cmd);
        DiffSelDoublySparseAppArgParse ddargs(cmd);
        cmd.parse();
        chain_driver =
            new ChainDriver(cmd.chain_name(), args.every.getValue(), args.until.getValue());
        model = unique_ptr<DiffSelDoublySparseModel>(new DiffSelDoublySparseModel(
            args.alignment.getValue(), args.treefile.getValue(), ddargs.ncond.getValue(),
            ddargs.nlevel.getValue(), codonmodel, ddargs.epsilon.getValue(),
            ddargs.fitnessshape.getValue(), ddargs.pihypermean.getValue(),
            ddargs.shiftprobmean.getValue(), ddargs.shiftprobinvconc.getValue(),
            param_mode_t(ddargs.fitnesscentermode.getValue()), true));
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