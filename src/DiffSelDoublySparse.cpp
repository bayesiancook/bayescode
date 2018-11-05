#include <cmath>
#include <fstream>
#include "DiffSelDoublySparseModel.hpp"
#include "InferenceAppArgParse.hpp"
#include "components/ChainCheckpoint.hpp"
#include "components/ChainDriver.hpp"
#include "components/ConsoleLogger.hpp"
#include "components/StandardTracer.hpp"

using namespace std;


 class DiffSelDoublySparseAppArgParse : public BaseArgParse {
   public:
     DiffSelDoublySparseAppArgParse(ChainCmdLine &cmd) : BaseArgParse(cmd) {}
     ValueArg<int> ncond{"", "ncond", "Number of conditions", false, 1, "int", cmd};
     ValueArg<int> nlevel{"", "nlevel", "Number of levels", false, 1, "int", cmd};
     ValueArg<double> fitnessshape{"", "shape", "Shape of the fitness distribution", false, 1, "double", cmd}; // We don't do the "free" option that was available in the previous version
     ValueArg<double> fitnesscentermode{"", "center", "Center mode of the fitness distribution", false, 3, "double", cmd};
     ValueArg<double> epsilon{"", "epsilon", "Epsilon value for rare amino acids", false, -1, "double", cmd};
     ValueArg<double> pihypermean{"", "pihypermean", "Probability that there is a differential effect in a condition", false, 0.1, "double", cmd};
     ValueArg<double> shiftprobmean{"", "shiftprobmean", "Mean of the Beta probability of profile change", false, 0.1, "double", cmd};
     ValueArg<double> shiftprobinvconc{"", "shiftprobinvconc", "Inverse concentration of the Beta probability of profile change (0 means very pointed probability)", false, 0.1, "double", cmd};
//     ValueArg<int> burnin{"b", "burnin", "Burnin iterations that will be discarded", false, 0, "int", cmd}; No more burnin: if the model has difficulty starting, we'll think about it again.

 };


int main(int argc, char *argv[]) {


    string message =
                      "\t\tDiffSelDoublySparse\n\tA doubly-sparse version of the differential selection model\nthe model defines K conditions: background (k=0) and alternative (k=1..K-1);\nbranches are a priori allocated to any one of the K conditions\nThe model is based on the following system of random variables:\n - G_kia > 0     : an array of site- and condition-specific pre-fitness parameters, for condition k=0..K-1, site i and amino-acid a;\n - m_ia =  0 or 1: an array of site-specific masks\n - d_kia = 0 or 1: an array of site and condition-specific toggles (for conditions k=1..K-1).\nthe site-specific masks determine whether an amino-acid is or is not allowed at a given site;\nif allowed, then the toggles determine whether the fitness for the amino-acid at that site changes upon going from condition 0 to condition k;\nquantitatively, the fitness vector at site i under condition k F_kia is defined as:\n  F_0ia = m_ia * G_0ia + (1 - m_ia) * epsilon\nF_kia = m_ia * d_kia * F_kia + m_ia * (1 - m_ia * d_kia) * F_0ia, for k=1..K-1\nwhere epsilon << 1 is a residual fitness for background amino-acids\npriors and distributions:\nG_kia ~ Exponential(1.0)\n(m_ia)_{a=1..20} ~ IID_Bernoulli(maskprob), conditional on at least one m_ia being equal to 1\nd_kia ~ Bernoulli(shiftprob_k)\n\nmaskprob and epsilon are estimated by default (uniform prior in both cases)\nshiftprob_k is estimated and is under a mixed prior:\n - with probability 1-pi_k, shiftprob_k = 0\n - with probability pi_k, shiftprob_k ~ Beta(shiftprobhypermean_k, shiftprobhyperinvconc_k)\nfor each condition, the model has 3 thus three key hyperparameters: pi_k, shiftprobhypermean_k and shiftprobhyperinvconc_k;\nin single gene analyses, these parameters are always the same across conditions and can be specified by the user;\nin multi-gene analyses, these parameters can be estimated across genes (separately for each condition, see multigenediffseldsparse).\nStatistical support for a differential selection effect in condition k is quantified:\n - for a given site i and a given amino acid a, by the posterior probability that d_kia is equal to 1;\n - for a given site i, by the posterior probability that d_kia == 1 for at least one amino-acid a;\n - for the gene: by the posterior probability that shiftprob_k > 0;\nprogram options:\n\t-f: force overwrite of already existing chain\n\t-x <every> <until>: saving frequency and stopping time (default: every = 1, until = -1)\n model options:\n\t-ncond <ncond>:  specify number of conditions\n\t-pi <pi>: specify the value for pi_k (same for all k=1..K-1, default value is 0.1)\n\t-shiftprob <mean> <invconc>: specify the values for shiftprobhypermean_k and shiftprobhyperinvconc_k (same for all k=1..K-1, default: 0.1 amd 0.1\n";


    ChainCmdLine cmd{argc, argv, message, ' ', "0.1"};

    ChainDriver *chain_driver = nullptr;
    DiffSelDoublySparseModel *model = nullptr;

    // Default values, as in the original version:
    int codonmodel = 1;

    if (cmd.resume_from_checkpoint()) {
        std::ifstream is = cmd.checkpoint_file();
        chain_driver = new ChainDriver(is);
      //  model = new DiffSelDoublySparseModel(is); TODO : IMPLEMENT
    } else {
        InferenceAppArgParse args(cmd);
        DiffSelDoublySparseAppArgParse ddargs (cmd) ;
        cmd.parse();
        chain_driver =
            new ChainDriver(cmd.chain_name(), args.every.getValue(), args.until.getValue());
        model = new DiffSelDoublySparseModel(args.alignment.getValue(), args.treefile.getValue(),
        ddargs.ncond.getValue(), ddargs.nlevel.getValue(), codonmodel, ddargs.epsilon.getValue(),
        ddargs.fitnessshape.getValue(), ddargs.pihypermean.getValue(), ddargs.shiftprobmean.getValue(), ddargs.shiftprobinvconc.getValue() );


        model->SetWithToggles(true);

        model->SetFitnessCenterMode(ddargs.fitnesscentermode.getValue());
        model->Allocate();
        model->Update();


    }

    ConsoleLogger console_logger;
    ChainCheckpoint chain_checkpoint(cmd.chain_name() + ".param", *chain_driver, *model);
    //StandardTracer trace(*model, cmd.chain_name());
    chain_driver->add(*model);
    chain_driver->add(console_logger);
    chain_driver->add(chain_checkpoint);
    //chain_driver->add(trace);
    chain_driver->go();

    delete chain_driver;
    delete model;
}




















#include <cmath>
#include <fstream>
#include "Chain.hpp"
#include "DiffSelDoublySparseModel.hpp"
using namespace std;

/**
 * \brief Chain object for running an MCMC under the DiffSelDoublySparseModel
 */

class DiffSelDoublySparseChain : public Chain {
  private:
    // Chain parameters
    string modeltype, datafile, treefile;
    int ncond, nlevel, codonmodel;
    double fitnessshape;
    int fitnesscentermode;
    // -1: estimated
    // 1 : no mask
    double epsilon;
    double pihypermean;
    double shiftprobmean;
    double shiftprobinvconc;
    int burnin;

  public:
    DiffSelDoublySparseModel *GetModel() { return static_cast<DiffSelDoublySparseModel *>(model); }

    string GetModelType() override { return modeltype; }

    DiffSelDoublySparseChain(string indata, string intree, int inncond, int innlevel,
        int incodonmodel, double infitnessshape, int infitnesscentermode, double inepsilon,
        double inpihypermean, double inshiftprobmean, double inshiftprobinvconc, int inburnin,
        int inevery, int inuntil, int insaveall, string inname, int force)
        : modeltype("DIFFSELSPARSE"),
          datafile(indata),
          treefile(intree),
          ncond(inncond),
          nlevel(innlevel),
          codonmodel(incodonmodel),
          fitnessshape(infitnessshape),
          fitnesscentermode(infitnesscentermode),
          epsilon(inepsilon),
          pihypermean(inpihypermean),
          shiftprobmean(inshiftprobmean),
          shiftprobinvconc(inshiftprobinvconc) {
        burnin = inburnin;
        every = inevery;
        until = inuntil;
        saveall = insaveall;
        name = inname;
        New(force);
    }

    DiffSelDoublySparseChain(string filename) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new DiffSelDoublySparseModel(datafile, treefile, ncond, nlevel, codonmodel, epsilon,
            fitnessshape, pihypermean, shiftprobmean, shiftprobinvconc);
        if (burnin) {
            GetModel()->SetWithToggles(0);
        } else {
            GetModel()->SetWithToggles(1);
        }
        GetModel()->SetFitnessCenterMode(fitnesscentermode);
        GetModel()->Allocate();
        GetModel()->Update();
        Reset(force);
        cerr << "-- initial ln prob = " << GetModel()->GetLogProb() << "\n";
        model->Trace(cerr);
    }

    void Open() override {
        ifstream is((name + ".param").c_str());
        if (!is) {
            cerr << "-- Error : cannot find file : " << name << ".param\n";
            exit(1);
        }
        is >> modeltype;
        is >> datafile >> treefile >> ncond >> nlevel;
        is >> codonmodel;
        is >> fitnessshape;
        is >> fitnesscentermode;
        is >> epsilon;
        is >> pihypermean >> shiftprobmean >> shiftprobinvconc;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> burnin;
        is >> every >> until >> saveall >> size;

        if (modeltype == "DIFFSELSPARSE") {
            model = new DiffSelDoublySparseModel(datafile, treefile, ncond, nlevel, codonmodel,
                epsilon, fitnessshape, pihypermean, shiftprobmean, shiftprobinvconc);
        } else {
            cerr << "-- Error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }
        if (size < burnin) {
            GetModel()->SetWithToggles(0);
        } else {
            GetModel()->SetWithToggles(1);
        }
        GetModel()->SetFitnessCenterMode(fitnesscentermode);
        GetModel()->Allocate();
        GetModel()->FromStream(is);
        GetModel()->Update();
        cerr << size << " points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
        model->Trace(cerr);
    }

    void Save() override {
        if (size == burnin) { GetModel()->SetWithToggles(1); }

        ofstream param_os((name + ".param").c_str());
        param_os << GetModelType() << '\n';
        param_os << datafile << '\t' << treefile << '\t' << ncond << '\t' << nlevel << '\n';
        param_os << codonmodel << '\n';
        param_os << fitnessshape << '\n';
        param_os << fitnesscentermode << '\n';
        param_os << epsilon << '\n';
        param_os << pihypermean << '\t' << shiftprobmean << '\t' << shiftprobinvconc << '\n';
        param_os << 0 << '\n';
        param_os << burnin << '\n';
        param_os << every << '\t' << until << '\t' << saveall << '\t' << size << '\n';

        model->ToStream(param_os);
    }

    void SavePoint() override {
        Chain::SavePoint();
        for (int k = 0; k < ncond; k++) {
            ostringstream s;
            s << name << "_" << k;
            if (k) {
                ofstream tos((s.str() + ".shifttoggle").c_str(), ios_base::app);
                GetModel()->TraceToggle(k, tos);
            }
            ofstream fos((s.str() + ".fitness").c_str(), ios_base::app);
            GetModel()->TraceFitness(k, fos);
        }
    }

    void MakeFiles(int force) override {
        Chain::MakeFiles(force);
        for (int k = 0; k < ncond; k++) {
            ostringstream s;
            s << name << "_" << k;
            if (k) { ofstream tos((s.str() + ".shifttoggle").c_str()); }
            ofstream fos((s.str() + ".fitness").c_str());
        }
    }
};



/*
int main(int argc, char *argv[]) {
    string name = "";
    DiffSelDoublySparseChain *chain = 0;

    if (argc == 1) {
        cerr << '\n';
        cerr << "A doubly-sparse version of the differential selection model\n";
        cerr << "the model defines K conditions: background (k=0) and alternative (k=1..K-1);\n";
        cerr << "branches are a priori allocated to any one of the K conditions\n";
        cerr << "The model is based on the following system of random variables:\n";
        cerr << " - G_kia > 0     : an array of site- and condition-specific pre-fitness "
                "parameters, for condition k=0..K-1, site i and amino-acid a;\n";
        cerr << " - m_ia =  0 or 1: an array of site-specific masks\n";
        cerr << " - d_kia = 0 or 1: an array of site and condition-specific toggles (for "
                "conditions k=1..K-1).\n";
        cerr << "the site-specific masks determine whether an amino-acid is or is not allowed at a "
                "given site;\n";
        cerr << "if allowed, then the toggles determine whether the fitness for the amino-acid at "
                "that site changes upon going from condition 0 to condition k;\n";
        cerr << "quantitatively, the fitness vector at site i under condition k F_kia is defined "
                "as:\n";
        cerr << "  F_0ia = m_ia * G_0ia + (1 - m_ia) * epsilon\n";
        cerr
            << "  F_kia = m_ia * d_kia * F_kia + m_ia * (1 - m_ia * d_kia) * F_0ia, for k=1..K-1\n";
        cerr << "where epsilon << 1 is a residual fitness for background amino-acids\n";
        cerr << '\n';
        cerr << "priors and distributions:\n";
        cerr << "G_kia ~ Exponential(1.0)\n";
        cerr << "(m_ia)_{a=1..20} ~ IID_Bernoulli(maskprob), conditional on at least one m_ia "
                "being equal to 1\n";
        cerr << "d_kia ~ Bernoulli(shiftprob_k)\n";
        cerr << "\n";
        cerr << "maskprob and epsilon are estimated by default (uniform prior in both cases)\n";
        cerr << "shiftprob_k is estimated and is under a mixed prior:\n";
        cerr << " - with probability 1-pi_k, shiftprob_k = 0\n";
        cerr << " - with probability pi_k, shiftprob_k ~ Beta(shiftprobhypermean_k, "
                "shiftprobhyperinvconc_k)\n";
        cerr << "for each condition, the model has 3 thus three key hyperparameters: pi_k, "
                "shiftprobhypermean_k and shiftprobhyperinvconc_k;\n";
        cerr << "in single gene analyses, these parameters are always the same across conditions "
                "and can be specified by the user;\n";
        cerr << "in multi-gene analyses, these parameters can be estimated across genes "
                "(separately for each condition, see multigenediffseldsparse).\n";
        cerr << '\n';
        cerr << "Statistical support for a differential selection effect in condition k is "
                "quantified:\n";
        cerr << " - for a given site i and a given amino acid a, by the posterior probability that "
                "d_kia is equal to 1;\n";
        cerr << " - for a given site i, by the posterior probability that d_kia == 1 for at least "
                "one amino-acid a;\n";
        cerr << " - for the gene: by the posterior probability that shiftprob_k > 0;\n";
        // cerr << "of note, setting pi_k = 0.5, and noting pp_k the post prob that shiftprob_k > 0,
        // then the Bayes factor in favor of the alternative against the null is given by BF = pp_k
        // / (1 - pp_k)\n";
        cerr << "\n";
        cerr << "command: diffseldsparse -d <alignment> -t <treefile> <chainname>\n";
        cerr << '\n';
        cerr << "program options:\n";
        cerr << "\t-f: force overwrite of already existing chain\n";
        cerr << "\t-x <every> <until>: saving frequency and stopping time (default: every = 1, "
                "until = -1)\n";
        cerr << '\n';
        cerr << "model options:\n";
        cerr << "\t-ncond <ncond>:  specify number of conditions\n";
        cerr << "\t-pi <pi>: specify the value for pi_k (same for all k=1..K-1, default value is "
                "0.1)\n";
        cerr << "\t-shiftprob <mean> <invconc>: specify the values for shiftprobhypermean_k and "
                "shiftprobhyperinvconc_k (same for all k=1..K-1, default: 0.1 amd 0.1\n";
        cerr << '\n';

        exit(0);
    }

    // this is an already existing chain on the disk; reopen and restart
    if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        chain = new DiffSelDoublySparseChain(name);
    }

    // this is a new chain
    else {
        string datafile = "";
        string treefile = "";
        int ncond = 2;
        int nlevel = 2;
        int codonmodel = 1;
        double epsilon = -1;

        double fitnessshape = 20;
        int fitnesscentermode = 3;

        double pihypermean = 0.1;
        double shiftprobmean = 0.1;
        double shiftprobinvconc = 0.1;

        name = "";
        int burnin = 0;
        int every = 1;
        int until = -1;
        int saveall = 1;
        int force = 0;

        try {
            if (argc == 1) { throw(0); }

            int i = 1;
            while (i < argc) {
                string s = argv[i];

                if (s == "-d") {
                    i++;
                    datafile = argv[i];
                } else if ((s == "-t") || (s == "-T")) {
                    i++;
                    treefile = argv[i];
                } else if (s == "-f") {
                    force = 1;
                } else if (s == "+s") {
                    saveall = 1;
                } else if (s == "-s") {
                    saveall = 0;
                } else if (s == "-ncond") {
                    i++;
                    ncond = atoi(argv[i]);
                } else if (s == "-nlevel") {
                    i++;
                    nlevel = atoi(argv[i]);
                } else if (s == "-shape") {
                    i++;
                    string tmp = argv[i];
                    if (tmp == "free") {
                        fitnessshape = 0;
                    } else {
                        fitnessshape = atof(argv[i]);
                    }
                } else if (s == "-center") {
                    i++;
                    string tmp = argv[i];
                    if (tmp == "free") {
                        fitnesscentermode = 0;
                    } else if ((tmp == "fixed") || (tmp == "uniform")) {
                        fitnesscentermode = 3;
                    }
                } else if ((s == "-eps") || (s == "-epsilon")) {
                    i++;
                    string tmp = argv[i];
                    if (tmp == "free") {
                        epsilon = -1;
                    } else {
                        epsilon = atof(argv[i]);
                    }
                } else if (s == "-pi") {
                    i++;
                    pihypermean = atof(argv[i]);
                } else if (s == "-shiftprob") {
                    i++;
                    shiftprobmean = atof(argv[i]);
                    i++;
                    shiftprobinvconc = atof(argv[i]);
                } else if (s == "-b") {
                    i++;
                    burnin = atoi(argv[i]);
                } else if ((s == "-x") || (s == "-extract")) {
                    i++;
                    if (i == argc) throw(0);
                    every = atoi(argv[i]);
                    i++;
                    if (i == argc) throw(0);
                    until = atoi(argv[i]);
                } else {
                    if (i != (argc - 1)) { throw(0); }
                    name = argv[i];
                }
                i++;
            }
        } catch (...) {
            cerr << "error in command\n";
            exit(1);
        }
        chain = new DiffSelDoublySparseChain(datafile, treefile, ncond, nlevel, codonmodel,
            fitnessshape, fitnesscentermode, epsilon, pihypermean, shiftprobmean, shiftprobinvconc,
            burnin, every, until, saveall, name, force);
    }
    cerr << "chain " << name << " started\n";
    chain->Start();
    cerr << "chain " << name << " stopped\n";
    cerr << chain->GetSize()
         << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
    chain->GetModel()->Trace(cerr);
}*/
