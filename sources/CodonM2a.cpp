#include <cmath>
#include <fstream>
#include "Chain.hpp"
#include "CodonM2aModel.hpp"
using namespace std;

/**
 * \brief Chain object for running an MCMC under CodonM2aModel
 */

class CodonM2aChain : public Chain {
  private:
    // Chain parameters
    string modeltype, datapath, datafile, treefile;
    double pi;
    double puromhypermean, puromhyperinvconc, dposomhypermean, dposomhyperinvshape;
    double purwhypermean, purwhyperinvconc, poswhypermean, poswhyperinvconc;

  public:
    //! return the model, with its derived type (unlike ProbModel::GetModel)
    CodonM2aModel *GetModel() { return static_cast<CodonM2aModel *>(model); }

    string GetModelType() override { return modeltype; }

    //! constructor for a new chain: datapath datafile, treefile, pi (fraction of sites
    //! under positive selection) saving frequency, final chain size, chain name
    //! and overwrite flag -- calls New
    CodonM2aChain(string indatapath, string indatafile, string intreefile, double inpi, 
                    double inpuromhypermean, double inpuromhyperinvconc,
                    double indposomhypermean, double indposomhyperinvshape,
                    double inpurwhypermean, double inpurwhyperinvconc,
                    double inposwhypermean, double inposwhyperinvconc,
                    int inevery, int inuntil,
                    string inname, int force)
        : modeltype("CODONM2A"), datapath(indatapath), datafile(indatafile), treefile(intreefile), pi(inpi) {
        puromhypermean = inpuromhypermean;
        puromhyperinvconc = inpuromhyperinvconc;
        dposomhypermean = indposomhypermean;
        dposomhyperinvshape = indposomhyperinvshape;
        purwhypermean = inpurwhypermean;
        purwhyperinvconc = inpurwhyperinvconc;
        poswhypermean = inposwhypermean;
        poswhyperinvconc = inposwhyperinvconc;
        every = inevery;
        until = inuntil;
        name = inname;
        New(force);
    }

    //! constructor for re-opening an already existing chain from file -- calls
    //! Open
    CodonM2aChain(string filename) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new CodonM2aModel(datapath, datafile, treefile, pi);
        GetModel()->SetMixtureHyperParameters(puromhypermean, puromhyperinvconc, dposomhypermean,
                                              dposomhyperinvshape, pi, purwhypermean, purwhyperinvconc,
                                              poswhypermean, poswhyperinvconc);
        GetModel()->Allocate();
        GetModel()->Update();
        cerr << "-- Reset" << endl;
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
        is >> datapath >> datafile >> treefile >> pi;
        is >> puromhypermean >> puromhyperinvconc;
        is >> dposomhypermean >> dposomhyperinvshape;
        is >> purwhypermean >> purwhyperinvconc;
        is >> poswhypermean >> poswhyperinvconc;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "CODONM2A") {
            model = new CodonM2aModel(datapath, datafile, treefile, pi);
        } else {
            cerr << "-- Error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }
        GetModel()->SetMixtureHyperParameters(puromhypermean, puromhyperinvconc, dposomhypermean,
                                              dposomhyperinvshape, pi, purwhypermean, purwhyperinvconc,
                                              poswhypermean, poswhyperinvconc);
        GetModel()->Allocate();
        GetModel()->FromStream(is);
        GetModel()->Update();
        cerr << size << " points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
        model->Trace(cerr);
    }

    void Save() override {
        ofstream param_os((name + ".param").c_str());
        param_os << GetModelType() << '\n';
        param_os << datapath << '\t' << datafile << '\t' << treefile << '\t' << pi << '\n';
        param_os << puromhypermean << '\t' << puromhyperinvconc << '\n';
        param_os << dposomhypermean << '\t' << dposomhyperinvshape << '\n';
        param_os << purwhypermean << '\t' << purwhyperinvconc << '\n';
        param_os << poswhypermean << '\t' << poswhyperinvconc << '\n';
        param_os << 0 << '\n';
        param_os << every << '\t' << until << '\t' << size << '\n';
        model->ToStream(param_os);
    }

    void SavePoint() override {
        Chain::SavePoint();
        ofstream pos((name + ".sitepp").c_str(), ios_base::app);
        GetModel()->TracePostProb(pos);
    }

    void MakeFiles(int force) override {
        Chain::MakeFiles(force);
        ofstream pos((name + ".sitepp").c_str());
    }
};

int main(int argc, char *argv[]) {

    if (argc == 1)  {

        cerr << "\n";
        cerr << "the M2a model of codeml (Muse and Gaut version)\n";
        cerr << "the omega_i's across sites are a mixture with 3 components\n";
        cerr << " - omega0 < 1, with weight w0\n";
        cerr << " - omega1 = 1, with weight w1\n";
        cerr << " - omega2 > 1, with weight w2\n";
        cerr << '\n';
        cerr << "Here, the model is parameterized as follows:\n";
        cerr << " - omega0 = purom,\n";
        cerr << " - omega1 = 1,\n";
        cerr << " - omega2 = 1 + dposom,\n";
        cerr << "where 0 < purom < 1 and dposom > 0;\n";
        cerr << "purom has a beta prior (hyperparams: puromhypermean and puromhyperinvconc);\n";
        cerr << "dposom has a gamma prior (hyperparams: dposomhypermean and dposomhyperinvshape).\n";
        cerr << '\n';
        cerr << "The weights of the mixture are parameterized as follows:\n";
        cerr << " - w0 = purw * (1 - posw)\n";
        cerr << " - w1 = (1-purw) * (1-posw)\n";
        cerr << " - w2 = posw\n";
        cerr << "where 0<purw<1 and 0<=posw<1;\n";
        cerr << "purw has a beta prior (hyperparams: purwhypermean and purwhyperinvconc);\n";
        cerr << "the prior on posw is a mixture:\n";
        cerr << " - with probability 1-pi, posw = 0\n";
        cerr << " - with probability pi, 0 < posw < 1, in which case is it from a beta prior\n";
        cerr << " (hyperparams: poswhypermean and poswhyperinvconc). Thus, setting pi = 0 imposes a model without positive selection.\n";
        cerr << '\n';
        cerr << "In total, the 9 hyperparameters of the mixture of omegas are as follows:\n";
        cerr << "puromhypermean, puromhyperinvconc, dposomhypermean, dposomhyperinvshape,\n";
        cerr << "purwhypermean, purwhyperinvconc, pi, poswhypermean, poswhyperinvconc.\n";
        cerr << "In a single-gene context, these hyperparameters are fixed;\n";
        cerr << "in a multigene context, they can be either fixed or estimated across genes (see multigenecodonm2a).\n";
        cerr << '\n';
        cerr << "command: codonm2a -d <alignment_list> -t <tree> <chainname>\n";
        cerr << '\n';
        cerr << "program options:\n";
        cerr << "\t-f: force overwrite of already existing chain\n";
        cerr << "\t-x <every> <until>: saving frequency and stopping time (default: every = 1, until = -1)\n";
        cerr << "\t-pi <pi>: specify value for pi (default pi = 0.1)\n";
        cerr << '\n';
        exit(0);
    }

    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        string name = argv[1];
        CodonM2aChain *chain = new CodonM2aChain(name);
        cerr << "chain " << name << " started\n";
        chain->Start();
        cerr << "chain " << name << " stopped\n";
        cerr << chain->GetSize()
             << " points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
        chain->GetModel()->Trace(cerr);
    }

    // new chain
    else {
        string datapath = "./";
        string datafile = "";
        string treefile = "";
        double pi = 0.1;
        double puromhypermean = 0.5;
        double puromhyperinvconc = 0.5;
        double purwhypermean = 0.5;
        double purwhyperinvconc = 0.5;
        double poswhypermean = 0.5;
        double poswhyperinvconc = 0.1;
        double dposomhypermean = 1.0;
        double dposomhyperinvshape = 0.5;
        string name = "";
        int force = 1;
        int every = 1;
        int until = -1;

        try {
            if (argc == 1) {
                throw(0);
            }

            int i = 1;
            while (i < argc) {
                string s = argv[i];

                if (s == "-d") {
                    i++;
                    datafile = argv[i];
                } else if (s == "-p") {
                    i++;
                    datapath = argv[i];
                } else if ((s == "-t") || (s == "-T")) {
                    i++;
                    treefile = argv[i];
                } else if (s == "-f") {
                    force = 1;
                } else if (s == "-hyper")   {
                    i++;
                    ifstream is(argv[i]);
                    string temp;
                    is >> temp >> puromhypermean;
                    is >> temp >> puromhyperinvconc;
                    is >> temp >> dposomhypermean;
                    is >> temp >> dposomhyperinvshape;
                    is >> temp >> purwhypermean;
                    is >> temp >> purwhyperinvconc;
                    is >> temp >> poswhypermean;
                    is >> temp >> poswhyperinvconc;
                } else if (s == "-pi") {
                    i++;
                    pi = atof(argv[i]);
                } else if (s == "-purom") {
                    i++;
                    puromhypermean = atof(argv[i]);
                    i++;
                    puromhyperinvconc = atof(argv[i]);
                } else if (s == "-dposom") {
                    i++;
                    dposomhypermean = atof(argv[i]);
                    i++;
                    dposomhyperinvshape = atof(argv[i]);
                } else if (s == "-purw") {
                    i++;
                    purwhypermean = atof(argv[i]);
                    i++;
                    purwhyperinvconc = atof(argv[i]);
                } else if (s == "-posw") {
                    i++;
                    poswhypermean = atof(argv[i]);
                    i++;
                    poswhyperinvconc = atof(argv[i]);
                } else if ((s == "-x") || (s == "-extract")) {
                    i++;
                    if (i == argc) throw(0);
                    every = atoi(argv[i]);
                    i++;
                    if (i == argc) throw(0);
                    until = atoi(argv[i]);
                } else {
                    if (i != (argc - 1)) {
                        throw(0);
                    }
                    name = argv[i];
                }
                i++;
            }
            if ((datafile == "") || (treefile == "") || (name == "")) {
                throw(0);
            }
        } catch (...) {
            cerr << "error in command\n";
            cerr << '\n';
            exit(0);
        }

        CodonM2aChain *chain = new CodonM2aChain(datapath, datafile, treefile, pi, puromhypermean, puromhyperinvconc, dposomhypermean, dposomhyperinvshape, purwhypermean, purwhyperinvconc, poswhypermean, poswhyperinvconc, every, until, name, force);
        cerr << "chain " << name << " started\n";
        chain->Start();
        cerr << "chain " << name << " stopped\n";
        cerr << chain->GetSize()
             << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
        chain->GetModel()->Trace(cerr);
    }
}
