
#include <cmath>
#include <fstream>
#include "Chain.hpp"
#include "CodonM8Model.hpp"
using namespace std;

/**
 * \brief Chain object for running an MCMC under CodonM8Model
 */

class CodonM8Chain : public Chain {
  private:
    // Chain parameters
    string modeltype, datapath, datafile, treefile;
    int ncat;
    double pi;
    double purifmeanhypermean, purifmeanhyperinvconc;
    double purifinvconchypermean, purifinvconchyperinvshape;
    vector<double> purifweighthypercenter;
    double purifweighthyperinvconc;
    double dposomhypermean, dposomhyperinvshape;
    double poswhypermean, poswhyperinvconc;

  public:
    //! return the model, with its derived type (unlike ProbModel::GetModel)
    CodonM8Model *GetModel() { return static_cast<CodonM8Model *>(model); }

    string GetModelType() override { return modeltype; }

    //! constructor for a new chain: datapath datafile, treefile, pi (fraction of sites
    //! under positive selection) saving frequency, final chain size, chain name
    //! and overwrite flag -- calls New
    CodonM8Chain(string indatapath, string indatafile, string intreefile,
            int inncat, double inpi, 
            double inpurifmeanhypermean, double inpurifmeanhyperinvconc,
            double inpurifinvconchypermean, double inpurifinvconchyperinvshape,
            vector<double> inpurifweighthypercenter, double inpurifweighthyperinvconc,
            double indposomhypermean, double indposomhyperinvshape,
            double inposwhypermean, double inposwhyperinvconc,
            int inevery, int inuntil,
            string inname, int force)
        : modeltype("CODONM8"), datapath(indatapath), datafile(indatafile), treefile(intreefile),
            ncat(inncat),
            pi(inpi),
            purifmeanhypermean(inpurifmeanhypermean),
            purifmeanhyperinvconc(inpurifmeanhyperinvconc),
            purifinvconchypermean(inpurifinvconchypermean),
            purifinvconchyperinvshape(inpurifinvconchyperinvshape),
            purifweighthypercenter(inpurifweighthypercenter),
            purifweighthyperinvconc(inpurifweighthyperinvconc),
            dposomhypermean(indposomhypermean),
            dposomhyperinvshape(indposomhyperinvshape),
            poswhypermean(inposwhypermean),
            poswhyperinvconc(inposwhyperinvconc)    {

            every = inevery;
            until = inuntil;
            name = inname;
            New(force);
    }

    //! constructor for re-opening an already existing chain from file -- calls
    //! Open
    CodonM8Chain(string filename) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new CodonM8Model(datapath, datafile, treefile, ncat);

        GetModel()->SetMixtureHyperParameters(
                pi,
                purifmeanhypermean, purifmeanhyperinvconc,
                purifinvconchypermean, purifinvconchyperinvshape,
                purifweighthypercenter, purifweighthyperinvconc,
                dposomhypermean, dposomhyperinvshape,
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
        is >> datapath >> datafile >> treefile >> ncat >> pi;
        is >> purifmeanhypermean >> purifmeanhyperinvconc;
        is >> purifinvconchypermean >> purifinvconchyperinvshape;
        is >> purifweighthypercenter >> purifweighthyperinvconc;
        is >> dposomhypermean >> dposomhyperinvshape;
        is >> poswhypermean >> poswhyperinvconc;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "CODONM8") {
            model = new CodonM8Model(datapath, datafile, treefile, ncat);
        } else {
            cerr << "-- Error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }

        GetModel()->SetMixtureHyperParameters(
                pi,
                purifmeanhypermean, purifmeanhyperinvconc,
                purifinvconchypermean, purifinvconchyperinvshape,
                purifweighthypercenter, purifweighthyperinvconc,
                dposomhypermean, dposomhyperinvshape,
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
        param_os << datapath << '\t' << datafile << '\t' << treefile << '\t' << ncat << '\t' << pi << '\n';
        param_os << purifmeanhypermean << '\t' << purifmeanhyperinvconc << '\n';
        param_os << purifinvconchypermean << '\t' << purifinvconchyperinvshape << '\n';
        param_os << purifweighthypercenter << '\t' << purifweighthyperinvconc << '\n';
        param_os << dposomhypermean << '\t' << dposomhyperinvshape << '\n';
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
        cerr << "the M8 model of codeml (Muse and Gaut version)\n";
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
        CodonM8Chain *chain = new CodonM8Chain(name);
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

        int ncat = 4;
        double pi = 0.1;

        double purifmeanhypermean = 0.5;
        double purifmeanhyperinvconc = 0.5;
        double purifinvconchypermean = 1.0;
        double purifinvconchyperinvshape = 1.0;
        vector<double> purifweighthypercenter(3,1.0/3);
        double purifweighthyperinvconc = 1.0;

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
                } else if (s == "-ncat")    {
                    i++;
                    ncat = atoi(argv[i]);
                } else if (s == "-pi") {
                    i++;
                    pi = atof(argv[i]);
                } else if (s == "-dposom") {
                    i++;
                    dposomhypermean = atof(argv[i]);
                    i++;
                    dposomhyperinvshape = atof(argv[i]);
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

        CodonM8Chain *chain = new CodonM8Chain(datapath, datafile, treefile, ncat, pi, purifmeanhypermean, purifmeanhyperinvconc, purifinvconchypermean, purifinvconchyperinvshape, purifweighthypercenter, purifweighthyperinvconc, dposomhypermean, dposomhyperinvshape, poswhypermean, poswhyperinvconc, every, until, name, force);
        cerr << "chain " << name << " started\n";
        chain->Start();
        cerr << "chain " << name << " stopped\n";
        cerr << chain->GetSize()
             << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
        chain->GetModel()->Trace(cerr);
    }
}
