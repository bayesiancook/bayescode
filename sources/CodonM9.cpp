
#include <cmath>
#include <fstream>
#include "Chain.hpp"
#include "CodonM9Model.hpp"
using namespace std;

/**
 * \brief Chain object for running an MCMC under CodonM9Model
 */

class CodonM9Chain : public Chain {
  private:
    // Chain parameters
    string modeltype, datafile, treefile;
    double pi;
    double purifmeanhypermean, purifmeanhyperinvconc;
    double purifinvconchypermean, purifinvconchyperinvshape;
    vector<double> purifweighthypercenter;
    double purifweighthyperinvconc;
    double posmeanhypermean, posmeanhyperinvshape;
    double posinvshapehypermean, posinvshapehyperinvshape;
    double poswhypermean, poswhyperinvconc;

  public:
    //! return the model, with its derived type (unlike ProbModel::GetModel)
    CodonM9Model *GetModel() { return static_cast<CodonM9Model *>(model); }

    string GetModelType() override { return modeltype; }

    //! constructor for a new chain: datafile, treefile, pi (fraction of sites
    //! under positive selection) saving frequency, final chain size, chain name
    //! and overwrite flag -- calls New
    CodonM9Chain(string indatafile, string intreefile,
            double inpi, 
            double inpurifmeanhypermean, double inpurifmeanhyperinvconc,
            double inpurifinvconchypermean, double inpurifinvconchyperinvshape,
            vector<double> inpurifweighthypercenter, double inpurifweighthyperinvconc,
            double inposmeanhypermean, double inposmeanhyperinvshape,
            double inposinvshapehypermean, double inposinvshapehyperinvshape,
            double inposwhypermean, double inposwhyperinvconc,
            int inevery, int inuntil,
            string inname, int force)
        : modeltype("CODONM9"), datafile(indatafile), treefile(intreefile),
            pi(inpi),
            purifmeanhypermean(inpurifmeanhypermean),
            purifmeanhyperinvconc(inpurifmeanhyperinvconc),
            purifinvconchypermean(inpurifinvconchypermean),
            purifinvconchyperinvshape(inpurifinvconchyperinvshape),
            purifweighthypercenter(inpurifweighthypercenter),
            purifweighthyperinvconc(inpurifweighthyperinvconc),
            posmeanhypermean(inposmeanhypermean),
            posmeanhyperinvshape(inposmeanhyperinvshape),
            posinvshapehypermean(inposinvshapehypermean),
            posinvshapehyperinvshape(inposinvshapehyperinvshape),
            poswhypermean(inposwhypermean),
            poswhyperinvconc(inposwhyperinvconc)    {

            every = inevery;
            until = inuntil;
            name = inname;
            New(force);
    }

    //! constructor for re-opening an already existing chain from file -- calls
    //! Open
    CodonM9Chain(string filename) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new CodonM9Model(datafile, treefile, pi);

        GetModel()->SetMixtureHyperParameters(
                pi,
                purifmeanhypermean, purifmeanhyperinvconc,
                purifinvconchypermean, purifinvconchyperinvshape,
                purifweighthypercenter, purifweighthyperinvconc,
                posmeanhypermean, posmeanhyperinvshape,
                posinvshapehypermean, posinvshapehyperinvshape,
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
        is >> datafile >> treefile >> pi;
        is >> purifmeanhypermean >> purifmeanhyperinvconc;
        is >> purifinvconchypermean >> purifinvconchyperinvshape;
        is >> purifweighthypercenter >> purifweighthyperinvconc;
        is >> posmeanhypermean >> posmeanhyperinvshape;
        is >> posinvshapehypermean >> posinvshapehyperinvshape;
        is >> poswhypermean >> poswhyperinvconc;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "CODONM9") {
            model = new CodonM9Model(datafile, treefile, pi);
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
                posmeanhypermean, posmeanhyperinvshape,
                posinvshapehypermean, posinvshapehyperinvshape,
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
        param_os << datafile << '\t' << treefile << '\t' << pi << '\n';
        param_os << purifmeanhypermean << '\t' << purifmeanhyperinvconc << '\n';
        param_os << purifinvconchypermean << '\t' << purifinvconchyperinvshape << '\n';
        param_os << purifweighthypercenter << '\t' << purifweighthyperinvconc << '\n';
        param_os << posmeanhypermean << '\t' << posmeanhyperinvshape << '\n';
        param_os << posinvshapehypermean << '\t' << posinvshapehyperinvshape << '\n';
        param_os << poswhypermean << '\t' << poswhyperinvconc << '\n';
        param_os << 0 << '\n';
        param_os << every << '\t' << until << '\t' << size << '\n';
        model->ToStream(param_os);
    }

    void SavePoint() override {
        Chain::SavePoint();
        ofstream os((name + ".siteom").c_str(), ios_base::app);
        GetModel()->TraceOmega(os);
    }

    void MakeFiles(int force) override {
        Chain::MakeFiles(force);
        ofstream os((name + ".siteom").c_str());
    }
};

int main(int argc, char *argv[]) {

    if (argc == 1)  {

        cerr << "\n";
        cerr << "the M9 model of codeml (Muse and Gaut version)\n";
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
        CodonM9Chain *chain = new CodonM9Chain(name);
        cerr << "chain " << name << " started\n";
        chain->Start();
        cerr << "chain " << name << " stopped\n";
        cerr << chain->GetSize()
             << " points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
        chain->GetModel()->Trace(cerr);
    }

    // new chain
    else {
        string datafile = "";
        string treefile = "";

        double pi = 0.1;

        double purifmeanhypermean = 0.5;
        double purifmeanhyperinvconc = 0.5;
        double purifinvconchypermean = 1.0;
        double purifinvconchyperinvshape = 1.0;
        vector<double> purifweighthypercenter(3,1.0/3);
        double purifweighthyperinvconc = 1.0;

        double poswhypermean = 0.5;
        double poswhyperinvconc = 0.1;

        double posmeanhypermean = 1.0;
        double posmeanhyperinvshape = 1.0;
        double posinvshapehypermean = 1.0;
        double posinvshapehyperinvshape = 1.0;

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
                } else if ((s == "-t") || (s == "-T")) {
                    i++;
                    treefile = argv[i];
                } else if (s == "-f") {
                    force = 1;
                } else if (s == "-pi") {
                    i++;
                    pi = atof(argv[i]);
                } else if (s == "-dposom") {
                    i++;
                    posmeanhypermean = atof(argv[i]);
                    i++;
                    posmeanhyperinvshape = atof(argv[i]);
                    i++;
                    posinvshapehypermean = atof(argv[i]);
                    i++;
                    posinvshapehyperinvshape = atof(argv[i]);
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

        CodonM9Chain *chain = new CodonM9Chain(datafile, treefile, pi, purifmeanhypermean, purifmeanhyperinvconc, purifinvconchypermean, purifinvconchyperinvshape, purifweighthypercenter, purifweighthyperinvconc, posmeanhypermean, posmeanhyperinvshape, posinvshapehypermean, posinvshapehyperinvshape, poswhypermean, poswhyperinvconc, every, until, name, force);
        cerr << "chain " << name << " started\n";
        chain->Start();
        cerr << "chain " << name << " stopped\n";
        cerr << chain->GetSize()
             << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
        chain->GetModel()->Trace(cerr);
    }
}
