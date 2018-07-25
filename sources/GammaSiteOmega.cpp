#include <cmath>
#include <fstream>
#include "Chain.hpp"
#include "GammaSiteOmegaModel.hpp"
using namespace std;

class GammaSiteOmegaChain : public Chain {
  private:
    // Chain parameters
    string modeltype, datafile, treefile;

  public:
    GammaSiteOmegaModel *GetModel() { return static_cast<GammaSiteOmegaModel *>(model); }

    string GetModelType() override { return modeltype; }

    GammaSiteOmegaChain(string indatafile, string intreefile, int inevery, int inuntil,
                        string inname, int force)
        : modeltype("GAMMASITEOMEGA"), datafile(indatafile), treefile(intreefile) {
        every = inevery;
        until = inuntil;
        name = inname;
        New(force);
    }

    GammaSiteOmegaChain(string filename) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new GammaSiteOmegaModel(datafile, treefile);
        GetModel()->Allocate();
        GetModel()->Unfold();
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
        is >> datafile >> treefile;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "GAMMASITEOMEGA") {
            model = new GammaSiteOmegaModel(datafile, treefile);
        } else {
            cerr << "-- Error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }
        GetModel()->Allocate();
        model->FromStream(is);
        model->Update();
        GetModel()->Unfold();
        cerr << size << " points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
        model->Trace(cerr);
    }

    void Save() override {
        ofstream param_os((name + ".param").c_str());
        param_os << GetModelType() << '\n';
        param_os << datafile << '\t' << treefile << '\n';
        param_os << 0 << '\n';
        param_os << every << '\t' << until << '\t' << size << '\n';
        model->ToStream(param_os);

        ofstream pos((name + ".siteom").c_str(), ios_base::app);
        GetModel()->TraceSiteOmega(pos);
    }

    void MakeFiles(int force) override {
        Chain::MakeFiles(force);
        ofstream pos((name + ".sitepp").c_str());
    }
};

int main(int argc, char *argv[]) {
    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        string name = argv[1];
        GammaSiteOmegaChain *chain = new GammaSiteOmegaChain(name);
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
            cerr << "globom -d <alignment> -t <tree> <chainname> \n";
            cerr << '\n';
            exit(1);
        }

        GammaSiteOmegaChain *chain =
            new GammaSiteOmegaChain(datafile, treefile, every, until, name, force);
        cerr << "chain " << name << " started\n";
        chain->Start();
        cerr << "chain " << name << " stopped\n";
        cerr << chain->GetSize()
             << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
        chain->GetModel()->Trace(cerr);
    }
}
