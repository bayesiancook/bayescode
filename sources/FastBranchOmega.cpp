#include <cmath>
#include <fstream>
#include "Chain.hpp"
#include "FastBranchOmegaModel.hpp"
using namespace std;

/**
 * \brief Chain object for running an MCMC under FastBranchOmegaModel
 */

class FastBranchOmegaChain : public Chain {
  private:
    // Chain parameters
    string modeltype;
    string taxonfile, treefile, dsomsuffstatfile;

  public:
    //! constructor for a new chain: datafile, treefile, saving frequency, final
    //! chain size, chain name and overwrite flag -- calls New
    FastBranchOmegaChain(string intaxonfile, string intreefile, string indsomsuffstatfile, int inevery, int inuntil, string inname,
                     int force)
        : modeltype("FASTBRANCHOMEGA"), taxonfile(intaxonfile), treefile(intreefile), dsomsuffstatfile(indsomsuffstatfile) {
        every = inevery;
        until = inuntil;
        name = inname;
        New(force);
    }

    //! constructor for re-opening an already existing chain from file -- calls
    //! Open
    FastBranchOmegaChain(string filename) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new FastBranchOmegaModel(taxonfile, treefile, dsomsuffstatfile);
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
        is >> taxonfile >> treefile >> dsomsuffstatfile;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "FASTBRANCHOMEGA") {
            model = new FastBranchOmegaModel(taxonfile, treefile, dsomsuffstatfile);
        } else {
            cerr << "-- Error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }
        GetModel()->Allocate();
        model->FromStream(is);
        model->Update();
        cerr << size << " points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
        model->Trace(cerr);
    }

    void Save() override {
        ofstream param_os((name + ".param").c_str());
        param_os << GetModelType() << '\n';
        param_os << taxonfile << '\t' << treefile << '\t' << dsomsuffstatfile << '\n';
        param_os << 0 << '\n';
        param_os << every << '\t' << until << '\t' << size << '\n';
        model->ToStream(param_os);
    }

    void SavePoint() override {
        Chain::SavePoint();
        ofstream pos((name + ".branchomega").c_str(), ios_base::app);
        GetModel()->TraceBranchOmega(pos);
    }

    void MakeFiles(int force) override {
        Chain::MakeFiles(force);
        ofstream pos((name + ".branchomega").c_str());
    }

    //! return the model, with its derived type (unlike ProbModel::GetModel)
    FastBranchOmegaModel *GetModel() { return static_cast<FastBranchOmegaModel *>(model); }

    //! return model type
    string GetModelType() override { return modeltype; }
};

int main(int argc, char *argv[]) {
    string name = "";
    FastBranchOmegaChain *chain = 0;

    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        chain = new FastBranchOmegaChain(name);
    }

    // new chain
    else {
        string taxonfile = "";
        string treefile = "";
        string dsomsuffstatfile = "";
        name = "";
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

                if (s == "-tax") {
                    i++;
                    taxonfile = argv[i];
                } else if ((s == "-t") || (s == "-T")) {
                    i++;
                    treefile = argv[i];
                } else if (s == "-f") {
                    force = 1;
                } else if (s == "-ss") {
                    i++;
                    dsomsuffstatfile = argv[i];
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
            if ((taxonfile == "") || (treefile == "") || (name == "")) {
                throw(0);
            }
        } catch (...) {
            cerr << "globom -d <alignment> -t <tree> <chainname> \n";
            cerr << '\n';
            exit(1);
        }

        chain = new FastBranchOmegaChain(taxonfile, treefile, dsomsuffstatfile, every, until, name, force);
    }

    cerr << "chain " << name << " started\n";
    chain->Start();
    cerr << "chain " << name << " stopped\n";
    cerr << chain->GetSize()
         << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
    chain->GetModel()->Trace(cerr);
}
