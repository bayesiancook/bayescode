#include <cmath>
#include <fstream>
#include "Chain.hpp"
#include "BranchStrandSymmetricCodonModel.hpp"
using namespace std;

/**
 * \brief Chain object for running an MCMC under BranchStrandSymmetricCodonModel
 */

class BranchStrandSymmetricChain : public Chain {
  private:
    // Chain parameters
    string modeltype;
    string datafile, treefile, taxonfile;
    int syn_devmode, om_devmode;

  public:
    //! constructor for a new chain: datafile, treefile, saving frequency, final
    //! chain size, chain name and overwrite flag -- calls New
    BranchStrandSymmetricChain(string indatafile, string intreefile, string intaxonfile, int insyn_devmode, int inom_devmode, int inevery, int inuntil, string inname,
                     int force)
        : modeltype("BRANCHCODON"), datafile(indatafile), treefile(intreefile), taxonfile(intaxonfile) {
        every = inevery;
        until = inuntil;
        syn_devmode = insyn_devmode;
        om_devmode = inom_devmode;
        name = inname;
        New(force);
    }

    //! constructor for re-opening an already existing chain from file -- calls
    //! Open
    BranchStrandSymmetricChain(string filename) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new BranchStrandSymmetricCodonModel(datafile, treefile, taxonfile, syn_devmode, om_devmode);
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
        is >> datafile >> treefile >> taxonfile;
        is >> syn_devmode >> om_devmode;
        int tmp;
        is >> tmp;
        if (tmp)    {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "BRANCHCODON") {
            model = new BranchStrandSymmetricCodonModel(datafile, treefile, taxonfile, syn_devmode, om_devmode);
        } else {
            cerr << "-- Error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }
        model->FromStream(is);
        model->Update();
        cerr << size << " points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
        model->Trace(cerr);
    }

    void Save() override {
        ofstream param_os((name + ".param").c_str());
        param_os << GetModelType() << '\n';
        param_os << datafile << '\t' << treefile << '\t' << taxonfile << '\n';
        param_os << syn_devmode << '\t' << om_devmode << '\n';
        param_os << 0 << '\n';
        param_os << every << '\t' << until << '\t' << size << '\n';
        model->ToStream(param_os);
    }

    //! return the model, with its derived type (unlike ProbModel::GetModel)
    BranchStrandSymmetricCodonModel *GetModel() { return static_cast<BranchStrandSymmetricCodonModel *>(model); }

    //! return model type
    string GetModelType() override { return modeltype; }
};

int main(int argc, char *argv[]) {
    string name = "";
    BranchStrandSymmetricChain *chain = 0;

    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        chain = new BranchStrandSymmetricChain(name);
    }

    // new chain
    else {
        string datafile = "";
        string treefile = "";
        string taxonfile = "";
        int syn_devmode = 1;
        int om_devmode = 1;
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

                if (s == "-d") {
                    i++;
                    datafile = argv[i];
                } else if ((s == "-t") || (s == "-T")) {
                    i++;
                    treefile = argv[i];
                } else if (s == "-tax") {
                    i++;
                    taxonfile = argv[i];
                } else if (s == "-f") {
                    force = 1;
                } else if (s == "-gamdev") {
                    syn_devmode = 1;
                    om_devmode = 1;
                } else if (s == "-mixdev") {
                    syn_devmode = 2;
                    om_devmode = 2;
                } else if (s == "-nodev")   {
                    syn_devmode = 0;
                    om_devmode = 0;
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
            if ((datafile == "") || (treefile == "") || (taxonfile == "") || (name == "")) {
                throw(0);
            }
        } catch (...) {
            cerr << "genebranchdnds -d <alignment> -t <tree> <chainname> \n";
            cerr << '\n';
            exit(1);
        }

        chain = new BranchStrandSymmetricChain(datafile, treefile, taxonfile, syn_devmode, om_devmode, every, until, name, force);
    }

    cerr << "chain " << name << " started\n";
    chain->Start();
    cerr << "chain " << name << " stopped\n";
    cerr << chain->GetSize()
         << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
    chain->GetModel()->Trace(cerr);
}
