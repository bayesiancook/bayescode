#include <cmath>
#include <fstream>
#include "Chain.hpp"
#include "GeneBranchNucModel.hpp"
using namespace std;

/**
 * \brief Chain object for running an MCMC under GeneBranchNucModel
 */

class GeneBranchStrandSymmetricChain : public Chain {
  private:
    // Chain parameters
    string modeltype;
    string datafile, treefile, taxonfile;
    int syn_devmode, nuc_devmode;

  public:
    //! constructor for a new chain: datafile, treefile, saving frequency, final
    //! chain size, chain name and overwrite flag -- calls New
    GeneBranchStrandSymmetricChain(string indatafile, string intreefile, string intaxonfile, int insyn_devmode, int innuc_devmode, int inevery, int inuntil, string inname,
                     int force)
        : modeltype("GENEBRANCHNUC"), datafile(indatafile), treefile(intreefile), taxonfile(intaxonfile) {
        every = inevery;
        until = inuntil;
        syn_devmode = insyn_devmode;
        nuc_devmode = innuc_devmode;
        name = inname;
        New(force);
    }

    //! constructor for re-opening an already existing chain from file -- calls
    //! Open
    GeneBranchStrandSymmetricChain(string filename) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new GeneBranchNucModel(datafile, treefile, taxonfile, syn_devmode, nuc_devmode);
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
        is >> syn_devmode >> nuc_devmode;
        int tmp;
        is >> tmp;
        if (tmp)    {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "GENEBRANCHNUC") {
            model = new GeneBranchNucModel(datafile, treefile, taxonfile, syn_devmode, nuc_devmode);
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
        param_os << syn_devmode << '\t' << nuc_devmode << '\n';
        param_os << 0 << '\n';
        param_os << every << '\t' << until << '\t' << size << '\n';
        model->ToStream(param_os);
    }

    //! return the model, with its derived type (unlike ProbModel::GetModel)
    GeneBranchNucModel *GetModel() { return static_cast<GeneBranchNucModel *>(model); }

    //! return model type
    string GetModelType() override { return modeltype; }
};

int main(int argc, char *argv[]) {
    string name = "";
    GeneBranchStrandSymmetricChain *chain = 0;

    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        chain = new GeneBranchStrandSymmetricChain(name);
    }

    // new chain
    else {
        string datafile = "";
        string treefile = "";
        string taxonfile = "";
        int syn_devmode = 1;
        int nuc_devmode = 1;
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
                } else if (s == "-tax")   {
                    i++;
                    taxonfile = argv[i];
                } else if (s == "-f") {
                    force = 1;
                } else if (s == "-gamdev") {
                    syn_devmode = 1;
                } else if (s == "-mixdev") {
                    syn_devmode = 2;
                } else if (s == "-nodev")   {
                    syn_devmode = 0;
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
            cerr << "genebranchdnds -d <alignment> -t <tree> -tax <taxfile> <chainname> \n";
            cerr << '\n';
            exit(1);
        }

        chain = new GeneBranchStrandSymmetricChain(datafile, treefile, taxonfile, syn_devmode, nuc_devmode, every, until, name, force);
    }

    cerr << "chain " << name << " started\n";
    chain->Start();
    cerr << "chain " << name << " stopped\n";
    cerr << chain->GetSize()
         << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
    chain->GetModel()->Trace(cerr);
}
