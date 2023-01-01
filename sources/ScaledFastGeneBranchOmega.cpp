#include <cmath>
#include <fstream>
#include "Chain.hpp"
#include "ScaledFastGeneBranchOmegaModel.hpp"
using namespace std;

/**
 * \brief Chain object for running an MCMC under ScaledFastGeneBranchOmegaModel
 */

class ScaledFastGeneBranchOmegaChain : public Chain {
  private:
    // Chain parameters
    string modeltype;
    string datafile, treefile;
    int syn_devmode, om_devmode;
    int fixed_invshape;

  public:
    //! constructor for a new chain: datafile, treefile, saving frequency, final
    //! chain size, chain name and overwrite flag -- calls New
    ScaledFastGeneBranchOmegaChain(string indatafile, string intreefile, int insyn_devmode, int inom_devmode, int infixed_invshape, int inevery, int inuntil, string inname,
                     int force)
        : modeltype("FASTGENEBRANCHOMEGA"), datafile(indatafile), treefile(intreefile) {
        every = inevery;
        until = inuntil;
        syn_devmode = insyn_devmode;
        om_devmode = inom_devmode;
        fixed_invshape = infixed_invshape;
        name = inname;
        New(force);
    }

    //! constructor for re-opening an already existing chain from file -- calls
    //! Open
    ScaledFastGeneBranchOmegaChain(string filename) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new ScaledFastGeneBranchOmegaModel(datafile, treefile, 
                syn_devmode, om_devmode, fixed_invshape);
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
        is >> datafile >> treefile;
        is >> syn_devmode >> om_devmode >> fixed_invshape;
        int tmp;
        is >> tmp;
        if (tmp)    {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "FASTGENEBRANCHOMEGA") {
            model = new ScaledFastGeneBranchOmegaModel(datafile, treefile,
                    syn_devmode, om_devmode, fixed_invshape);
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
        param_os << datafile << '\t' << treefile << '\n';
        param_os << syn_devmode << '\t' << om_devmode << '\n';
        param_os << fixed_invshape << '\n';
        param_os << 0 << '\n';
        param_os << every << '\t' << until << '\t' << size << '\n';
        model->ToStream(param_os);
    }

    /*
    void SavePoint() override {
        Chain::SavePoint();
        ofstream pos((name + ".branchdsomsuffstat").c_str(), ios_base::app);
        GetModel()->TracePostProb(pos);
    }

    void MakeFiles(int force) override {
        Chain::MakeFiles(force);
        ofstream dsomos((name + ".branchdsomsuffstat").c_str());
    }
    */

    //! return the model, with its derived type (unlike ProbModel::GetModel)
    ScaledFastGeneBranchOmegaModel *GetModel() { return static_cast<ScaledFastGeneBranchOmegaModel *>(model); }

    //! return model type
    string GetModelType() override { return modeltype; }
};

int main(int argc, char *argv[]) {
    string name = "";
    ScaledFastGeneBranchOmegaChain *chain = 0;

    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        chain = new ScaledFastGeneBranchOmegaChain(name);
    }

    // new chain
    else {
        string datafile = "";
        string treefile = "";
        int syn_devmode = 1;
        int om_devmode = 1;
        int fixed_invshape = 0;
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
                } else if (s == "-f") {
                    force = 1;
                } else if (s == "-nonscaled")   {
                    fixed_invshape = 1;
                } else if (s == "-scaled")  {
                    fixed_invshape = 0;
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
            if ((datafile == "") || (treefile == "") || (name == "")) {
                throw(0);
            }
        } catch (...) {
            cerr << "genebranchdnds -d <alignment> -t <tree> <chainname> \n";
            cerr << '\n';
            exit(1);
        }

        chain = new ScaledFastGeneBranchOmegaChain(datafile, treefile, syn_devmode, om_devmode, fixed_invshape, every, until, name, force);
    }

    cerr << "chain " << name << " started\n";
    chain->Start();
    cerr << "chain " << name << " stopped\n";
    cerr << chain->GetSize()
         << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
    chain->GetModel()->Trace(cerr);
}
