#include <cmath>
#include <fstream>
#include "Chain.hpp"
#include "SingleOmegaModel.hpp"
using namespace std;

/**
 * \brief Chain object for running an MCMC under SingleOmegaModel
 */

class SingleOmegaChain : public Start {
  private:
    // Chain parameters
    string datafile, treefile;
    string name;  // name for files
    Chain* chain;

  public:
    //! constructor for a new chain: datafile, treefile, saving frequency, final chain size, chain name and overwrite flag --
    //! calls New
    SingleOmegaChain(string indatafile, string intreefile, int inevery, int inuntil, string inname, int force)
        : modeltype("SINGLEOMEGA"), datafile(indatafile), treefile(intreefile) {
        every = inevery;
        until = inuntil;
        name = inname;
        New(force);
    }

    //! constructor for re-opening an already existing chain from file -- calls Open
    SingleOmegaChain(string filename) {
        name = filename;
        Open();
        Save();
    }

    void start() override { chain->start(); }

    void New(int force) override {
        model = new SingleOmegaModel(datafile, treefile);
        GetModel()->DeclareModel();  // instead of Allocate
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

        if (modeltype == "SINGLEOMEGA") {
            model = new SingleOmegaModel(datafile, treefile);
        } else {
            cerr << "-- Error when opening file " << name << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }
        // GetModel()->Allocate();
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
    }

    //! return the model, with its derived type (unlike ProbModel::GetModel)
    SingleOmegaModel* GetModel() { return static_cast<SingleOmegaModel*>(model); }

    //! return model type
    string GetModelType() override { return modeltype; }
};

int main(int argc, char* argv[]) {
    Model model;
    string name;

    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        cerr << "ERROR: restarting chain not implemented in component version!\n";
        exit(1);
        // SingleOmegaChain* chain = new SingleOmegaChain(name);
        // cerr << "chain " << name << " started\n";
        // chain->start();
        // cerr << "chain " << name << " stopped\n";
        // cerr << chain->GetSize() << " points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
        // chain->GetModel()->Trace(cerr);
    }

    // new chain
    else {
        string datafile = "";
        string treefile = "";
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

        model.component<SingleOmegaModel>("model", datafile, treefile);  // remove files form this constructor

        model.component<Chain>("chain")
            .connect<Use<ProbModel>>("model")
            .connect<Use<TraceFile>>("chainfile", "chainfile")
            .connect<Use<TraceFile>>("fitnessfile", "fitnessfile")
            .connect<Use<TraceFile>>("monitorfile", "monitorfile")
            .connect<Use<TraceFile>>("paramfile", "paramfile")
            .connect<Use<TraceFile>>("runfile", "runfile")
            .connect<Use<TraceFile>>("tracefile", "tracefile");

        model.component<TraceFile>("chainfile", name + ".chain");
        model.component<TraceFile>("fitnessfile", name + ".fitness");
        model.component<TraceFile>("monitorfile", name + ".monitor");
        model.component<TraceFile>("paramfile", name + ".param");
        model.component<TraceFile>("runfile", name + ".run");
        model.component<TraceFile>("tracefile", name + ".trace");
    }
    Assembly assembly(model);  // instantiating assembly!

    // SingleOmegaChain* chain = new SingleOmegaChain(datafile, treefile, every, until, name, force);
    cerr << "-- Chain " << name << " starting...\n";
    // chain->start();
    assembly.call("SingleOmegaDriver", "start");
    cerr << "-- Chain " << name << " stopped\n";
    // cerr << chain->GetSize() << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
    // chain->GetModel()->Trace(cerr);
}
