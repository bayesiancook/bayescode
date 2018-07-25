#include <cmath>
#include <fstream>
#include "Chain.hpp"
#include "SparseAASubSelModel.hpp"
using namespace std;

/**
 * \brief Chain object for running an MCMC under SparseAASubSelModel
 */

class SparseAASubSelChain : public Chain {
  private:
    // Chain parameters
    string modeltype, datafile, treefile;
    // -1: estimated
    // 1 : no mask
    int ratemode;
    int profilemode;
    double epsilon;

  public:
    SparseAASubSelModel *GetModel() { return static_cast<SparseAASubSelModel *>(model); }

    string GetModelType() override { return modeltype; }

    SparseAASubSelChain(string indatafile, string intreefile, int inratemode, int inprofilemode,
                        double inepsilon, int inevery, int inuntil, string inname, int force)
        : modeltype("SPARSEAASUBSEL"),
          datafile(indatafile),
          treefile(intreefile),
          ratemode(inratemode),
          profilemode(inprofilemode),
          epsilon(inepsilon) {
        every = inevery;
        until = inuntil;
        name = inname;
        New(force);
    }

    SparseAASubSelChain(string filename) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new SparseAASubSelModel(datafile, treefile, ratemode, profilemode, epsilon);
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
        is >> datafile >> treefile;
        is >> ratemode;
        is >> profilemode;
        is >> epsilon;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "SPARSEAASUBSEL") {
            model = new SparseAASubSelModel(datafile, treefile, ratemode, profilemode, epsilon);
        } else {
            cerr << "-- Error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }
        GetModel()->Allocate();
        model->FromStream(is);
        GetModel()->Update();
        cerr << size << " points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
        model->Trace(cerr);
    }

    void Save() override {
        ofstream param_os((name + ".param").c_str());
        param_os << GetModelType() << '\n';
        param_os << datafile << '\t' << treefile << '\n';
        param_os << ratemode << '\n';
        param_os << profilemode << '\n';
        param_os << epsilon << '\n';
        param_os << 0 << '\n';
        param_os << every << '\t' << until << '\t' << size << '\n';
        model->ToStream(param_os);
    }
};

int main(int argc, char *argv[]) {
    string name = "";
    SparseAASubSelChain *chain = 0;

    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        chain = new SparseAASubSelChain(name);
    }

    // new chain
    else {
        string datafile = "";
        string treefile = "";
        double epsilon = -1;
        int ratemode = 0;
        int profilemode = 0;
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
                } else if ((s == "-eps") || (s == "-epsilon")) {
                    i++;
                    string tmp = argv[i];
                    if (tmp == "free") {
                        epsilon = -1;
                    } else {
                        epsilon = atof(argv[i]);
                    }
                } else if (s == "-uni") {
                    ratemode = 0;
                } else if (s == "-cgam") {
                    ratemode = 1;
                } else if (s == "-flat") {
                    profilemode = 3;
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
            cerr << "sparsecatgtr -d <alignment> -t <tree> <chainname> \n";
            cerr << '\n';
            exit(1);
        }

        chain = new SparseAASubSelChain(datafile, treefile, ratemode, profilemode, epsilon, every,
                                        until, name, force);
    }

    cerr << "chain " << name << " started\n";
    chain->Start();
    cerr << "chain " << name << " stopped\n";
    cerr << chain->GetSize()
         << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
    chain->GetModel()->Trace(cerr);
}
