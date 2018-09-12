#include <cmath>
#include <fstream>
#include "Chain.hpp"
#include "DiffSelModel.hpp"
using namespace std;

/**
 * \brief Chain object for running an MCMC under the DiffSelModel
 */

class DiffSelChain : public Chain {
  private:
    // Chain parameters
    string modeltype, datafile, treefile;
    int codonmodel, ncond, nlevel, fixglob, fixvar;

  public:
    DiffSelModel *GetModel() { return static_cast<DiffSelModel *>(model); }

    string GetModelType() override { return modeltype; }

    DiffSelChain(string indata, string intree, int inncond, int innlevel, int inevery, int inuntil,
        int insaveall, int infixglob, int infixvar, int incodonmodel, string inname, int force)
        : modeltype("DIFFSEL"),
          datafile(indata),
          treefile(intree),
          codonmodel(incodonmodel),
          ncond(inncond),
          nlevel(innlevel),
          fixglob(infixglob),
          fixvar(infixvar) {
        every = inevery;
        until = inuntil;
        saveall = insaveall;
        name = inname;
        New(force);
    }

    DiffSelChain(string filename) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new DiffSelModel(datafile, treefile, ncond, nlevel, fixglob, fixvar, codonmodel);
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
        is >> datafile >> treefile >> ncond >> nlevel;
        is >> fixglob >> fixvar >> codonmodel;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> saveall >> size;

        if (modeltype == "DIFFSEL") {
            model =
                new DiffSelModel(datafile, treefile, ncond, nlevel, fixglob, fixvar, codonmodel);
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
        param_os << datafile << '\t' << treefile << '\t' << ncond << '\t' << nlevel << '\n';
        param_os << fixglob << '\t' << fixvar << '\t' << codonmodel << '\n';
        param_os << 0 << '\n';
        param_os << every << '\t' << until << '\t' << saveall << '\t' << size << '\n';

        model->ToStream(param_os);
    }

    void SavePoint() override {
        Chain::SavePoint();
        ofstream os((name + ".baseline").c_str(), ios_base::app);
        GetModel()->TraceBaseline(os);
        for (int k = 1; k < ncond; k++) {
            ostringstream s;
            s << name << "_" << k;
            ofstream os((s.str() + ".delta").c_str(), ios_base::app);
            GetModel()->TraceDelta(k, os);
        }
    }

    void MakeFiles(int force) override {
        Chain::MakeFiles(force);
        ofstream os((name + ".baseline").c_str());
        for (int k = 1; k < ncond; k++) {
            ostringstream s;
            s << name << "_" << k;
            ofstream os((s.str() + ".delta").c_str());
        }
    }
};

int main(int argc, char *argv[]) {
    string name = "";
    DiffSelChain *chain = 0;

    if (argc == 1) {
        cerr << "command: diffsel -d <alignment> -t <treefile> <chainname>\n";
        cerr << '\n';
        cerr << "program options:\n";
        cerr << "\t-f: force overwrite of already existing chain\n";
        cerr << "\t-x <every> <until>: saving frequency and stopping time (default: every = 1, "
                "until = -1)\n";
        cerr << "\t-ncond <ncond>:  specify number of conditions\n";
        cerr << '\n';
        exit(0);
    }

    // this is an already existing chain on the disk; reopen and restart
    if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        chain = new DiffSelChain(name);
    }

    // this is a new chain
    else {
        string datafile = "";
        string treefile = "";
        int ncond = 2;
        int nlevel = 2;
        int fixglob = 1;
        int fixvar = 1;
        int codonmodel = 1;

        int every = 1;
        int until = -1;
        int saveall = 1;
        int force = 0;

        try {
            if (argc == 1) { throw(0); }

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
                } else if (s == "+s") {
                    saveall = 1;
                } else if (s == "-s") {
                    saveall = 0;
                } else if (s == "-ncond") {
                    i++;
                    ncond = atoi(argv[i]);
                } else if (s == "-nlevel") {
                    i++;
                    nlevel = atoi(argv[i]);
                } else if ((s == "-x") || (s == "-extract")) {
                    i++;
                    if (i == argc) throw(0);
                    every = atoi(argv[i]);
                    i++;
                    if (i == argc) throw(0);
                    until = atoi(argv[i]);
                } else {
                    if (i != (argc - 1)) { throw(0); }
                    name = argv[i];
                }
                i++;
            }
        } catch (...) {
            cerr << "error in command\n";
            exit(0);
        }
        chain = new DiffSelChain(datafile, treefile, ncond, nlevel, every, until, saveall, fixglob,
            fixvar, codonmodel, name, force);
    }

    cerr << "chain " << name << " started\n";
    chain->Start();
    cerr << "chain " << name << " stopped\n";
    cerr << chain->GetSize()
         << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
    chain->GetModel()->Trace(cerr);
}
