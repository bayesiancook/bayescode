#include <cmath>
#include <fstream>
#include "DiscSelACModel.hpp"
#include "Chain.hpp"
using namespace std;

/**
 * \brief Chain object for running an MCMC under DiscSelACModel
 */

class SelACChain : public Chain {
  private:
    // Chain parameters
    string modeltype, datafile, treefile;
    int aadistmodel;
    int Gcat;
    double xmin, xmax;

  public:
    DiscSelACModel *GetModel() { return static_cast<DiscSelACModel *>(model); }

    string GetModelType() override { return modeltype; }

    SelACChain(string indatafile, string intreefile, int inaadistmodel,
                            int inGcat, double inxmin, double inxmax, int inevery,
                            int inuntil, string inname, int force)
        : modeltype("DISCSELAC"),
          datafile(indatafile),
          treefile(intreefile),
          aadistmodel(inaadistmodel),
          Gcat(inGcat), xmin(inxmin), xmax(inxmax)  {
        every = inevery;
        until = inuntil;
        name = inname;
        New(force);
    }

    SelACChain(string filename) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        cerr << "new model\n";
        model =
            new DiscSelACModel(datafile, treefile, aadistmodel, Gcat, xmin, xmax);
        cerr << "allocate\n";
        GetModel()->Allocate();
        cerr << "update\n";
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
        is >> aadistmodel;
        is >> Gcat;
        is >> xmin >> xmax;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "DISCSELAC") {
            model = new DiscSelACModel(datafile, treefile, aadistmodel, Gcat, xmin, xmax);
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

    void MakeFiles(int force) override {
        Chain::MakeFiles(force);
        ofstream os((name + ".gweights").c_str());
    }

    void SavePoint() override {
        Chain::SavePoint();
        ofstream os((name + ".gweights").c_str(), ios_base::app);
        GetModel()->TraceGweights(os);
    }

    void Save() override {
        ofstream param_os((name + ".param").c_str());
        param_os << GetModelType() << '\n';
        param_os << datafile << '\t' << treefile << '\n';
        param_os << aadistmodel << '\n';
        param_os << Gcat << '\t' << xmin << '\t' << xmax << '\n';
        param_os << 0 << '\n';
        param_os << every << '\t' << until << '\t' << size << '\n';
        model->ToStream(param_os);
    }
};

int main(int argc, char *argv[]) {
    string name = "";
    SelACChain *chain = 0;

    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        chain = new SelACChain(name);
    }

    // new chain
    else {
        string datafile = "";
        string treefile = "";
        int Gcat = 8;
        double xmin = -5;
        double xmax = 5;

        // uncons by default
        int aadistmodel = 1;
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
                } else if (s == "-gcat") {
                    i++;
                    Gcat = atoi(argv[i]);
                    i++;
                    xmin = atof(argv[i]);
                    i++;
                    xmax = atof(argv[i]);
                } else if (s == "-grantham")    {
                    aadistmodel = 0;
                } else if (s == "-uncons")  {
                    aadistmodel = 1;
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
            cerr << "selac -d <alignment> -t <tree> -gcat <gcat> <chainname> \n";
            cerr << '\n';
            exit(1);
        }

        chain = new SelACChain(datafile, treefile, aadistmodel, Gcat, xmin, xmax,
                                            every, until, name, force);
    }

    cerr << "chain " << name << " started\n";
    chain->Start();
    cerr << "chain " << name << " stopped\n";
    cerr << chain->GetSize()
         << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
    chain->GetModel()->Trace(cerr);
}
