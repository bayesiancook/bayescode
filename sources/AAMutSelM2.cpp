#include <cmath>
#include <fstream>
using namespace std;
#include "AAMutSelM2Model.hpp"
#include "Chain.hpp"

/**
 * \brief Chain object for running an MCMC under AAMutSelM2Model
 */

class AAMutSelM2Chain : public Chain {
  private:
    // Chain parameters
    string modeltype, datafile, treefile;

    // prior probability for the gene to be under positive selection (i.e. prior
    // prob that posw > 0)
    double pi;

    // Beta prior for posw (assuming posw>0)
    double poswhypermean;
    double poswhyperinvconc;

    // Gamma prior for dposom = omega_pos - 1 (with hyper mean and inverse shape
    // parameter)
    double dposomhypermean;
    double dposomhyperinvshape;

    int Ncat, baseNcat;

  public:
    AAMutSelM2Model *GetModel() { return static_cast<AAMutSelM2Model *>(model); }

    string GetModelType() override { return modeltype; }

    AAMutSelM2Chain(string indatafile, string intreefile,
                            double inpi, double inposwhypermean, double inposwhyperinvconc,
                            double indposomhypermean, double indposomhyperinvshape,
                            int inNcat, int inbaseNcat, int inevery,
                            int inuntil, string inname, int force)
        : modeltype("AAMUTSELDSBDPOMEGA"),
          datafile(indatafile),
          treefile(intreefile),
          pi(inpi),
          poswhypermean(inposwhypermean),
          poswhyperinvconc(inposwhyperinvconc),
          dposomhypermean(indposomhypermean),
          dposomhyperinvshape(indposomhyperinvshape),
          Ncat(inNcat),
          baseNcat(inbaseNcat) {
        every = inevery;
        until = inuntil;
        name = inname;
        New(force);
    }

    AAMutSelM2Chain(string filename) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        cerr << "new model\n";
        model = new AAMutSelM2Model(datafile, treefile, Ncat, baseNcat);
        GetModel()->SetOmegaMixtureHyperParameters(pi, poswhypermean, poswhyperinvconc, dposomhypermean, dposomhyperinvshape);
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
        is >> pi >> poswhypermean >> poswhyperinvconc;
        is >> dposomhypermean >> dposomhyperinvshape;
        is >> Ncat >> baseNcat;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "AAMUTSELDSBDPOMEGA") {
            model = new AAMutSelM2Model(datafile, treefile, Ncat, baseNcat);
            GetModel()->SetOmegaMixtureHyperParameters(pi, poswhypermean, poswhyperinvconc, dposomhypermean, dposomhyperinvshape);
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
        param_os << pi << '\t' << poswhypermean << '\t' << poswhyperinvconc << '\t';
        param_os << dposomhypermean << '\t' << dposomhyperinvshape << '\t';
        param_os << Ncat << '\t' << baseNcat << '\n';
        param_os << 0 << '\n';
        param_os << every << '\t' << until << '\t' << size << '\n';
        model->ToStream(param_os);
    }
};

int main(int argc, char *argv[]) {
    string name = "";
    AAMutSelM2Chain *chain = 0;

    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        chain = new AAMutSelM2Chain(name);
    }

    // new chain
    else {
        string datafile = "";
        string treefile = "";
        int Ncat = 100;
        int baseNcat = 1;
        double pi = 0.1;
        double poswhypermean = 0.5;
        double poswhyperinvconc = 0.1;
        double dposomhypermean = 1.0;
        double dposomhyperinvshape = 0.5;
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
                } else if (s == "-ncat") {
                    i++;
                    Ncat = atoi(argv[i]);
                } else if (s == "-basencat") {
                    i++;
                    baseNcat = atoi(argv[i]);
                } else if (s == "-mixomega") {
                    i++;
                    pi = atof(argv[i]);
                    i++;
                    poswhypermean = atof(argv[i]);
                    i++;
                    poswhyperinvconc = atof(argv[i]);
                    i++;
                    dposomhypermean = atof(argv[i]);
                    i++;
                    dposomhyperinvshape = atof(argv[i]);
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
            cerr << "aamutseldp -d <alignment> -t <tree> -ncat <ncat> <chainname> \n";
            cerr << '\n';
            exit(1);
        }

        chain = new AAMutSelM2Chain(datafile, treefile,
                pi, poswhypermean, poswhyperinvconc,
                dposomhypermean, dposomhyperinvshape,
                Ncat, baseNcat,
                every, until, name, force);
    }

    cerr << "chain " << name << " started\n";
    chain->Start();
    cerr << "chain " << name << " stopped\n";
    cerr << chain->GetSize()
         << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
    chain->GetModel()->Trace(cerr);
}
