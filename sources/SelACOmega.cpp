#include <cmath>
#include <fstream>
#include "SelACOmegaModel.hpp"
#include "Chain.hpp"
using namespace std;

/**
 * \brief Chain object for running an MCMC under SelACOmegaModel
 */

class SelACOmegaChain : public Chain {
  private:
    // Chain parameters
    string modeltype, datafile, treefile;
    int aadistmodel;
    int aadistmode;
    int omegamode, omegaprior;
    double dposompi, dposomhypermean, dposomhyperinvshape;
    int Gcat;

  public:
    SelACOmegaModel *GetModel() { return static_cast<SelACOmegaModel *>(model); }

    string GetModelType() override { return modeltype; }

    SelACOmegaChain(string indatafile, string intreefile, int inaadistmodel, int inaadistmode, int inomegamode, int inomegaprior,
                            double indposompi, double indposomhypermean,
                            double indposomhyperinvshape, int inGcat, int inevery,
                            int inuntil, string inname, int force)
        : modeltype("AAMUTSELDSBDPOMEGA"),
          datafile(indatafile),
          treefile(intreefile),
          aadistmodel(inaadistmodel),
          aadistmode(inaadistmode),
          omegamode(inomegamode),
          omegaprior(inomegaprior),
          dposompi(indposompi),
          dposomhypermean(indposomhypermean),
          dposomhyperinvshape(indposomhyperinvshape),
          Gcat(inGcat)  {
        every = inevery;
        until = inuntil;
        name = inname;
        New(force);
    }

    SelACOmegaChain(string filename) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        cerr << "new model\n";
        model =
            new SelACOmegaModel(datafile, treefile, aadistmodel, aadistmode, omegamode, omegaprior, Gcat);
        if (omegaprior == 1) {
            GetModel()->SetDPosOmHyperParameters(dposompi, dposomhypermean, dposomhyperinvshape);
        }
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
        is >> aadistmode;
        is >> omegamode >> omegaprior >> dposompi >> dposomhypermean >> dposomhyperinvshape;
        is >> Gcat;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "AAMUTSELDSBDPOMEGA") {
            model = new SelACOmegaModel(datafile, treefile, aadistmodel, aadistmode, omegamode, omegaprior, Gcat);
            if (omegaprior == 1) {
                GetModel()->SetDPosOmHyperParameters(dposompi, dposomhypermean,
                                                     dposomhyperinvshape);
            }
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
        param_os << aadistmodel << '\n';
        param_os << aadistmode << '\n';
        param_os << omegamode << '\t' << omegaprior << '\t' << dposompi << '\t' << dposomhypermean
                 << '\t' << dposomhyperinvshape << '\n';
        param_os << Gcat << '\n';
        param_os << 0 << '\n';
        param_os << every << '\t' << until << '\t' << size << '\n';
        model->ToStream(param_os);
    }
};

int main(int argc, char *argv[]) {
    string name = "";
    SelACOmegaChain *chain = 0;

    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        chain = new SelACOmegaChain(name);
    }

    // new chain
    else {
        string datafile = "";
        string treefile = "";
        int Gcat = 4;
        // uncons by default
        int aadistmode = 0;
        int aadistmodel = 1;
        int omegamode = 3;
        int omegaprior = 0;
        double dposompi = 0.1;
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
                } else if (s == "-gcat") {
                    i++;
                    Gcat = atoi(argv[i]);
                } else if (s == "-grantham")    {
                    aadistmodel = 0;
                } else if (s == "-uncons")  {
                    aadistmodel = 1;
                } else if (s == "-fixomega") {
                    omegamode = 3;
                } else if (s == "-freeomega") {
                    omegamode = 1;
                } else if (s == "-gamomega") {
                    omegaprior = 0;
                } else if (s == "-mixomega") {
                    omegaprior = 1;
                    i++;
                    dposompi = atof(argv[i]);
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
            cerr << "selac -d <alignment> -t <tree> -gcat <gcat> <chainname> \n";
            cerr << '\n';
            exit(1);
        }

        chain = new SelACOmegaChain(datafile, treefile, aadistmodel, aadistmode, omegamode, omegaprior, dposompi,
                                            dposomhypermean, dposomhyperinvshape, Gcat,
                                            every, until, name, force);
    }

    cerr << "chain " << name << " started\n";
    chain->Start();
    cerr << "chain " << name << " stopped\n";
    cerr << chain->GetSize()
         << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
    chain->GetModel()->Trace(cerr);
}
