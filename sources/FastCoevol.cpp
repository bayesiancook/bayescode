#include <cmath>
#include <fstream>
#include "Chain.hpp"
#include "FastCoevolModel.hpp"
using namespace std;

/**
 * \brief Chain object for running an MCMC under FastCoevolModel
 */

class FastCoevolChain : public Chain {
  private:
    // Chain parameters
    string modeltype;
    string contdatafile, treefile, rootfile;
    string dsomsuffstatfile;
    int wndsmode, wnommode;
    double beta, dbeta;

  public:
    FastCoevolChain(string incontdatafile, string intreefile, string inrootfile, string indsomsuffstatfile, int inwndsmode, int inwnommode, double inbeta, double indbeta, int inevery, int inuntil, string inname,
                     int force)
        : modeltype("FASTCOEVOLDNDS"), contdatafile(incontdatafile), treefile(intreefile), rootfile(inrootfile), dsomsuffstatfile(indsomsuffstatfile), wndsmode(inwndsmode), wnommode(inwnommode), beta(inbeta), dbeta(indbeta) {
        every = inevery;
        until = inuntil;
        name = inname;
        New(force);
    }

    //! constructor for re-opening an already existing chain from file -- calls
    //! Open
    FastCoevolChain(string filename) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new FastCoevolModel(contdatafile, treefile, rootfile, dsomsuffstatfile, wndsmode, wnommode);
        GetModel()->Allocate();
        GetModel()->Update();
        cerr << "-- Reset" << endl;
        Reset(force);
        GetModel()->SetBeta(beta);
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
        is >> contdatafile >> treefile >> rootfile;
        is >> dsomsuffstatfile;
        is >> wndsmode >> wnommode;
        int tmp;
        is >> tmp;
        if (tmp) {
            is >> beta >> dbeta;
            is >> tmp;
            if (tmp)    {
                cerr << "-- Error when reading model\n";
                exit(1);
            }
        }
        is >> every >> until >> size;

        if (modeltype == "FASTCOEVOLDNDS") {
            model = new FastCoevolModel(contdatafile, treefile, rootfile, dsomsuffstatfile, wndsmode, wnommode);
        } else {
            cerr << "-- Error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }
        GetModel()->Allocate();
        model->FromStream(is);
        GetModel()->SetBeta(beta);
        model->Update();
        cerr << size << " points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
        model->Trace(cerr);
    }

    void Save() override {
        ofstream param_os((name + ".param").c_str());
        param_os << GetModelType() << '\n';
        param_os << contdatafile << '\t' << treefile << '\t' << rootfile << '\n';
        param_os << dsomsuffstatfile << '\n';
        param_os << wndsmode << '\t' << wnommode << '\n';
        param_os << 1 << '\n';
        param_os << beta << '\t' << dbeta << '\n';
        param_os << 0 << '\n';
        param_os << every << '\t' << until << '\t' << size << '\n';
        model->ToStream(param_os);
    }

    void Move() override    {
        for (int i = 0; i < every; i++) {
            model->Move();
        }
        if (beta < 1.0) {
            beta += dbeta;
            if (beta > 1.0) {
                beta = 1.0;
            }
            GetModel()->SetBeta(beta);
        }
        SavePoint();
        Save();
        Monitor();
    }

    //! return the model, with its derived type (unlike ProbModel::GetModel)
    FastCoevolModel *GetModel() { return static_cast<FastCoevolModel *>(model); }

    //! return model type
    string GetModelType() override { return modeltype; }
};

int main(int argc, char *argv[]) {
    string name = "";
    FastCoevolChain *chain = 0;

    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        chain = new FastCoevolChain(name);
    }

    // new chain
    else {
        string contdatafile = "None";
        string treefile = "";
        string rootfile = "";
        string dsomsuffstatfile = "None";
        name = "";
        int force = 1;
        int every = 1;
        int until = -1;
        int wndsmode = 0;
        int wnommode = 0;

        double beta = 1.0;
        double dbeta = 0;

        try {
            if (argc == 1) {
                throw(0);
            }

            int i = 1;
            while (i < argc) {
                string s = argv[i];

                if (s == "-ss")  {
                    i++;
                    dsomsuffstatfile = argv[i];
                } else if (s == "-c")   {
                    i++;
                    contdatafile = argv[i];
                } else if ((s == "-t") || (s == "-T")) {
                    i++;
                    treefile = argv[i];
                } else if (s == "-r")   {
                    i++;
                    rootfile = argv[i];
                } else if (s == "-ugam")  {
                    wndsmode = wnommode = 1;
                } else if (s == "-wn")    {
                    wndsmode = wnommode = 2;
                } else if (s == "-vlin")  {
                    wndsmode = wnommode = 3;
                } else if (s == "-v3")  {
                    wndsmode = wnommode = 4;
                } else if (s == "-annealing")   {
                    i++;
                    beta = atof(argv[i]);
                    i++;
                    dbeta = atof(argv[i]);
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
            if ((dsomsuffstatfile == "") || (treefile == "") || (rootfile == "") || (name == "")) {
                throw(0);
            }
        } catch (...) {
            cerr << "fastcoevol -ss <suffstat> -c <contdata> -t <tree> -r <rootprior> <chainname> \n";
            cerr << '\n';
            exit(1);
        }

        chain = new FastCoevolChain(contdatafile, treefile, rootfile, dsomsuffstatfile, wndsmode, wnommode, beta, dbeta, every, until, name, force);
    }

    cerr << "chain " << name << " started\n";
    chain->Start();
    cerr << "chain " << name << " stopped\n";
    cerr << chain->GetSize()
         << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
    chain->GetModel()->Trace(cerr);
}
