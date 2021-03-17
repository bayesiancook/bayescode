#include <cmath>
#include <fstream>
#include "Chain.hpp"
#include "CoevolModel.hpp"
using namespace std;

/**
 * \brief Chain object for running an MCMC under CoevolModel
 */

class CoevolChain : public Chain {
  private:
    // Chain parameters
    string modeltype;
    string datafile, contdatafile, treefile, rootfile;
    string suffstatfile, dsomsuffstatfile;
    GeneticCodeType codetype;

  public:
    CoevolChain(string indatafile, string incontdatafile, string intreefile, string inrootfile, string insuffstatfile, string indsomsuffstatfile, GeneticCodeType incodetype, int inevery, int inuntil, string inname,
                     int force)
        : modeltype("COEVOLDNDS"), datafile(indatafile), contdatafile(incontdatafile), treefile(intreefile), rootfile(inrootfile), suffstatfile(insuffstatfile), dsomsuffstatfile(indsomsuffstatfile), codetype(incodetype) {
        every = inevery;
        until = inuntil;
        name = inname;
        New(force);
    }

    //! constructor for re-opening an already existing chain from file -- calls
    //! Open
    CoevolChain(string filename) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new CoevolModel(datafile, contdatafile, treefile, rootfile, suffstatfile, dsomsuffstatfile, codetype);
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
        is >> datafile >> contdatafile >> treefile >> rootfile;
        is >> suffstatfile >> dsomsuffstatfile;
        is >> codetype;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "COEVOLDNDS") {
            model = new CoevolModel(datafile, contdatafile, treefile, rootfile, suffstatfile, dsomsuffstatfile, codetype);
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
        param_os << datafile << '\t' << contdatafile << '\t' << treefile << '\t' << rootfile << '\n';
        param_os << suffstatfile << '\t' << dsomsuffstatfile << '\n';
        param_os << codetype << '\n';
        param_os << 0 << '\n';
        param_os << every << '\t' << until << '\t' << size << '\n';
        model->ToStream(param_os);
    }

    //! return the model, with its derived type (unlike ProbModel::GetModel)
    CoevolModel *GetModel() { return static_cast<CoevolModel *>(model); }

    void SavePoint() override {
        Chain::SavePoint();
        ofstream pos((name + ".branchdsomsuffstat").c_str(), ios_base::app);
        GetModel()->TracedSOmegaPathSuffStat(pos);
        ofstream sos((name + ".branchsuffstat").c_str(), ios_base::app);
        GetModel()->TraceRelativePathSuffStat(sos);
    }

    void MakeFiles(int force) override {
        Chain::MakeFiles(force);
        // ofstream pos((name + ".branchomega").c_str());
        ofstream pos((name + ".branchdsomsuffstat").c_str());
        ofstream sos((name + ".branchsuffstat").c_str());
    }

    //! return model type
    string GetModelType() override { return modeltype; }
};

int main(int argc, char *argv[]) {
    string name = "";
    CoevolChain *chain = 0;

    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        chain = new CoevolChain(name);
    }

    // new chain
    else {
        string datafile = "";
        string contdatafile = "None";
        string treefile = "";
        string rootfile = "";
        string suffstatfile = "None";
        string dsomsuffstatfile = "None";
        GeneticCodeType codetype = Universal;
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
                } else if (s == "-dsomss")  {
                    i++;
                    dsomsuffstatfile = argv[i];
                } else if (s == "-ss")  {
                    i++;
                    suffstatfile = argv[i];
                } else if (s == "-c")   {
                    i++;
                    contdatafile = argv[i];
                } else if ((s == "-t") || (s == "-T")) {
                    i++;
                    treefile = argv[i];
                } else if (s == "-r")   {
                    i++;
                    rootfile = argv[i];
                } else if (s == "-f") {
                    force = 1;
                } else if ((s == "-mtvert") || (s == "-mtmam")) {
                    codetype = MtMam;
                } else if (s == "-universal")   {
                    codetype = Universal;
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
            if ((datafile == "") || (treefile == "") || (rootfile == "") || (name == "")) {
                throw(0);
            }
        } catch (...) {
            cerr << "coevol -d <alignment> -c <contdata> -t <tree> <chainname> \n";
            cerr << '\n';
            exit(1);
        }

        chain = new CoevolChain(datafile, contdatafile, treefile, rootfile, suffstatfile, dsomsuffstatfile, codetype, every, until, name, force);
    }

    cerr << "chain " << name << " started\n";
    chain->Start();
    cerr << "chain " << name << " stopped\n";
    cerr << chain->GetSize()
         << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
    chain->GetModel()->Trace(cerr);
}
