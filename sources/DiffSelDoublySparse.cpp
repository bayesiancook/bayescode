#include <cmath>
#include <fstream>
#include "Chain.hpp"
#include "DiffSelDoublySparseModel.hpp"
using namespace std;

/**
 * \brief Chain object for running an MCMC under the DiffSelDoublySparseModel
 */

class DiffSelDoublySparseChain : public Chain {
  private:
    // Chain parameters
    string modeltype, datafile, treefile;
    int ncond, nlevel, codonmodel;
    double fitnessshape;
    int fitnesscentermode;
    // -1: estimated
    // 1 : no mask
    double epsilon;
    double pihypermean;
    double shiftprobmean;
    double shiftprobinvconc;
    int burnin;

  public:
    DiffSelDoublySparseModel *GetModel() { return static_cast<DiffSelDoublySparseModel *>(model); }

    string GetModelType() override { return modeltype; }

    DiffSelDoublySparseChain(string indata, string intree, int inncond, int innlevel,
                             int incodonmodel, double infitnessshape, int infitnesscentermode,
                             double inepsilon, double inpihypermean, double inshiftprobmean, double inshiftprobinvconc, 
                             int inburnin, int inevery, int inuntil,
                             int insaveall, string inname, int force)
        : modeltype("DIFFSELSPARSE"),
          datafile(indata),
          treefile(intree),
          ncond(inncond),
          nlevel(innlevel),
          codonmodel(incodonmodel),
          fitnessshape(infitnessshape),
          fitnesscentermode(infitnesscentermode),
          epsilon(inepsilon),
          pihypermean(inpihypermean), shiftprobmean(inshiftprobmean), shiftprobinvconc(inshiftprobinvconc) {
        burnin = inburnin;
        every = inevery;
        until = inuntil;
        saveall = insaveall;
        name = inname;
        New(force);
    }

    DiffSelDoublySparseChain(string filename) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new DiffSelDoublySparseModel(datafile, treefile, ncond, nlevel, codonmodel, epsilon, fitnessshape,
                                             pihypermean, shiftprobmean, shiftprobinvconc);
        if (burnin) {
            GetModel()->SetWithToggles(0);
        } else {
            GetModel()->SetWithToggles(1);
        }
        GetModel()->SetFitnessCenterMode(fitnesscentermode);
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
        is >> datafile >> treefile >> ncond >> nlevel;
        is >> codonmodel;
        is >> fitnessshape;
        is >> fitnesscentermode;
        is >> epsilon;
        is >> pihypermean >> shiftprobmean >> shiftprobinvconc;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> burnin;
        is >> every >> until >> saveall >> size;

        if (modeltype == "DIFFSELSPARSE") {
            model = new DiffSelDoublySparseModel(datafile, treefile, ncond, nlevel, codonmodel,
                                                 epsilon, fitnessshape, pihypermean, shiftprobmean, shiftprobinvconc);
        } else {
            cerr << "-- Error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }
        if (size < burnin) {
            GetModel()->SetWithToggles(0);
        } else {
            GetModel()->SetWithToggles(1);
        }
        GetModel()->SetFitnessCenterMode(fitnesscentermode);
        GetModel()->Allocate();
        GetModel()->FromStream(is);
        GetModel()->Update();
        cerr << size << " points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
        model->Trace(cerr);
    }

    void Save() override {
        if (size == burnin) {
            GetModel()->SetWithToggles(1);
        }

        ofstream param_os((name + ".param").c_str());
        param_os << GetModelType() << '\n';
        param_os << datafile << '\t' << treefile << '\t' << ncond << '\t' << nlevel << '\n';
        param_os << codonmodel << '\n';
        param_os << fitnessshape << '\n';
        param_os << fitnesscentermode << '\n';
        param_os << epsilon << '\n';
        param_os << pihypermean << '\t' << shiftprobmean << '\t' << shiftprobinvconc << '\n';
        param_os << 0 << '\n';
        param_os << burnin << '\n';
        param_os << every << '\t' << until << '\t' << saveall << '\t' << size << '\n';

        model->ToStream(param_os);
    }

    void SavePoint() override {
        Chain::SavePoint();
        for (int k = 0; k < ncond; k++) {
            ostringstream s;
            s << name << "_" << k;
            if (k) {
                ofstream tos((s.str() + ".shifttoggle").c_str(), ios_base::app);
                GetModel()->TraceToggle(k, tos);
            }
            ofstream fos((s.str() + ".fitness").c_str(), ios_base::app);
            GetModel()->TraceFitness(k, fos);
        }
    }

    void MakeFiles(int force) override {
        Chain::MakeFiles(force);
        for (int k = 0; k < ncond; k++) {
            ostringstream s;
            s << name << "_" << k;
            if (k) {
                ofstream tos((s.str() + ".shifttoggle").c_str());
            }
            ofstream fos((s.str() + ".fitness").c_str());
        }
    }
};

int main(int argc, char *argv[]) {
    string name = "";
    DiffSelDoublySparseChain *chain = 0;

    // this is an already existing chain on the disk; reopen and restart
    if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        chain = new DiffSelDoublySparseChain(name);
    }

    // this is a new chain
    else {
        string datafile = "";
        string treefile = "";
        int ncond = 2;
        int nlevel = 2;
        int codonmodel = 1;
        double epsilon = -1;

        double fitnessshape = 20;
        int fitnesscentermode = 3;

        double pihypermean = 0.1;
        double shiftprobmean = 0.1;
        double shiftprobinvconc = 0.1;

        name = "";
        int burnin = 0;
        int every = 1;
        int until = -1;
        int saveall = 1;
        int force = 0;

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
                } else if (s == "-shape") {
                    i++;
                    string tmp = argv[i];
                    if (s == "free") {
                        fitnessshape = 0;
                    } else {
                        fitnessshape = atof(argv[i]);
                    }
                } else if (s == "-center") {
                    i++;
                    string tmp = argv[i];
                    if (s == "free") {
                        fitnesscentermode = 0;
                    } else if ((s == "fixed") || (s == "uniform")) {
                        fitnesscentermode = 3;
                    }
                } else if ((s == "-eps") || (s == "-epsilon")) {
                    i++;
                    string tmp = argv[i];
                    if (tmp == "free") {
                        epsilon = -1;
                    } else {
                        epsilon = atof(argv[i]);
                    }
                } else if (s == "-pi")  {
                    i++;
                    pihypermean = atof(argv[i]);
                } else if (s == "-shiftprob")   {
                    i++;
                    shiftprobmean = atof(argv[i]);
                    i++;
                    shiftprobinvconc = atof(argv[i]);
                } else if (s == "-b") {
                    i++;
                    burnin = atoi(argv[i]);
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
        } catch (...) {
            cerr << "error in command\n";
            exit(1);
        }
        chain = new DiffSelDoublySparseChain(datafile, treefile, ncond, nlevel, codonmodel,
                                             fitnessshape, fitnesscentermode, epsilon,
                                             pihypermean, shiftprobmean, shiftprobinvconc,
                                             burnin, every, until, saveall, name, force);
    }
    cerr << "chain " << name << " started\n";
    chain->Start();
    cerr << "chain " << name << " stopped\n";
    cerr << chain->GetSize()
         << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
    chain->GetModel()->Trace(cerr);
}
