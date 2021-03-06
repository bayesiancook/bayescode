#include <cmath>
#include <fstream>
#include "Chrono.hpp"
#include "MultiGeneDiscSelACModel.hpp"
#include "MultiGeneChain.hpp"

using namespace std;

MPI_Datatype Propagate_arg;

class MultiGeneDiscSelACChain : public MultiGeneChain {
  private:
    // Chain parameters
    string modeltype, datafile, treefile, initfile;
    int writegenedata;
    int Gcat;
    double xmin, xmax;
    int aadistmodel;
    int blmode, nucmode, aadistmode, gvarmode, psimode;

  public:
    MultiGeneDiscSelACModel *GetModel() {
        return static_cast<MultiGeneDiscSelACModel *>(model);
    }

    string GetModelType() override { return modeltype; }

    MultiGeneDiscSelACChain(string indatafile, string intreefile, string ininitfile, int inGcat, double inxmin, double inxmax, int inaadistmodel,
                                     int inblmode, int innucmode, int inaadistmode, int ingvarmode, int inpsimode, 
                                     int inevery, int inuntil,
                                     int inwritegenedata, string inname, int force, int inmyid,
                                     int innprocs)
        : MultiGeneChain(inmyid, innprocs),
          modeltype("MULTIGENEDISCSELAC"),
          datafile(indatafile),
          treefile(intreefile),
          initfile(ininitfile),
          Gcat(inGcat),
          xmin(inxmin), xmax(inxmax),
          aadistmodel(inaadistmodel),
          blmode(inblmode),
          nucmode(innucmode),
          aadistmode(inaadistmode),
          gvarmode(ingvarmode),
          psimode(inpsimode)    {
        every = inevery;
        until = inuntil;
        writegenedata = inwritegenedata;
        name = inname;
        New(force);
    }

    MultiGeneDiscSelACChain(string filename, int inmyid, int innprocs)
        : MultiGeneChain(inmyid, innprocs) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new MultiGeneDiscSelACModel(datafile, treefile, initfile, Gcat, xmin, xmax, aadistmodel,
                                                blmode, nucmode, aadistmode, gvarmode, psimode, 
                                                     myid, nprocs);
        if (!myid) {
            cerr << " -- allocate\n";
        }
        GetModel()->Allocate();
        if (!myid) {
            cerr << " -- update\n";
        }
        GetModel()->Update();
        Reset(force);
        if (!myid) {
            model->Trace(cerr);
        }
    }

    void Open() override {
        ifstream is((name + ".param").c_str());
        if (!is) {
            cerr << "-- Error : cannot find file : " << name << ".param\n";
            exit(1);
        }
        is >> modeltype;
        is >> datafile >> treefile >> initfile;
        is >> writegenedata;
        is >> Gcat >> xmin >> xmax >> aadistmodel;
        is >> blmode >> nucmode >> aadistmode >> gvarmode >> psimode;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "MULTIGENEAAMUTSELDSBDPOMEGA") {
            model = new MultiGeneDiscSelACModel(datafile, treefile, initfile, Gcat, xmin, xmax, aadistmodel, blmode,
                                                         nucmode, aadistmode, gvarmode, psimode, 
                                                         myid, nprocs);
        } else {
            cerr << "-- Error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }

        GetModel()->Allocate();
        model->FromStream(is);
        GetModel()->Update();

        if (!myid) {
            cerr << size << " points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
            model->Trace(cerr);
        }
    }

    void Save() override {
        if (!myid) {
            ofstream param_os((name + ".param").c_str());
            param_os << GetModelType() << '\n';
            param_os << datafile << '\t' << treefile << '\t' << initfile << '\n';
            param_os << writegenedata << '\n';
            param_os << Gcat << '\t' << xmin << '\t' << xmax << aadistmodel << '\n';
            param_os << blmode << '\t' << nucmode << '\t' << aadistmode << '\t';
            param_os << gvarmode << '\t' << psimode << '\n';
            param_os << 0 << '\n';
            param_os << every << '\t' << until << '\t' << size << '\n';
            GetModel()->MasterToStream(param_os);
        } else {
            GetModel()->SlaveToStream();
        }
    }

    void MakeFiles(int force) override {
        if (myid) {
            cerr << "error: in specialized makefiles\n";
            exit(1);
        }
        MultiGeneChain::MakeFiles(force);
        ofstream os((name + ".geneom").c_str());
        ofstream psios((name + ".genepsi").c_str());
        ofstream gos((name + ".genegvar").c_str());
        ofstream aaos((name + ".aadist").c_str());
        GetModel()->TraceAADistHeader(aaos);
    }

    void SavePoint() override {
        MultiGeneChain::SavePoint();
        if (! myid) {
            ofstream os((name + ".aadist").c_str(), ios_base::app);
            GetModel()->TraceAADist(os);
            if (writegenedata) {
                ofstream psios((name + ".genepsi").c_str(), ios_base::app);
                GetModel()->TracePsi(psios);
                ofstream gos((name + ".genegvar").c_str(), ios_base::app);
                GetModel()->TraceGvar(gos);
                ofstream os((name + ".geneom").c_str(), ios_base::app);
                GetModel()->TracePredictedDNDS(os);
            }
        }
    }
};

int main(int argc, char *argv[]) {
    Chrono chrono;
    chrono.Start();

    int myid = 0;
    int nprocs = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    string name = "";
    MultiGeneDiscSelACChain *chain = 0;

    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        chain = new MultiGeneDiscSelACChain(name, myid, nprocs);
    }

    // new chain
    else {
        string datafile = "";
        string treefile = "";
        string initfile = "None";
        int Gcat = 8;
	double xmin = -5.0;
	double xmax = 5.0;
        int force = 1;
        int every = 1;
        int until = -1;

        int blmode = 1;
        int nucmode = 2;
        int aadistmode = 2;
        int gvarmode = 1;
        int psimode = 1;
        int aadistmodel = 1;

        int writegenedata = 1;

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
                } else if (s == "-i")   {
                    i++;
                    initfile = argv[i];
                } else if (s == "-g") {
                    writegenedata = 0;
                } else if (s == "+g") {
                    writegenedata = 1;
                } else if (s == "-bl") {
                    i++;
                    string tmp = argv[i];
                    if (tmp == "shared") {
                        blmode = 2;
                    } else if (tmp == "shrunken") {
                        blmode = 1;
                    } else if ((tmp == "ind") || (tmp == "independent")) {
                        blmode = 0;
                    } else {
                        cerr << "error: does not recongnize command after -bl\n";
                        exit(1);
                    }
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
                } else if (s == "-fixaadist")   {
                    aadistmode = 3;
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
            cerr << "multigenediscselac -d <list> -t <tree> -gcat <gcat> <xmin> <xmax>"
                    "<chainname> \n";
            cerr << '\n';
            exit(1);
        }

        chain = new MultiGeneDiscSelACChain(
            datafile, treefile, initfile, Gcat, xmin, xmax, aadistmodel, blmode, nucmode, aadistmode, gvarmode, psimode,
            every, until, writegenedata, name, force, myid, nprocs);
    }

    chrono.Stop();
    if (!myid) {
        cout << "total time to set things up: " << chrono.GetTime() << '\n';
    }
    chrono.Reset();
    chrono.Start();
    if (!myid) {
        cerr << "chain " << name << " started\n";
    }
    chain->Start();
    if (!myid) {
        cerr << "chain " << name << " stopped\n";
        cerr << chain->GetSize()
             << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
        chain->GetModel()->Trace(cerr);
    }
    chrono.Stop();
    if (!myid) {
        cout << "total time in MCMC: " << chrono.GetTime() << '\n';
        cout << "total time in master moves: " << chain->GetModel()->GetMasterMoveTime() << '\n';
        cout << "mean total time in slave moves: " << chain->GetModel()->GetSlaveMoveTime() << '\n';
        cout << "mean total time in substitution mapping: " << chain->GetModel()->GetSlaveMapTime()
             << '\n';
    }

    MPI_Finalize();
}
