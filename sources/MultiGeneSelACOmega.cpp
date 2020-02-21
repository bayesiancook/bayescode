#include <cmath>
#include <fstream>
#include "Chrono.hpp"
#include "MultiGeneSelACOmegaModel.hpp"
#include "MultiGeneChain.hpp"

using namespace std;

MPI_Datatype Propagate_arg;

class MultiGeneSelACOmegaChain : public MultiGeneChain {
  private:
    // Chain parameters
    string modeltype, datafile, treefile, initfile;
    int writegenedata;
    int Gcat;
    double Ginvshape;
    int aadistmodel;
    int blmode, nucmode, aadistmode, omegamode, omegaprior, modalprior;
    double pihypermean, pihyperinvconc;

  public:
    MultiGeneSelACOmegaModel *GetModel() {
        return static_cast<MultiGeneSelACOmegaModel *>(model);
    }

    string GetModelType() override { return modeltype; }

    MultiGeneSelACOmegaChain(string indatafile, string intreefile, string ininitfile, int inGcat, double inGinvshape, int inaadistmodel,
                                     int inblmode, int innucmode, int inaadistmode,
                                     int inomegamode, int inomegaprior, int inmodalprior, double inpihypermean,
                                     double inpihyperinvconc, int inevery, int inuntil,
                                     int inwritegenedata, string inname, int force, int inmyid,
                                     int innprocs)
        : MultiGeneChain(inmyid, innprocs),
          modeltype("MULTIGENEAAMUTSELDSBDPOMEGA"),
          datafile(indatafile),
          treefile(intreefile),
          initfile(ininitfile),
          Gcat(inGcat),
          Ginvshape(inGinvshape),
          aadistmodel(inaadistmodel),
          blmode(inblmode),
          nucmode(innucmode),
          aadistmode(inaadistmode),
          omegamode(inomegamode),
          omegaprior(inomegaprior),
          modalprior(inmodalprior),
          pihypermean(inpihypermean),
          pihyperinvconc(inpihyperinvconc) {
        every = inevery;
        until = inuntil;
        writegenedata = inwritegenedata;
        name = inname;
        New(force);
    }

    MultiGeneSelACOmegaChain(string filename, int inmyid, int innprocs)
        : MultiGeneChain(inmyid, innprocs) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new MultiGeneSelACOmegaModel(datafile, treefile, initfile, Gcat, Ginvshape, aadistmodel, blmode,
                                                     nucmode, aadistmode, omegamode, omegaprior, modalprior,
                                                     pihypermean, pihyperinvconc, myid, nprocs);
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
        is >> Gcat >> Ginvshape >> aadistmodel;
        is >> blmode >> nucmode >> aadistmode >> omegamode >> omegaprior >> modalprior;
        is >> pihypermean >> pihyperinvconc;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "MULTIGENEAAMUTSELDSBDPOMEGA") {
            model = new MultiGeneSelACOmegaModel(datafile, treefile, initfile, Gcat, Ginvshape, aadistmodel, blmode,
                                                         nucmode, aadistmode, omegamode, omegaprior, modalprior,
                                                         pihypermean, pihyperinvconc, myid, nprocs);
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
            param_os << Gcat << '\t' << Ginvshape << '\t' << aadistmodel << '\n';
            param_os << blmode << '\t' << nucmode << '\t' << aadistmode << '\t' << omegamode << '\t'
                     << omegaprior << '\t' << modalprior << '\n';
            param_os << pihypermean << '\t' << pihyperinvconc << '\n';
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
                ofstream os((name + ".geneom").c_str(), ios_base::app);
                if (omegamode == 3) {
                    GetModel()->TracePredictedDNDS(os);
                }
                else    {
                    GetModel()->TraceOmega(os);
                }
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
    MultiGeneSelACOmegaChain *chain = 0;

    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        chain = new MultiGeneSelACOmegaChain(name, myid, nprocs);
    }

    // new chain
    else {
        string datafile = "";
        string treefile = "";
        string initfile = "None";
        int Gcat = 4;
        double Ginvshape = 0.2;
        int force = 1;
        int every = 1;
        int until = -1;

        int blmode = 1;
        int nucmode = 1;
        int aadistmode = 2;
        int aadistmodel = 1;
        int omegamode = 3;
        int omegaprior = 0;
        int modalprior = 1;

        double pihypermean = 0.1;
        double pihyperinvconc = 0.2;

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
                } else if (s == "-pi") {
                    i++;
                    pihypermean = atof(argv[i]);
                    i++;
                    pihyperinvconc = atof(argv[i]);
                } else if (s == "-nucrates") {
                    i++;
                    string tmp = argv[i];
                    if (tmp == "shared") {
                        nucmode = 2;
                    } else if (tmp == "shrunken") {
                        nucmode = 1;
                    } else if ((tmp == "ind") || (tmp == "independent")) {
                        nucmode = 0;
                    } else {
                        cerr << "error: does not recongnize command after -nucrates\n";
                        exit(1);
                    }
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
                } else if (s == "-initg")   {
                    i++;
                    Ginvshape = atof(argv[i]);
                } else if (s == "-grantham")    {
                    aadistmodel = 0;
                } else if (s == "-uncons")  {
                    aadistmodel = 1;
                } else if (s == "-fixaadist")   {
                    aadistmode = 3;
                } else if (s == "-fixomega") {
                    omegamode = 3;
                } else if (s == "-freeomega") {
                    omegamode = 1;
                } else if (s == "-gamomega") {
                    omegaprior = 0;
                } else if (s == "-mixomega") {
                    omegaprior = 1;
                } else if (s == "-modalprior")  {
                    modalprior = 1;
                } else if (s == "-unconsprior") {
                    modalprior = 0;
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
            cerr << "multigeneselac -d <list> -t <tree> -gcat <gcat> "
                    "<chainname> \n";
            cerr << '\n';
            exit(1);
        }

        chain = new MultiGeneSelACOmegaChain(
            datafile, treefile, initfile, Gcat, Ginvshape, aadistmodel, blmode, nucmode, aadistmode, omegamode, omegaprior, modalprior,
            pihypermean, pihyperinvconc, every, until, writegenedata, name, force, myid, nprocs);
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
