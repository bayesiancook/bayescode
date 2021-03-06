#include <cmath>
#include <fstream>
#include "Chrono.hpp"
#include "MultiGeneAAMutSelDSBDPOmegaModel.hpp"
#include "MultiGeneChain.hpp"

using namespace std;

MPI_Datatype Propagate_arg;

class MultiGeneAAMutSelDSBDPOmegaChain : public MultiGeneChain {
  private:
    // Chain parameters
    string modeltype, datafile, treefile;
    int writegenedata;
    int Ncat;
    int baseNcat;
    int blmode, nucmode, basemode, omegamode, omegaprior, modalprior;
    double pihypermean, pihyperinvconc;
    double maxdposom;

  public:
    MultiGeneAAMutSelDSBDPOmegaModel *GetModel() {
        return static_cast<MultiGeneAAMutSelDSBDPOmegaModel *>(model);
    }

    string GetModelType() override { return modeltype; }

    MultiGeneAAMutSelDSBDPOmegaChain(string indatafile, string intreefile, int inNcat,
                                     int inbaseNcat, int inblmode, int innucmode, int inbasemode,
                                     int inomegamode, int inomegaprior, int inmodalprior, double inpihypermean,
                                     double inpihyperinvconc, double inmaxdposom, int inevery, int inuntil,
                                     int inwritegenedata, string inname, int force, int inmyid,
                                     int innprocs)
        : MultiGeneChain(inmyid, innprocs),
          modeltype("MULTIGENEAAMUTSELDSBDPOMEGA"),
          datafile(indatafile),
          treefile(intreefile),
          Ncat(inNcat),
          baseNcat(inbaseNcat),
          blmode(inblmode),
          nucmode(innucmode),
          basemode(inbasemode),
          omegamode(inomegamode),
          omegaprior(inomegaprior),
          modalprior(inmodalprior),
          pihypermean(inpihypermean),
          pihyperinvconc(inpihyperinvconc),
          maxdposom(inmaxdposom)  {
        every = inevery;
        until = inuntil;
        writegenedata = inwritegenedata;
        name = inname;
        New(force);
    }

    MultiGeneAAMutSelDSBDPOmegaChain(string filename, int inmyid, int innprocs)
        : MultiGeneChain(inmyid, innprocs) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new MultiGeneAAMutSelDSBDPOmegaModel(datafile, treefile, Ncat, baseNcat, blmode,
                                                     nucmode, basemode, omegamode, omegaprior, modalprior,
                                                     pihypermean, pihyperinvconc, maxdposom, myid, nprocs);
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
        is >> datafile >> treefile;
        is >> writegenedata;
        is >> Ncat >> baseNcat;
        is >> blmode >> nucmode >> basemode >> omegamode >> omegaprior >> modalprior;
        is >> pihypermean >> pihyperinvconc;
        maxdposom = 0;
        int tmp;
        is >> tmp;
        if (tmp) {
            is >> maxdposom;
            is >> tmp;
            if (tmp)    {
                cerr << "-- Error when reading model\n";
                exit(1);
            }
        }
        is >> every >> until >> size;

        if (modeltype == "MULTIGENEAAMUTSELDSBDPOMEGA") {
            model = new MultiGeneAAMutSelDSBDPOmegaModel(datafile, treefile, Ncat, baseNcat, blmode,
                                                         nucmode, basemode, omegamode, omegaprior, modalprior,
                                                         pihypermean, pihyperinvconc, maxdposom, myid, nprocs);
        } else {
            cerr << "-- Error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }

        GetModel()->Allocate();
        GetModel()->SetChainSize(GetSize());
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
            param_os << datafile << '\t' << treefile << '\n';
            param_os << writegenedata << '\n';
            param_os << Ncat << '\t' << baseNcat << '\n';
            param_os << blmode << '\t' << nucmode << '\t' << basemode << '\t' << omegamode << '\t'
                     << omegaprior << '\t' << modalprior << '\n';
            param_os << pihypermean << '\t' << pihyperinvconc << '\n';
            param_os << 1 << '\n';
            param_os << maxdposom << '\n';
            param_os << 0 << '\n';
            param_os << every << '\t' << until << '\t' << size << '\n';
            GetModel()->MasterToStream(param_os);
        } else {
            GetModel()->SlaveToStream();
        }
    }

    /*
    void Monitor() override {
        if (myid)   {
            cerr << "error: in specialized monitor\n";
            exit(1);
        }
        MultiGeneChain::Monitor();
        // ofstream trace_os((name + ".basemix").c_str(), ios_base::app);
        ofstream trace_os((name + ".basemix").c_str());
        GetModel()->TraceMixture(trace_os);
        ofstream logo_os((name + ".basemixlogo").c_str());
        GetModel()->PrintBaseMixtureLogo(logo_os);
        ofstream sample_os((name + ".basesamplelogo").c_str());
        GetModel()->PrintBaseSampleLogo(sample_os);
    }
    */

    void MakeFiles(int force) override {
        if (myid) {
            cerr << "error: in specialized makefiles\n";
            exit(1);
        }
        MultiGeneChain::MakeFiles(force);
        if (writegenedata)  {
            if (omegamode != 3) {
                ofstream os((name + ".geneom").c_str());
            }
            ofstream os((name + ".genednds").c_str());
            if (writegenedata == 2) {
                ofstream os((name + ".sitednds").c_str());
            }
        }
    }

    void SavePoint() override {
        MultiGeneChain::SavePoint();
        if (writegenedata) {
            if (!myid) {
                if (omegamode != 3) {
                    ofstream os((name + ".geneom").c_str(), ios_base::app);
                    GetModel()->TraceOmega(os);
                }
                ofstream dos((name + ".genednds").c_str(), ios_base::app);
                GetModel()->TracePredictedDNDS(dos);
            }
        }
        if (writegenedata == 2) {
            if (!myid) {
                ofstream os((name + ".sitednds").c_str(), ios_base::app);
                GetModel()->MasterTraceSitePredictedDNDS(os);
            } else {
                GetModel()->SlaveTraceSitePredictedDNDS();
            }
        }
    }
};

int main(int argc, char *argv[]) {
    Chrono chrono;
    chrono.Start();

    int myid = 0;
    int nprocs = 0;

    double maxtime = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    string name = "";
    MultiGeneAAMutSelDSBDPOmegaChain *chain = 0;

    string datafile = "";
    string treefile = "";
    int Ncat = 100;
    int baseNcat = 1;
    int force = 0;
    int every = 1;
    int until = -1;

    int blmode = 1;
    int nucmode = 1;
    int basemode = 0;
    int omegamode = 3;
    int omegaprior = 0;
    int modalprior = 1;

    double pihypermean = 0.1;
    double pihyperinvconc = 0.2;

    double maxdposom = 0;

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
            } else if (s == "-g") {
                writegenedata = 0;
            } else if (s == "+g") {
                writegenedata = 1;
            } else if (s == "+G") {
                writegenedata = 2;
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
            } else if (s == "-ncat") {
                i++;
                Ncat = atoi(argv[i]);
            } else if (s == "-basencat") {
                i++;
                baseNcat = atoi(argv[i]);
            } else if (s == "-basemix") {
                i++;
                string tmp = argv[i];
                if (tmp == "shared") {
                    basemode = 2;
                } else if ((tmp == "ind") || (tmp == "independent")) {
                    basemode = 0;
                } else {
                    cerr << "error: does not recognize command after -basemix\n";
                    throw(0);
                }
            } else if (s == "-maxdposom")    {
                i++;
                maxdposom = atof(argv[i]);
            } else if (s == "-fixomega") {
                omegamode = 3;
            } else if (s == "-freeomega") {
                omegamode = 1;
            } else if (s == "-gamomega") {
                omegaprior = 0;
            } else if ((s == "-mixomega") || (s == "-gammixomega")) {
                omegaprior = 1;
                modalprior = 0;
            } else if (s == "-loggammixomega") {
                omegaprior = 2;
                modalprior = 0;
            } else if (s == "-cauchymixomega") {
                omegaprior = 3;
                modalprior = 0;
            } else if (s == "-modalprior")  {
                modalprior = 1;
            } else if (s == "-unconsprior") {
                modalprior = 0;
            } else if (s == "-maxtime") {
                i++;
                maxtime = atof(argv[i]);
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
        cerr << "multigeneaamutselddp -d <list> -t <tree> -ncat <ncat> "
                "<chainname> \n";
        cerr << '\n';
        MPI_Finalize();
        exit(0);
    }

    if ((datafile == "") && (treefile == ""))   {
        // existing chain
        // check whether chain is already existing (and running) before extending it
        if (! force)    {
            ifstream is((name + ".run").c_str());
            int tmp;
            is >> tmp;
            if (tmp)    {
                if (! myid) {
                    cerr << "error: chain still running\n";
                }
                MPI_Finalize();
                exit(0);
            }
        }
        chain = new MultiGeneAAMutSelDSBDPOmegaChain(name, myid, nprocs);
    }
    else    {
        // new chain
        chain = new MultiGeneAAMutSelDSBDPOmegaChain(
            datafile, treefile, Ncat, baseNcat, blmode, nucmode, basemode, omegamode, omegaprior, modalprior,
            pihypermean, pihyperinvconc, maxdposom, every, until, writegenedata, name, force, myid, nprocs);
    }

    chrono.Stop();
    if (!myid) {
        cout << "total time to set things up: " << chrono.GetTime() << '\n';
        if (maxtime > 0)    {
            maxtime -= chrono.GetTime() / 3600000;
            cout << "remaining time: " << maxtime << '\n';
            if (maxtime < 0)    {
                cerr << "error: maxtime already exceeded\n";
                MPI_Finalize();
                exit(1);
            }
        }
        chain->SetMaxTime(maxtime);
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
