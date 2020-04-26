#include <cmath>
#include <fstream>
#include "Chrono.hpp"
#include "MultiGeneAAMutSelSparseOmegaModel.hpp"
#include "MultiGeneChain.hpp"

using namespace std;

MPI_Datatype Propagate_arg;

class MultiGeneAAMutSelSparseOmegaChain : public MultiGeneChain {
  private:
    // Chain parameters
    string modeltype, datafile, treefile;
    int writegenedata;
    int blmode, nucmode, omegamode, omegaprior, modalprior;
    double dposompihypermean, dposompihyperinvconc;
    double maxdposom;
    double epsilonhypermean, epsilonhyperinvconc;
    double pihypermean, pihyperinvconc;

  public:
    MultiGeneAAMutSelSparseOmegaModel *GetModel() {
        return static_cast<MultiGeneAAMutSelSparseOmegaModel *>(model);
    }

    string GetModelType() override { return modeltype; }

    MultiGeneAAMutSelSparseOmegaChain(string indatafile, string intreefile, 
                                     int inblmode, int innucmode,
                                     int inomegamode, int inomegaprior, int inmodalprior,
                                     double indposompihypermean, double indposompihyperinvconc, double inmaxdposom, 
                                     double inepsilonhypermean, double inepsilonhyperinvconc,
                                     double inpihypermean, double inpihyperinvconc,
                                     int inevery, int inuntil,
                                     int inwritegenedata, string inname, int force, int inmyid,
                                     int innprocs)
        : MultiGeneChain(inmyid, innprocs),
          modeltype("MULTIGENEAAMUTSELSparseOMEGA"),
          datafile(indatafile),
          treefile(intreefile),
          blmode(inblmode),
          nucmode(innucmode),
          omegamode(inomegamode),
          omegaprior(inomegaprior),
          modalprior(inmodalprior),
          dposompihypermean(indposompihypermean),
          dposompihyperinvconc(indposompihyperinvconc),
          maxdposom(inmaxdposom),
          epsilonhypermean(inepsilonhypermean),
          epsilonhyperinvconc(inepsilonhyperinvconc),
          pihypermean(inpihypermean),
          pihyperinvconc(inpihyperinvconc)  {
        every = inevery;
        until = inuntil;
        writegenedata = inwritegenedata;
        name = inname;
        New(force);
    }

    MultiGeneAAMutSelSparseOmegaChain(string filename, int inmyid, int innprocs)
        : MultiGeneChain(inmyid, innprocs) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new MultiGeneAAMutSelSparseOmegaModel(datafile, treefile, 
                blmode, nucmode, omegamode, omegaprior, modalprior,
                dposompihypermean, dposompihyperinvconc, maxdposom, 
                epsilonhypermean, epsilonhyperinvconc, pihypermean, pihyperinvconc,
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
        is >> datafile >> treefile;
        is >> writegenedata;
        is >> blmode >> nucmode >> omegamode >> omegaprior >> modalprior;
        is >> dposompihypermean >> dposompihyperinvconc >> maxdposom;
        is >> epsilonhypermean >> epsilonhyperinvconc;
        is >> pihypermean >> pihyperinvconc;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "MULTIGENEAAMUTSELSparseOMEGA") {
            model = new MultiGeneAAMutSelSparseOmegaModel(datafile, treefile,
                    blmode, nucmode, omegamode, omegaprior, modalprior,
                    dposompihypermean, dposompihyperinvconc, maxdposom, 
                    epsilonhypermean, epsilonhyperinvconc, pihypermean, pihyperinvconc,
                    myid, nprocs);

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
            param_os << blmode << '\t' << nucmode << '\t' << omegamode << '\t'
                     << omegaprior << '\t' << modalprior << '\n';
            param_os << dposompihypermean << '\t' << dposompihyperinvconc << '\n';
            param_os << maxdposom << '\n';
            param_os << epsilonhypermean << '\t' << epsilonhyperinvconc << '\n';
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
        if (writegenedata)  {
            if (omegamode != 3) {
                ofstream os((name + ".geneom").c_str());
            }
            ofstream os((name + ".genednds").c_str());
            ofstream eos((name + ".geneeps").c_str());
            ofstream pos((name + ".genepi").c_str());
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
                ofstream eos((name + ".geneeps").c_str(), ios_base::app);
                GetModel()->TraceEpsilon(eos);
                ofstream pos((name + ".genepi").c_str(), ios_base::app);
                GetModel()->TracePi(pos);
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
    MultiGeneAAMutSelSparseOmegaChain *chain = 0;

    string datafile = "";
    string treefile = "";

    int force = 0;
    int every = 1;
    int until = -1;

    int blmode = 1;
    int nucmode = 1;
    int omegamode = 3;
    int omegaprior = 0;
    int modalprior = 1;

    double dposompihypermean = 0.1;
    double dposompihyperinvconc = 0.2;
    double maxdposom = 0;

    double epsilonhypermean = 0.01;
    double epsilonhyperinvconc = 1.0;
    double pihypermean = 0.1;
    double pihyperinvconc = 0.5;

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
            } else if ((s == "-eps") || (s == "-epsilon"))   {
                i++;
                epsilonhypermean = atof(argv[i]);
                i++;
                epsilonhyperinvconc = atof(argv[i]);
            } else if (s == "-pi")  {
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
            } else if (s == "-dposompi") {
                i++;
                dposompihypermean = atof(argv[i]);
                i++;
                dposompihyperinvconc = atof(argv[i]);
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
        chain = new MultiGeneAAMutSelSparseOmegaChain(name, myid, nprocs);
    }
    else    {
        // new chain
        chain = new MultiGeneAAMutSelSparseOmegaChain(
            datafile, treefile, blmode, nucmode, omegamode, omegaprior, modalprior,
            dposompihypermean, dposompihyperinvconc, maxdposom, 
            epsilonhypermean, epsilonhyperinvconc, pihypermean, pihyperinvconc,
            every, until, writegenedata, name, force, myid, nprocs);
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
