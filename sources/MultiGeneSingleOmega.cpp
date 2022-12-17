#include <cmath>
#include <fstream>
#include "MultiGeneChain.hpp"
#include "MultiGeneSingleOmegaModel.hpp"

using namespace std;

MPI_Datatype Propagate_arg;

/**
 * \brief Chain object for running an MCMC under MultiGeneSingleOmegaModel
 */

class MultiGeneSingleOmegaChain : public MultiGeneChain {
  public:
    // Chain parameters
    string modeltype, datafile, treefile;
    int blmode, nucmode, omegamode;
    double omegahypermean, omegahyperinvshape;

    MultiGeneSingleOmegaModel *GetModel() {
        return static_cast<MultiGeneSingleOmegaModel *>(model);
    }

    string GetModelType() override { return modeltype; }

    //! \brief constructor for a new MCMC
    //!
    //! \param indatafile: name of file contanining sequence alignment
    //! \param intreefile: name of file contaning tree
    //! \param inevery: thinning factor
    //! \param inuntil: maximum MCMC sample size
    //! \param name: base name for all files related to this MCMC run
    //! \param force: overwrite existing files with same name
    //! \param inmyid, int innprocs: process id and total number of MPI processes
    MultiGeneSingleOmegaChain(string indatafile, string intreefile, int inblmode, int innucmode, int inomegamode, 
                              double inomegahypermean, double inomegahyperinvshape,
                              int inevery, int inuntil,
                              string inname, int force, int inmyid, int innprocs)
        : MultiGeneChain(inmyid, innprocs),
          modeltype("MULTIGENESINGLEOMEGA"),
          datafile(indatafile),
          treefile(intreefile) {
        blmode = inblmode;
        nucmode = innucmode;
        omegamode = inomegamode;
        omegahypermean = inomegahypermean;
        omegahyperinvshape = inomegahyperinvshape;
        every = inevery;
        until = inuntil;
        name = inname;
        New(force);
    }

    //! \brief constructor for opening and restarting an already existing chain
    MultiGeneSingleOmegaChain(string filename, int inmyid, int innprocs)
        : MultiGeneChain(inmyid, innprocs) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new MultiGeneSingleOmegaModel(datafile, treefile, myid, nprocs);
        GetModel()->SetAcrossGenesModes(blmode,nucmode,omegamode);
        GetModel()->SetOmegaHyperParameters(omegahypermean,omegahyperinvshape);
        if (!myid) {
            cerr << "allocate\n";
        }
        GetModel()->Allocate();
        if (!myid) {
            cerr << "update\n";
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
            cerr << "error : cannot find file : " << name << ".param\n";
            exit(1);
        }
        is >> modeltype;
        is >> datafile >> treefile;
        is >> blmode >> nucmode >> omegamode;
        is >> omegahypermean >> omegahyperinvshape;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "MULTIGENESINGLEOMEGA") {
            model = new MultiGeneSingleOmegaModel(datafile, treefile, myid, nprocs);
            GetModel()->SetAcrossGenesModes(blmode,nucmode,omegamode);
            GetModel()->SetOmegaHyperParameters(omegahypermean,omegahyperinvshape);
        } else {
            cerr << "error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }
        if (!myid) {
            cerr << "allocate\n";
        }
        GetModel()->Allocate();
        model->FromStream(is);
        if (!myid) {
            cerr << "update\n";
        }
        model->Update();
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
            param_os << blmode << '\t' << nucmode << '\t' << omegamode << '\n';
            param_os << omegahypermean << '\t' << omegahyperinvshape << '\n';
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
        if (GetModel()->GetBLMode() != 2)    {
            ofstream los((name + ".geneds").c_str());
        }
    }

    void SavePoint() override {
        MultiGeneChain::SavePoint();
        if (!myid) {
            ofstream os((name + ".geneom").c_str(), ios_base::app);
            GetModel()->TraceOmega(os);
            if (GetModel()->GetBLMode() != 2)    {
                ofstream los((name + ".geneds").c_str(), ios_base::app);
                GetModel()->TraceGeneTreeLength(los);
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

    MultiGeneSingleOmegaChain *chain = 0;
    string name = "";

    double maxtime = 0;

    string datafile = "";
    string treefile = "";
    int force = 0;
    int every = 1;
    int until = -1;
    int blmode = 1;
    int nucmode = 1;
    int omegamode = 1;
    double omegahypermean = 1.0;
    double omegahyperinvshape = 1.0;

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
            } else if (s == "-omega") {
                omegamode = 0;
                i++;
                string tmp = argv[i];
                if (tmp != "uninf") {
                    omegahypermean = atof(argv[i]);
                    i++;
                    omegahyperinvshape = atof(argv[i]);
                }
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
        cerr << "globom -d <alignment> -t <tree> <chainname> \n";
        cerr << '\n';
        exit(1);
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
        chain = new MultiGeneSingleOmegaChain(name, myid, nprocs);
    }
    else    {
        // new chain
        chain = new MultiGeneSingleOmegaChain(datafile, treefile, blmode, nucmode,
                omegamode, omegahypermean, omegahyperinvshape,
                every, until, name, force, myid,nprocs);
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
        /*
        cout << "total time in master moves: " << chain->GetModel()->GetMasterMoveTime() << '\n';
        cout << "mean total time in slave moves: " << chain->GetModel()->GetSlaveMoveTime() << '\n';
        cout << "mean total time in substitution mapping: " << chain->GetModel()->GetSlaveMapTime()
             << '\n';
        */
    }

    MPI_Finalize();
}

