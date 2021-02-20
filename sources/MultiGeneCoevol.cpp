#include <cmath>
#include <fstream>
#include "MultiGeneChain.hpp"
#include "MultiGeneCoevolModel.hpp"

using namespace std;

MPI_Datatype Propagate_arg;

/**
 * \brief Chain object for running an MCMC under MultiGeneCoevolModel
 */

class MultiGeneCoevolChain : public MultiGeneChain {
  private:
    // Chain parameters
    string modeltype, datafile, contdatafile, treefile, rootfile;
    GeneticCodeType codetype;
    int nucmode;

  public:
    MultiGeneCoevolModel *GetModel() {
        return static_cast<MultiGeneCoevolModel *>(model);
    }

    string GetModelType() override { return modeltype; }

    MultiGeneCoevolChain(string indatafile, string incontdatafile, string intreefile, string inrootfile, GeneticCodeType incodetype, int innucmode,
                              int inevery, int inuntil,
                              string inname, int force, int inmyid, int innprocs)
        : MultiGeneChain(inmyid, innprocs),
          modeltype("MULTIGENECOEVOLOMEGA"),
          datafile(indatafile),
	  contdatafile(incontdatafile),
          treefile(intreefile), rootfile(inrootfile), codetype(incodetype), nucmode(innucmode)	{
        every = inevery;
        until = inuntil;
        name = inname;
        New(force);
    }

    //! \brief constructor for opening and restarting an already existing chain
    MultiGeneCoevolChain(string filename, int inmyid, int innprocs)
        : MultiGeneChain(inmyid, innprocs) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new MultiGeneCoevolModel(datafile, contdatafile, treefile, rootfile, codetype, myid, nprocs);
	GetModel()->SetNucMode(nucmode);
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
        is >> datafile >> contdatafile >> treefile >> rootfile;
	is >> codetype;
        is >> nucmode;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "MULTIGENECOEVOLOMEGA") {
            model = new MultiGeneCoevolModel(datafile, contdatafile, treefile, rootfile, codetype, myid, nprocs);
	    GetModel()->SetNucMode(nucmode);
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
            param_os << datafile << '\t' << contdatafile << '\t' << treefile << '\t' << rootfile << '\n';
	    param_os << codetype << '\n';
            param_os << nucmode << '\n';
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
        // ofstream pos((name + ".branchomega").c_str());
    }

    void SavePoint() override {
        MultiGeneChain::SavePoint();
        if (!myid) {
            /*
            ofstream pos((name + ".branchomega").c_str(), ios_base::app);
            GetModel()->TraceOmega(pos);
            */
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

    MultiGeneCoevolChain *chain = 0;
    string name = "";

    double maxtime = 0;

    string datafile = "";
    string contdatafile = "None";
    string treefile = "";
    string rootfile = "";
    GeneticCodeType codetype = Universal;

    int force = 0;
    int every = 1;
    int until = -1;
    int nucmode = 2;

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
            } else if ((s == "-mtvert") || (s == "-mtmam")) {
                codetype = MtMam;
            } else if (s == "-universal")   {
                codetype = Universal;
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
        chain = new MultiGeneCoevolChain(name, myid, nprocs);
    }
    else    {
        // new chain
        chain = new MultiGeneCoevolChain(datafile, contdatafile, treefile, rootfile, codetype, nucmode,
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
    }

    MPI_Finalize();
}

