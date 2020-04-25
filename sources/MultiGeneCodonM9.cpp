
// this is the multi-gene version of CodonM9Model. See CodonM9Model.hpp for
// further information.

#include <cmath>
#include <fstream>
#include "Chrono.hpp"
#include "MultiGeneChain.hpp"
#include "MultiGeneCodonM9Model.hpp"

using namespace std;

MPI_Datatype Propagate_arg;

/**
 * \brief MultiGeneChain object for running an MCMC under MultiGeneCodonM9Model
 */

class MultiGeneCodonM9Chain : public MultiGeneChain {
  private:
    // Chain parameters
    string modeltype, datafile, treefile;
    int writegenedata;
    int blmode, blsamplemode, nucmode, omegamode;
    double pihypermean, pihyperinvconc;
    double purifmeanhypermean, purifmeanhyperinvconc;
    double purifinvconchypermean, purifinvconchyperinvshape;
    vector<double> purifweighthypercenter;
    double purifweighthyperinvconc;
    double posmeanhypermean, posmeanhyperinvshape;
    double posinvshapehypermean, posinvshapehyperinvshape;
    double poswhypermean, poswhyperinvconc;
    int modalprior;

  public:
    MultiGeneCodonM9Model *GetModel() { return static_cast<MultiGeneCodonM9Model *>(model); }

    string GetModelType() override { return modeltype; }

    //! \brief constructor for a new MCMC
    //!
    //! \param indatafile: name of file contanining list of sequence alignments (with total number of genes as header)
    //! \param intreefile: name of file contaning tree
    //! \param inevery: thinning factor
    //! \param inuntil: maximum MCMC sample size
    //! \param inwritegenedata: if 1, then trace gene- and condition-specific
    //! shift probabilities in separate files; if 2, then also trace site-specific
    //! shift probabilities \param name: base name for all files related to this
    //! MCMC run \param force: overwrite existing files with same name \param
    //! inmyid, int innprocs: process id and total number of MPI processes
    MultiGeneCodonM9Chain(string indatafile, string intreefile, 
            int inblmode, int inblsamplemode, int innucmode, int inomegamode,
            double inpihypermean, double inpihyperinvconc,
            double inpurifmeanhypermean, double inpurifmeanhyperinvconc,
            double inpurifinvconchypermean, double inpurifinvconchyperinvshape,
            vector<double> inpurifweighthypercenter, double inpurifweighthyperinvconc,
            double inposmeanhypermean, double inposmeanhyperinvshape,
            double inposinvshapehypermean, double inposinvshapehyperinvshape,
            double inposwhypermean, double inposwhyperinvconc,
            int inmodalprior,
            int inevery, int inuntil, int inwritegenedata,
            string inname, int force, int inmyid, int innprocs)
        : MultiGeneChain(inmyid, innprocs),
          modeltype("MULTIGENECODONM9"),
          datafile(indatafile),
          treefile(intreefile),
          blmode(inblmode),
          blsamplemode(inblsamplemode),
          nucmode(innucmode),
          omegamode(inomegamode),
            pihypermean(inpihypermean),
            pihyperinvconc(inpihyperinvconc),
            purifmeanhypermean(inpurifmeanhypermean),
            purifmeanhyperinvconc(inpurifmeanhyperinvconc),
            purifinvconchypermean(inpurifinvconchypermean),
            purifinvconchyperinvshape(inpurifinvconchyperinvshape),
            purifweighthypercenter(inpurifweighthypercenter),
            purifweighthyperinvconc(inpurifweighthyperinvconc),
            posmeanhypermean(inposmeanhypermean),
            posmeanhyperinvshape(inposmeanhyperinvshape),
            posinvshapehypermean(inposinvshapehypermean),
            posinvshapehyperinvshape(inposinvshapehyperinvshape),
            poswhypermean(inposwhypermean),
            poswhyperinvconc(inposwhyperinvconc),
            modalprior(inmodalprior)    {

        every = inevery;
        until = inuntil;
        writegenedata = inwritegenedata;
        name = inname;
        New(force);
    }

    //! \brief constructor for opening and restarting an already existing chain
    MultiGeneCodonM9Chain(string filename, int inmyid, int innprocs)
        : MultiGeneChain(inmyid, innprocs) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new MultiGeneCodonM9Model(datafile, treefile, pihypermean, pihyperinvconc, myid, nprocs);
        GetModel()->SetAcrossGenesModes(blmode, nucmode, omegamode);
        GetModel()->SetBLSamplingMode(blsamplemode);
        GetModel()->SetMixtureHyperParameters(
            pihypermean, 
            purifmeanhypermean, purifmeanhyperinvconc,
            purifinvconchypermean, purifinvconchyperinvshape,
            purifweighthypercenter, purifweighthyperinvconc,
            posmeanhypermean, posmeanhyperinvshape,
            posinvshapehypermean, posinvshapehyperinvshape,
            poswhypermean, poswhyperinvconc);

        GetModel()->SetModalMixturePrior(modalprior);

        if (!myid) {
            cerr << "allocate\n";
        }
        GetModel()->Allocate();
        GetModel()->SetChainSize(GetSize());
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
            cerr << "Error : cannot find file : " << name << ".param\n";
            exit(1);
        }


        is >> modeltype;
        is >> datafile >> treefile;
        is >> writegenedata;
        is >> blmode >> blsamplemode >> nucmode >> omegamode;
        is >> pihypermean >> pihyperinvconc;
        is >> purifmeanhypermean >> purifmeanhyperinvconc;
        is >> purifinvconchypermean >> purifinvconchyperinvshape;
        is >> purifweighthypercenter >> purifweighthyperinvconc;
        is >> posmeanhypermean >> posmeanhyperinvshape;
        is >> posinvshapehypermean >> posinvshapehyperinvshape;
        is >> poswhypermean >> poswhyperinvconc;
        is >> modalprior;

        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "MULTIGENECODONM9") {
            model = new MultiGeneCodonM9Model(datafile, treefile, pihypermean, pihyperinvconc, myid, nprocs);
            GetModel()->SetAcrossGenesModes(blmode, nucmode, omegamode);
            GetModel()->SetBLSamplingMode(blsamplemode);
            GetModel()->SetMixtureHyperParameters(
                pihypermean, 
                purifmeanhypermean, purifmeanhyperinvconc,
                purifinvconchypermean, purifinvconchyperinvshape,
                purifweighthypercenter, purifweighthyperinvconc,
                posmeanhypermean, posmeanhyperinvshape,
                posinvshapehypermean, posinvshapehyperinvshape,
                poswhypermean, poswhyperinvconc);
            GetModel()->SetModalMixturePrior(modalprior);
        } else {
            cerr << "Error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }

        if (!myid) {
            cerr << "allocate\n";
        }
        GetModel()->Allocate();
        GetModel()->SetChainSize(GetSize());

        if (!myid) {
            cerr << "read from file\n";
        }
        GetModel()->FromStream(is);

        if (!myid) {
            cerr << "update\n";
        }
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
            param_os << blmode << '\t' << blsamplemode << '\t' << nucmode << '\t' << omegamode << '\n';
            param_os << pihypermean << '\t' << pihyperinvconc << '\n';
            param_os << purifmeanhypermean << '\t' << purifmeanhyperinvconc << '\n';
            param_os << purifinvconchypermean << '\t' << purifinvconchyperinvshape << '\n';
            param_os << purifweighthypercenter << '\t' << purifweighthyperinvconc << '\n';
            param_os << posmeanhypermean << '\t' << posmeanhyperinvshape << '\n';
            param_os << posinvshapehypermean << '\t' << posinvshapehyperinvshape << '\n';
            param_os << poswhypermean << '\t' << poswhyperinvconc << '\n';
            param_os << modalprior << '\n';
            param_os << 0 << '\n';
            param_os << every << '\t' << until << '\t' << size << '\n';
            GetModel()->MasterToStream(param_os);
        } else {
            GetModel()->SlaveToStream();
        }
    }

    void MakeFiles(int force) override {
        MultiGeneChain::MakeFiles(force);

        if (writegenedata >= 1) {
            ofstream pos((name + ".posw").c_str());
        }
        if (writegenedata == 2) {
            ofstream siteos((name + ".siteom").c_str());
        }
    }

    void SavePoint() override {
        MultiGeneChain::SavePoint();
        if (writegenedata >= 1) {
            if (!myid) {
                ofstream posw_os((name + ".posw").c_str(), ios_base::app);
                GetModel()->TracePosWeight(posw_os);
            }
        }
        if (writegenedata == 2) {
            if (!myid) {
                ofstream pp_os((name + ".siteom").c_str(), ios_base::app);
                GetModel()->MasterTraceSiteOmega(pp_os);
            } else {
                GetModel()->SlaveTraceSiteOmega();
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

    MultiGeneCodonM9Chain *chain = 0;
    string name = "";

    if (argc == 1)  {
        if (! myid)	{
            cerr << '\n';
            cerr << "multigenecodonm9\n";
            cerr << '\n';
        }
        MPI_Finalize();
        exit(0);
    }

    string datafile = "";
    string treefile = "";

    double pihypermean = 0.1;
    double pihyperinvconc = 0.2;

    double purifmeanhypermean = 0.5;
    double purifmeanhyperinvconc = 0.5;
    double purifinvconchypermean = 1.0;
    double purifinvconchyperinvshape = 1.0;
    vector<double> purifweighthypercenter(3,1.0/3);
    double purifweighthyperinvconc = 1.0;

    double poswhypermean = 0.5;
    double poswhyperinvconc = 0.1;

    double posmeanhypermean = 1.0;
    double posmeanhyperinvshape = 1.0;
    double posinvshapehypermean = 1.0;
    double posinvshapehyperinvshape = 1.0;

    int modalprior = 1;

    int blmode = 1;
    int nucmode = 1;

    int blsamplemode = 0;
    int omegamode = 1;

    int writegenedata = 1;

    int force = 0;
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
            } else if (s == "-modalprior")  {
                modalprior = 1;
            } else if (s == "-unconsprior")    {
                modalprior = 0;
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
            } else if (s == "-blint")   {
                blsamplemode = 1;
            } else if (s == "-blnoint") {
                blsamplemode = 0;
            } else if (s == "-pi") {
                i++;
                pihypermean = atof(argv[i]);
                i++;
                pihyperinvconc = atof(argv[i]);
            } else if (s == "-g") {
                writegenedata = 0;
            } else if (s == "+g") {
                writegenedata = 1;
            } else if (s == "+G") {
                writegenedata = 2;
            } else if (s == "-f") {
                force = 1;
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
        if (! myid) {
            cerr << "error in command\n";
        }
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
        chain = new MultiGeneCodonM9Chain(name, myid, nprocs);
    }
    else    {
        // new chain
        chain = new MultiGeneCodonM9Chain(
            datafile, treefile, blmode, blsamplemode, nucmode, omegamode,
            pihypermean, pihyperinvconc, 
            purifmeanhypermean, purifmeanhyperinvconc, purifinvconchypermean, purifinvconchyperinvshape, 
            purifweighthypercenter, purifweighthyperinvconc, 
            posmeanhypermean, posmeanhyperinvshape, posinvshapehypermean, posinvshapehyperinvshape, 
            poswhypermean, poswhyperinvconc, 
            modalprior,
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
