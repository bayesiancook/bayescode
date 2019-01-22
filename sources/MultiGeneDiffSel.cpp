#include <cmath>
#include <fstream>
#include "Chrono.hpp"
#include "MultiGeneChain.hpp"
#include "MultiGeneDiffSelModel.hpp"

using namespace std;

MPI_Datatype Propagate_arg;

/**
 * \brief A MultiGeneChain object for running an MCMC under
 * MultiGeneDiffSelModel
 */
class MultiGeneDiffSelChain : public MultiGeneChain {
  private:
    // Chain parameters
    string modeltype, datafile, treefile;
    int ncond;
    int nlevel;
    int codonmodel;
    int blmode, nucmode;
    int writegenedata;

  public:
    MultiGeneDiffSelModel *GetModel() {
        return static_cast<MultiGeneDiffSelModel *>(model);
    }

    string GetModelType() override { return modeltype; }

    //! \brief constructor for a new MCMC
    //!
    //! \param indatafile: name of file contanining sequence alignment
    //! \param intreefile: name of file contaning tree (with branch names giving
    //! the allocation of branches to the conditions) \param inncond: number of
    //! conditions \param innlevel: number of levels of the model \param
    //! incodonmodel: type of codon substitution process (1: mutation-selection,
    //! 0: square-root) \param inevery: thinning factor \param inuntil: maximum
    //! MCMC sample size \param insaveall: if 1, then, save all information about
    //! each configuration visited during MCMC (into .chain file) \param
    //! inwritegenedata: if 1, then trace gene- and condition-specific shift
    //! probabilities in separate files; if 2, then also trace site-specific shift
    //! probabilities \param name: base name for all files related to this MCMC
    //! run \param force: overwrite existing files with same name \param inmyid,
    //! int innprocs: process id and total number of MPI processes
    MultiGeneDiffSelChain(string indatafile, string intreefile, int inncond, int innlevel,
                                int incodonmodel, int inblmode, int innucmode, int inevery, int inuntil, int insaveall,
                                int inwritegenedata, string inname, int force, int inmyid,
                                int innprocs)
        : MultiGeneChain(inmyid, innprocs),
          modeltype("MULTIGENEDIFFSEL"),
          datafile(indatafile),
          treefile(intreefile),
          ncond(inncond),
          nlevel(innlevel),
          codonmodel(incodonmodel), blmode(inblmode), nucmode(innucmode) {
        every = inevery;
        until = inuntil;
        saveall = insaveall;
        writegenedata = inwritegenedata;
        name = inname;
        New(force);
    }

    //! \brief constructor for opening and restarting an already existing chain
    MultiGeneDiffSelChain(string filename, int inmyid, int innprocs)
        : MultiGeneChain(inmyid, innprocs) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new MultiGeneDiffSelModel(datafile, treefile, ncond, nlevel, codonmodel, blmode, nucmode, myid,
                                                nprocs);
        if (!myid) {
            cerr << " -- master allocate\n";
        }
        GetModel()->Allocate();
        if (!myid) {
            cerr << " -- master unfold\n";
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
        is >> ncond >> nlevel >> codonmodel;
        is >> blmode >> nucmode;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> saveall >> writegenedata >> size;

        if (modeltype == "MULTIGENEDIFFSEL") {
            model = new MultiGeneDiffSelModel(datafile, treefile, ncond, nlevel, codonmodel, blmode, nucmode,
                                                    myid, nprocs);
        } else {
            cerr << "-- Error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }
        GetModel()->Allocate();
        GetModel()->FromStream(is);
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
            param_os << ncond << '\t' << nlevel << '\t' << codonmodel << '\n';
            param_os << blmode << '\t' << nucmode << '\n';
            param_os << 0 << '\n';
            param_os << every << '\t' << until << '\t' << saveall << '\t' << writegenedata << '\t'
                     << size << '\n';
            GetModel()->MasterToStream(param_os);
        } else {
            GetModel()->SlaveToStream();
        }
    }

    void MakeFiles(int force) override {
        MultiGeneChain::MakeFiles(force);
        if (writegenedata == 2) {
            ofstream os((name + ".fitness").c_str());
            for (int k = 1; k < ncond; k++) {
                ostringstream s;
                s << name << "_" << k;
                ofstream os((s.str() + ".delta").c_str());
            }
        }
    }

    void SavePoint() override {
        MultiGeneChain::SavePoint();
        if (writegenedata == 2) {
            if (!myid) {
                GetModel()->MasterTraceSiteStats(name, writegenedata);
            } else {
                GetModel()->SlaveTraceSiteStats(writegenedata);
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
    MultiGeneDiffSelChain *chain = 0;

    if (argc == 1)  {
        if (! myid)	{
            cerr << '\n';
            cerr << "The multi-gene version of the non-sparse differential selection model.\n";
            cerr << "see diffsel for the single-gene version.\n";
            cerr << '\n';
            cerr << "command: mpirun -np <n> multigenediffsel -d <alignment_list> -t <tree> -ncond <ncond> <chainname>\n";
            cerr << '\n';
            cerr << "chain options:\n";
            cerr << "\t-f: force overwrite of already existing chain\n";
            cerr << "\t-x <every> <until>: saving frequency and stopping time "
                    "(default: every = 1, until = -1)\n";
            cerr << "\t+G: with site-specific output files\n";
            cerr << '\n';
            cerr << "model options:\n";
            cerr << "\t-ncond <ncond>:  specify number of conditions\n";
            cerr << "\t-bl {shrunken|ind}: shrinkage mode for branch lengths\n";
            cerr << "\t-nucrates {shrunken|ind}: shrinkage mode for nucleotide substitution rates\n";
            cerr << '\n';
        }
        MPI_Finalize();
        exit(0);
    }

    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        chain = new MultiGeneDiffSelChain(name, myid, nprocs);
    }

    // new chain
    else {
        string datafile = "";
        string treefile = "";
        int ncond = 1;
        int nlevel = 1;
        int codonmodel = 1;
        int force = 1;
        int every = 1;
        int until = -1;
        int saveall = 1;
        int writegenedata = 0;
        int blmode = 1;
        int nucmode = 1;

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
                } else if (s == "-s") {
                    saveall = 0;
                } else if (s == "+s") {
                    saveall = 1;
                } else if (s == "-ncond") {
                    i++;
                    ncond = atoi(argv[i]);
                } else if (s == "-nlevel") {
                    i++;
                    nlevel = atoi(argv[i]);
                } else if (s == "-nucrates") {
                    i++;
                    string tmp = argv[i];
                    if (tmp == "shrunken") {
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
                    if (tmp == "shrunken") {
                        blmode = 1;
                    } else if ((tmp == "ind") || (tmp == "independent")) {
                        blmode = 0;
                    } else {
                        cerr << "error: does not recongnize command after -bl\n";
                        exit(1);
                    }
                } else if (s == "-g") {
                    writegenedata = 0;
                } else if (s == "+g") {
                    writegenedata = 1;
                } else if (s == "+G") {
                    writegenedata = 2;
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
            cerr << "error in command\n";
            cerr << '\n';
            exit(1);
        }

        chain = new MultiGeneDiffSelChain(datafile, treefile, ncond, nlevel, codonmodel, blmode, nucmode,
                                                every, until, saveall, writegenedata, name, force,
                                                myid, nprocs);
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
