
// this is the multi-gene version of CodonM2aModel. See CodonM2aModel.hpp for
// further information.

#include <cmath>
#include <fstream>
#include "Chrono.hpp"
#include "MultiGeneChain.hpp"
#include "MultiGeneCodonM2aModel.hpp"

using namespace std;

MPI_Datatype Propagate_arg;

/**
 * \brief MultiGeneChain object for running an MCMC under MultiGeneCodonM2aModel
 */

class MultiGeneCodonM2aChain : public MultiGeneChain {
  private:
    // Chain parameters
    string modeltype, datafile, treefile;
    int writegenedata;
    int blmode, blsamplemode, nucmode, purommode, dposommode, purwmode, poswmode;
    double pihypermean, pihyperinvconc;
    double puromhypermean, puromhyperinvconc;
    double dposomhypermean, dposomhyperinvshape;
    double purwhypermean, purwhyperinvconc;
    double poswhypermean, poswhyperinvconc;

  public:
    MultiGeneCodonM2aModel *GetModel() { return static_cast<MultiGeneCodonM2aModel *>(model); }

    string GetModelType() override { return modeltype; }

    //! \brief constructor for a new MCMC
    //!
    //! \param indatafile: name of file contanining sequence alignment
    //! \param intreefile: name of file contaning tree
    //! \param inevery: thinning factor
    //! \param inuntil: maximum MCMC sample size
    //! \param inwritegenedata: if 1, then trace gene- and condition-specific
    //! shift probabilities in separate files; if 2, then also trace site-specific
    //! shift probabilities \param name: base name for all files related to this
    //! MCMC run \param force: overwrite existing files with same name \param
    //! inmyid, int innprocs: process id and total number of MPI processes
    MultiGeneCodonM2aChain(string indatafile, string intreefile, int inblmode, int inblsamplemode, int innucmode,
                           int inpurommode, int indposommode, int inpurwmode, int inposwmode,
                           double inpihypermean, double inpihyperinvconc, double inpuromhypermean,
                           double inpuromhyperinvconc, double indposomhypermean,
                           double indposomhyperinvshape, double inpurwhypermean,
                           double inpurwhyperinvconc, double inposwhypermean,
                           double inposwhyperinvconc, int inevery, int inuntil, int inwritegenedata,
                           string inname, int force, int inmyid, int innprocs)
        : MultiGeneChain(inmyid, innprocs),
          modeltype("MULTIGENECODONM2A"),
          datafile(indatafile),
          treefile(intreefile) {
        blmode = inblmode;
        blsamplemode = inblsamplemode;
        nucmode = innucmode;
        purommode = inpurommode;
        dposommode = indposommode;
        purwmode = inpurwmode;
        poswmode = inposwmode;
        pihypermean = inpihypermean;
        pihyperinvconc = inpihyperinvconc;
        puromhypermean = inpuromhypermean;
        puromhyperinvconc = inpuromhyperinvconc;
        dposomhypermean = indposomhypermean;
        dposomhyperinvshape = indposomhyperinvshape;
        purwhypermean = inpurwhypermean;
        purwhyperinvconc = inpurwhyperinvconc;
        poswhypermean = inposwhypermean;
        poswhyperinvconc = inposwhyperinvconc;

        every = inevery;
        until = inuntil;
        writegenedata = inwritegenedata;
        name = inname;
        New(force);
    }

    //! \brief constructor for opening and restarting an already existing chain
    MultiGeneCodonM2aChain(string filename, int inmyid, int innprocs)
        : MultiGeneChain(inmyid, innprocs) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new MultiGeneCodonM2aModel(datafile, treefile, pihypermean, pihyperinvconc, myid,
                                           nprocs);
        GetModel()->SetAcrossGenesModes(blmode, nucmode, purommode, dposommode, purwmode, poswmode);
        GetModel()->SetBLSamplingMode(blsamplemode);
        GetModel()->SetMixtureHyperParameters(puromhypermean, puromhyperinvconc, dposomhypermean,
                                              dposomhyperinvshape, purwhypermean, purwhyperinvconc,
                                              poswhypermean, poswhyperinvconc);

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
            cerr << "Error : cannot find file : " << name << ".param\n";
            exit(1);
        }
        is >> modeltype;
        is >> datafile >> treefile;
        is >> writegenedata;
        is >> blmode >> nucmode >> dposommode >> purwmode >> poswmode;
        is >> pihypermean >> pihyperinvconc;
        is >> puromhypermean >> puromhyperinvconc;
        is >> dposomhypermean >> dposomhyperinvshape;
        is >> purwhypermean >> purwhyperinvconc;
        is >> poswhypermean >> poswhyperinvconc;

        blsamplemode = 0;
        int tmp;
        is >> tmp;
        if (tmp) {
            is >> blsamplemode;
            is >> tmp;
            if (tmp)    {
                cerr << "Error when reading model\n";
                exit(1);
            }
        }
        is >> every >> until >> size;

        if (modeltype == "MULTIGENECODONM2A") {
            model = new MultiGeneCodonM2aModel(datafile, treefile, pihypermean, pihyperinvconc,
                                               myid, nprocs);
            GetModel()->SetAcrossGenesModes(blmode, nucmode, purommode, dposommode, purwmode,
                                            poswmode);
            GetModel()->SetBLSamplingMode(blsamplemode);
            GetModel()->SetMixtureHyperParameters(
                puromhypermean, puromhyperinvconc, dposomhypermean, dposomhyperinvshape,
                purwhypermean, purwhyperinvconc, poswhypermean, poswhyperinvconc);
        } else {
            cerr << "Error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }

        if (!myid) {
            cerr << "allocate\n";
        }
        GetModel()->Allocate();

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
            param_os << blmode << '\t' << nucmode << '\t' << dposommode << '\t' << purwmode << '\t'
                     << poswmode << '\n';
            param_os << pihypermean << '\t' << pihyperinvconc << '\n';
            param_os << puromhypermean << '\t' << puromhyperinvconc << '\n';
            param_os << dposomhypermean << '\t' << dposomhyperinvshape << '\n';
            param_os << purwhypermean << '\t' << purwhyperinvconc << '\n';
            param_os << poswhypermean << '\t' << poswhyperinvconc << '\n';
            param_os << 1 << '\n';
            param_os << blsamplemode << '\n';
            param_os << 0 << '\n';
            param_os << every << '\t' << until << '\t' << size << '\n';
            GetModel()->MasterToStream(param_os);
        } else {
            GetModel()->SlaveToStream();
        }
    }

    void MakeFiles(int force) override {
        MultiGeneChain::MakeFiles(force);

        if (writegenedata) {
            ofstream pos((name + ".posw").c_str());
            ofstream omos((name + ".posom").c_str());
            if (writegenedata == 2) {
                ofstream siteos((name + ".sitepp").c_str());
            }
        }
    }

    void SavePoint() override {
        MultiGeneChain::SavePoint();
        if (writegenedata) {
            if (!myid) {
                ofstream posw_os((name + ".posw").c_str(), ios_base::app);
                GetModel()->TracePosWeight(posw_os);
                ofstream posom_os((name + ".posom").c_str(), ios_base::app);
                GetModel()->TracePosOm(posom_os);
                if (writegenedata == 2) {
                    ofstream pp_os((name + ".sitepp").c_str(), ios_base::app);
                    GetModel()->MasterTraceSitesPostProb(pp_os);
                }
            } else {
                if (writegenedata == 2) {
                    GetModel()->SlaveTraceSitesPostProb();
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

    int blockcounts[2] = {1, 3};
    MPI_Datatype types[2] = {MPI_DOUBLE, MPI_INT};
    MPI_Aint dtex, displacements[2];

    displacements[0] = (MPI_Aint)0;
    MPI_Type_extent(MPI_DOUBLE, &dtex);
    displacements[1] = dtex;
    MPI_Type_struct(2, blockcounts, displacements, types, &Propagate_arg);
    MPI_Type_commit(&Propagate_arg);

    MultiGeneCodonM2aChain *chain = 0;
    string name = "";

    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        chain = new MultiGeneCodonM2aChain(name, myid, nprocs);
    }

    // new chain
    else {
        string datafile = "";
        string treefile = "";

        double pihypermean = 0.1;
        double pihyperinvconc = 0.2;

        double puromhypermean = 0.5;
        double puromhyperinvconc = 0.5;
        int purommode = 1;

        double purwhypermean = 0.5;
        double purwhyperinvconc = 0.5;
        int purwmode = 1;

        double poswhypermean = 0.5;
        double poswhyperinvconc = 0.1;
        int poswmode = 1;

        double dposomhypermean = 1.0;
        double dposomhyperinvshape = 0.5;
        int dposommode = 1;

        int blmode = 1;
        int nucmode = 1;

        int blsamplemode = 0;

        int writegenedata = 1;

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
                } else if ((s == "-t") || (s == "-T")) {
                    i++;
                    treefile = argv[i];
                } else if (s == "-purom") {
                    purommode = 0;
                    i++;
                    string tmp = argv[i];
                    if (tmp != "uninf") {
                        puromhypermean = atof(argv[i]);
                        i++;
                        puromhyperinvconc = atof(argv[i]);
                    }
                } else if (s == "-dposom") {
                    dposommode = 0;
                    i++;
                    string tmp = argv[i];
                    if (tmp != "uninf") {
                        dposomhypermean = atof(argv[i]);
                        i++;
                        dposomhyperinvshape = atof(argv[i]);
                    }
                } else if (s == "-purw") {
                    purwmode = 0;
                    i++;
                    string tmp = argv[i];
                    if (tmp != "uninf") {
                        purwhypermean = atof(argv[i]);
                        i++;
                        purwhyperinvconc = atof(argv[i]);
                    }
                } else if (s == "-posw") {
                    poswmode = 0;
                    i++;
                    string tmp = argv[i];
                    if (tmp != "uninf") {
                        poswhypermean = atof(argv[i]);
                        i++;
                        poswhyperinvconc = atof(argv[i]);
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
            cerr << '\n';
            cerr << "mpirun -np <n> multigenecodonm8 -d <alignment_list> -t <tree> "
                    "[-bl {shared|shrunken|ind} -nucrates {shared|shrunken|ind}] "
                    "<chainname> \n";
            cerr << '\n';
            cerr << "chain options:\n";
            cerr << "\t-f: force overwrite of already existing chain\n";
            cerr << "\t-x <every> <until>: saving frequency and stopping time "
                    "(default: every = 1, until = -1)\n";
            cerr << "\t-g: without gene-specific output files (.posw and .posom)\n";
            cerr << "\t+g: with gene-specific output files (.posw and .posom)\n";
            cerr << "\t+G: with gene- and site-specific output files\n";
            cerr << "\tin all cases, complete information about chain state is saved "
                    "in .chain and .param files\n";
            cerr << "\t.chain: one line for each cycle\n";
            cerr << "\t.param: state at the end of last cycle\n";
            cerr << '\n';
            cerr << "model options:\n";
            cerr << "\t-bl {shared|shrunken|ind}: shrinkage mode for branch lengths\n";
            cerr << "\t-nucrates {shared|shrunken|ind}: shrinkage mode for "
                    "nucleotide substitution rates\n";
            cerr << '\n';
            exit(1);
        }

        chain = new MultiGeneCodonM2aChain(
            datafile, treefile, blmode, blsamplemode, nucmode, purommode, dposommode, purwmode, poswmode,
            pihypermean, pihyperinvconc, puromhypermean, puromhyperinvconc, dposomhypermean,
            dposomhyperinvshape, purwhypermean, purwhyperinvconc, poswhypermean, poswhyperinvconc,
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
