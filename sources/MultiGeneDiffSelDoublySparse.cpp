#include <cmath>
#include <fstream>
#include "Chrono.hpp"
#include "MultiGeneChain.hpp"
#include "MultiGeneDiffSelDoublySparseModel.hpp"

using namespace std;

MPI_Datatype Propagate_arg;

/**
 * \brief A MultiGeneChain object for running an MCMC under
 * MultiGeneDiffSelDoublySparseModel
 */
class MultiGeneDiffSelDoublySparseChain : public MultiGeneChain {
  private:
    // Chain parameters
    string modeltype, datafile, treefile;
    int ncond;
    int nlevel;
    int codonmodel;
    double epsilon;
    double fitnessshape;
    int fitnesscentermode;
    int blmode, nucmode, shiftmode;
    double pihypermean, pihyperinvconc;
    double shiftprobmean, shiftprobinvconc;
    int burnin;
    int writegenedata;

  public:
    MultiGeneDiffSelDoublySparseModel *GetModel() {
        return static_cast<MultiGeneDiffSelDoublySparseModel *>(model);
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
    MultiGeneDiffSelDoublySparseChain(string indatafile, string intreefile, int inncond,
                                      int innlevel, int incodonmodel, double inepsilon,
                                      double infitnessshape, int infitnesscentermode, int inblmode, int innucmode, int inshiftmode,
                                      double inpihypermean, double inpihyperinvconc, double inshiftprobmean, double inshiftprobinvconc,
                                      int inburnin, int inevery, int inuntil, int insaveall, int inwritegenedata,
                                      string inname, int force, int inmyid, int innprocs)
        : MultiGeneChain(inmyid, innprocs),
          modeltype("MULTIGENEDIFFSELDSPARSE"),
          datafile(indatafile),
          treefile(intreefile),
          ncond(inncond),
          nlevel(innlevel),
          codonmodel(incodonmodel),
          epsilon(inepsilon),
          fitnessshape(infitnessshape),
          fitnesscentermode(infitnesscentermode), blmode(inblmode), nucmode(innucmode),
          shiftmode(inshiftmode), pihypermean(inpihypermean), pihyperinvconc(inpihyperinvconc), shiftprobmean(inshiftprobmean), shiftprobinvconc(inshiftprobinvconc) {
        burnin = inburnin;
        every = inevery;
        until = inuntil;
        saveall = insaveall;
        writegenedata = inwritegenedata;
        name = inname;
        New(force);
    }

    //! \brief constructor for opening and restarting an already existing chain
    MultiGeneDiffSelDoublySparseChain(string filename, int inmyid, int innprocs)
        : MultiGeneChain(inmyid, innprocs) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new MultiGeneDiffSelDoublySparseModel(datafile, treefile, ncond, nlevel, codonmodel,
                                                      epsilon, fitnessshape, blmode, nucmode, shiftmode,
                                                      pihypermean, pihyperinvconc, shiftprobmean, shiftprobinvconc,
                                                      myid, nprocs);
        if (burnin) {
            GetModel()->SetWithToggles(0);
        } else {
            GetModel()->SetWithToggles(1);
        }
        GetModel()->SetFitnessCenterMode(fitnesscentermode);
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
        is >> epsilon >> fitnessshape;
        is >> fitnesscentermode;
        is >> blmode >> nucmode >> shiftmode;
        is >> pihypermean >> pihyperinvconc;
        is >> shiftprobmean >> shiftprobinvconc;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> burnin;
        is >> every >> until >> saveall >> writegenedata >> size;

        if (modeltype == "MULTIGENEDIFFSELDSPARSE") {
            model = new MultiGeneDiffSelDoublySparseModel(
                datafile, treefile, ncond, nlevel, codonmodel, epsilon, fitnessshape, blmode, nucmode, shiftmode, pihypermean, pihyperinvconc, shiftprobmean, shiftprobinvconc, myid, nprocs);
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
            param_os << epsilon << '\t' << fitnessshape << '\n';
            param_os << fitnesscentermode << '\n';
            param_os << blmode << '\t' << nucmode << '\t' << shiftmode << '\n';
            param_os << pihypermean << '\t' << pihyperinvconc << '\t' << shiftprobmean << '\t' << shiftprobinvconc << '\n';
            param_os << 0 << '\n';
            param_os << burnin << '\t';
            param_os << every << '\t' << until << '\t' << saveall << '\t' << writegenedata << '\t' << size << '\n';
            GetModel()->MasterToStream(param_os);
        } else {
            GetModel()->SlaveToStream();
        }
        if (size == burnin) {
            GetModel()->SetWithToggles(1);
        }
    }

    void MakeFiles(int force) override {
        MultiGeneChain::MakeFiles(force);
        cerr << writegenedata << '\t' << ncond << '\n';
        if (writegenedata) {
            if (ncond > 1)  {
                for (int k = 0; k < ncond; k++) {
                    ostringstream s;
                    s << name << "_" << k;
                    if (k) {
                        ofstream pos((s.str() + ".geneshiftprob").c_str());
                        ofstream cos((s.str() + ".geneshiftcounts").c_str());
                        if (writegenedata == 2) {
                            ofstream tos((s.str() + ".shifttoggle").c_str());
                        }
                    }
                    if (writegenedata == 2) {
                        ofstream fos((s.str() + ".fitness").c_str());
                    }
                }
                ofstream os((name + ".genemaskcounts").c_str());
            }
            else    {
                ofstream os((name + ".geneom").c_str());
            }
        }
    }

    void SavePoint() override {
        MultiGeneChain::SavePoint();
        if (writegenedata) {
            if (ncond > 1)  {
                if (!myid) {
                    GetModel()->MasterTraceSiteStats(name, writegenedata);
                } else {
                    GetModel()->SlaveTraceSiteStats(writegenedata);
                }
            }
            else    {
                if (! myid) {
                    ofstream os((name + ".geneom").c_str(), ios_base::app);
                    GetModel()->TracePredictedDNDS(os);
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
    MultiGeneDiffSelDoublySparseChain *chain = 0;

    // command syntax
    if (argc == 1)  {
        if (! myid)	{
            cerr << '\n';
            cerr << "The multi-gene version of the doubly-sparse differential selection model.\n";
            cerr << "see diffseldsparse for a more detailed description of the single-gene version.\n";
            cerr << "\n";
            cerr << "the key gene-specific parameters, for which shrinkage across genes is implemented, are:\n";
            cerr << " - branch lengths\n";
            cerr << " - nucleotide mutation rates\n";
            cerr << " - shiftprob_gk, for gene g, condition k=1..K-1 (specifying the probability that an amino-acid at any site in gene g undergoes a shift in condition k)\n";
            cerr << "concerning shiftprob_gk, the prior distribution is:\n";
            cerr << " - with prob (1-pi_k), shiftprob_gk = 0 (i.e. the gene does not have any site showing differential effect in condition k)\n";
            cerr << " - with prob pi_k, shiftprob_gk ~ Beta(shiftprobhypermean_k, shiftprobhyperinvconc_k)\n";
            cerr << "thus, the probability that a gene is under differential selection is given by the post prob that shiftprob_k > 0\n";
            cerr << "by default, the hyperparameters pi_k, shiftprobhypermean_k, shiftprobhyperinvconc_k are estimated across genes\n";
            cerr << "they can also be fixed a priori\n";
            cerr << '\n';
            cerr << "command: mpirun -np <n> multigenediffseldsparse -d <alignment_list> -t <tree> -ncond <ncond> <chainname>\n";
            cerr << '\n';
            cerr << "chain options:\n";
            cerr << "\t-f: force overwrite of already existing chain\n";
            cerr << "\t-x <every> <until>: saving frequency and stopping time "
                    "(default: every = 1, until = -1)\n";
            cerr << "\t-g: without gene-specific output files (.geneshiftprob and geneshiftcounts)\n";
            cerr << "\t+g: with gene-specific output files\n";
            cerr << "\t+G: with gene- and site-specific output files\n";
            // cerr << "\tin all cases, complete information about chain state is saved "
            //         "in .chain and .param files\n";
            // cerr << "\t.chain: one line for each cycle\n";
            // cerr << "\t.param: state at the end of last cycle\n";
            cerr << '\n';
            cerr << "model options:\n";
            cerr << "\t-ncond <ncond>:  specify number of conditions\n";
            cerr << "\t-bl {shrunken|ind}: shrinkage mode for branch lengths\n";
            cerr << "\t-nucrates {shrunken|ind}: shrinkage mode for nucleotide substitution rates\n";
            cerr << "\t-pi <hypermean> <hyperinvconc>: set parameters of beta hyperprior for pi_k for all k=1..K-1\n";
            cerr << "\t                                (default: hypermean = 0.1, hyperinvconc = 0.1)\n";
            cerr << "\t-shiftprob <hypermean> <hyperinvconc>: set values of shiftprobhypermean_k and shiftprobhyperinvconc_k for all k=1..K-1\n";
            cerr << "\t                                       (default: hypermean = 0.1, hyperinvconc = 0.1)\n";
            cerr << "\t                                       (if hyperinvconc == 0, then shiftprob_gk is fixed to hypermean for all genes and conditions)\n";
            cerr << '\n';
        }
        MPI_Finalize();
        exit(0);
    }

    // starting a chain from existing files
    else if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        chain = new MultiGeneDiffSelDoublySparseChain(name, myid, nprocs);
    }

    // new chain
    else {
        string datafile = "";
        string treefile = "";
        int ncond = 1;
        int nlevel = 1;
        int codonmodel = 1;
        int force = 1;
        int burnin = 0;
        int every = 1;
        int until = -1;
        int saveall = 1;
        int writegenedata = 1;
        double fitnessshape = 20;
        double epsilon = 0.001;
        int fitnesscentermode = 3;
        int blmode = 1;
        int nucmode = 1;
        int shiftmode = 1;
        double pihypermean = 0.1;
        double pihyperinvconc = 0.1;
        double shiftprobmean = 0.1;
        double shiftprobinvconc = 0.1;

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
                } else if (s == "-pi")    {
                    i++;
                    pihypermean = atof(argv[i]);
                    i++;
                    pihyperinvconc = atof(argv[i]);
                } else if (s == "-shiftprob")   {
                    i++;
                    string tmp = argv[i];
                    if (tmp == "shrunken") {
                        shiftmode = 1;
                    } else  {
                        shiftmode = 0;
                        if (tmp != "uninf")  {
                            shiftprobmean = atof(argv[i]);
                            i++;
                            shiftprobinvconc = atof(argv[i]);
                        }
                    }
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
                // burnin de-activated
                /*
                } else if (s == "-b") {
                    i++;
                    burnin = atoi(argv[i]);
                */
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
            if (! myid) {
                cerr << "error in command\n";
                cerr << '\n';
            }
            MPI_Finalize();
            exit(0);
        }

        chain = new MultiGeneDiffSelDoublySparseChain(
            datafile, treefile, ncond, nlevel, codonmodel, epsilon, fitnessshape, fitnesscentermode, blmode, nucmode, shiftmode,
            pihypermean, pihyperinvconc, shiftprobmean, shiftprobinvconc,
            burnin, every, until, saveall, writegenedata, name, force, myid, nprocs);
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
