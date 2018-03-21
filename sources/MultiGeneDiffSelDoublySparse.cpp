#include <cmath>
#include <fstream>
#include "MultiGeneChain.hpp"
#include "MultiGeneDiffSelDoublySparseModel.hpp"
#include "Chrono.hpp"

using namespace std;

MPI_Datatype Propagate_arg;

/**
 * \brief A MultiGeneChain object for running an MCMC under MultiGeneDiffSelDoublySparseModel
 */
class MultiGeneDiffSelDoublySparseChain : public MultiGeneChain  {

  private:
    // Chain parameters
    string modeltype, datafile, treefile;
    int ncond;
    int nlevel;
    int codonmodel;
    double epsilon;
    double fitnessshape;
    int burnin;
    int writegenedata;

  public:
    MultiGeneDiffSelDoublySparseModel* GetModel() {
        return static_cast<MultiGeneDiffSelDoublySparseModel*>(model);
    }

    string GetModelType() override { return modeltype; }

    //! \brief constructor for a new MCMC
    //!
    //! \param indatafile: name of file contanining sequence alignment
    //! \param intreefile: name of file contaning tree (with branch names giving the allocation of branches to the conditions)
    //! \param inncond: number of conditions
    //! \param innlevel: number of levels of the model
    //! \param incodonmodel: type of codon substitution process (1: mutation-selection, 0: square-root)
    //! \param inevery: thinning factor
    //! \param inuntil: maximum MCMC sample size
    //! \param insaveall: if 1, then, save all information about each configuration visited during MCMC (into .chain file)
    //! \param inwritegenedata: if 1, then trace gene- and condition-specific shift probabilities in separate files; if 2, then also trace site-specific shift probabilities
    //! \param name: base name for all files related to this MCMC run
    //! \param force: overwrite existing files with same name
    //! \param inmyid, int innprocs: process id and total number of MPI processes
    MultiGeneDiffSelDoublySparseChain(string indatafile, string intreefile, int inncond, int innlevel, int incodonmodel, double inepsilon, double infitnessshape, int inburnin, int inevery, int inuntil, int insaveall, int inwritegenedata, string inname, int force, int inmyid, int innprocs) : MultiGeneChain(inmyid,innprocs), modeltype("MULTIGENEDIFFSELSPARSE"), datafile(indatafile), treefile(intreefile), ncond(inncond), nlevel(innlevel), codonmodel(incodonmodel), epsilon(inepsilon), fitnessshape(infitnessshape) {
        burnin = inburnin;
        every = inevery;
        until = inuntil;
        saveall = insaveall;
        writegenedata = inwritegenedata;
        name = inname;
        New(force);
    }

    //! \brief constructor for opening and restarting an already existing chain
    MultiGeneDiffSelDoublySparseChain(string filename, int inmyid, int innprocs) : MultiGeneChain(inmyid,innprocs) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new MultiGeneDiffSelDoublySparseModel(datafile,treefile,ncond,nlevel,codonmodel,epsilon,fitnessshape,myid,nprocs);
        if (burnin) {
            GetModel()->SetWithToggles(0);
        }
        else    {
            GetModel()->SetWithToggles(1);
        }
        if (! myid) {
            cerr << " -- master allocate\n";
        }
        GetModel()->Allocate();
        if (! myid) {
            cerr << " -- master unfold\n";
        }
        GetModel()->Update();
        Reset(force);
        if (! myid) {
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
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> burnin;
        is >> every >> until >> saveall >> writegenedata >> size;

        if (modeltype == "MULTIGENEDIFFSELSPARSE") {
            model = new MultiGeneDiffSelDoublySparseModel(datafile,treefile,ncond,nlevel,codonmodel,epsilon,fitnessshape,myid,nprocs);
        } else {
            cerr << "-- Error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }
        if (size < burnin)  {
            GetModel()->SetWithToggles(0);
        }
        else    {
            GetModel()->SetWithToggles(1);
        }
        GetModel()->Allocate();
        GetModel()->FromStream(is);
        GetModel()->Update();
        if (! myid) {
            cerr << size << " points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
            model->Trace(cerr);
        }
    }

    void Save() override {
        if (!myid)   {
            ofstream param_os((name + ".param").c_str());
            param_os << GetModelType() << '\n';
            param_os << datafile << '\t' << treefile << '\n';
            param_os << ncond << '\t' << nlevel << '\t' << codonmodel << '\n';
            param_os << epsilon << '\t' << fitnessshape << '\n';
            param_os << 0 << '\n';
            param_os << burnin << '\t';
            param_os << every << '\t' << until << '\t' << saveall << '\t' << writegenedata << '\t' << size << '\n';
            GetModel()->MasterToStream(param_os);
        }
        else    {
            GetModel()->SlaveToStream();
        }
        if (size == burnin)  {
            GetModel()->SetWithToggles(1);
        }
    }

    void MakeFiles(int force) override  {
        MultiGeneChain::MakeFiles(force);
        if (writegenedata)  {
            for (int k=0; k<ncond; k++) {
                ostringstream s;
                s << name << "_" << k;
                if (k)  {
                    ofstream pos((s.str() + ".geneshiftprob").c_str());
                    if (writegenedata == 2) {
                        ofstream tos((s.str() + ".shifttoggle").c_str());
                    }
                }
                if (writegenedata == 2) {
                    ofstream fos((s.str() + ".fitness").c_str());
                }
            }
        }
    }

    void SavePoint() override   {
        MultiGeneChain::SavePoint();
        if (writegenedata)  {
            if (! myid) {
                GetModel()->MasterTraceSiteStats(name,writegenedata);
            }
            else    {
                GetModel()->SlaveTraceSiteStats(writegenedata);
            }
        }
    }
};

int main(int argc, char* argv[])	{

	Chrono chrono;
    chrono.Start();

	int myid  = 0;
	int nprocs = 0;

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

	int blockcounts[2] = {1,3};
	MPI_Datatype types[2] = {MPI_DOUBLE,MPI_INT};
	MPI_Aint dtex,displacements[2];
	
	displacements[0] = (MPI_Aint) 0;
	MPI_Type_extent(MPI_DOUBLE,&dtex);
	displacements[1] = dtex;
	MPI_Type_struct(2,blockcounts,displacements,types,&Propagate_arg);
	MPI_Type_commit(&Propagate_arg); 

    string name = "";
    MultiGeneDiffSelDoublySparseChain* chain = 0;

    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        chain = new MultiGeneDiffSelDoublySparseChain(name,myid,nprocs);
    }

    // new chain
    else    {
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
        double epsilon = -1;

        try	{

            if (argc == 1)	{
                throw(0);
            }

            int i = 1;
            while (i < argc)	{
                string s = argv[i];

                if (s == "-d")	{
                    i++;
                    datafile = argv[i];
                }
                else if ((s == "-t") || (s == "-T"))	{
                    i++;
                    treefile = argv[i];
                }
                else if (s == "-f")	{
                    force = 1;
                }
                else if (s == "-s") {
                    saveall = 0;
                }
                else if (s == "+s") {
                    saveall = 1;
                }
                else if (s == "-ncond")	{
                    i++;
                    ncond = atoi(argv[i]);
                }
                else if (s == "-nlevel")	{
                    i++;
                    nlevel = atoi(argv[i]);
                }
                else if (s == "-g")  {
                    writegenedata = 0;
                }
                else if (s == "+g")  {
                    writegenedata = 1;
                }
                else if (s == "+G")  {
                    writegenedata = 2;
                }
                else if (s == "-b") {
                    i++;
                    burnin = atoi(argv[i]);
                }
                else if ( (s == "-x") || (s == "-extract") )	{
                    i++;
                    if (i == argc) throw(0);
                    every = atoi(argv[i]);
                    i++;
                    if (i == argc) throw(0);
                    until = atoi(argv[i]);
                }
                else	{
                    if (i != (argc -1))	{
                        throw(0);
                    }
                    name = argv[i];
                }
                i++;
            }
            if ((datafile == "") || (treefile == "") || (name == ""))	{
                throw(0);
            }
        }
        catch(...)	{
            cerr << "error in command\n";
            cerr << '\n';
            exit(1);
        }

        chain = new MultiGeneDiffSelDoublySparseChain(datafile,treefile,ncond,nlevel,codonmodel,epsilon,fitnessshape,burnin,every,until,saveall,writegenedata,name,force,myid,nprocs);
    }

    chrono.Stop();
    if (! myid) {
        cout << "total time to set things up: " << chrono.GetTime() << '\n';
    }
    chrono.Reset();
    chrono.Start();
    if (! myid) {
        cerr << "chain " << name << " started\n";
    }
    chain->Start();
    if (! myid) {
        cerr << "chain " << name << " stopped\n";
        cerr << chain->GetSize() << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
        chain->GetModel()->Trace(cerr);
    }
    chrono.Stop();
    if (! myid) {
        cout << "total time in MCMC: " << chrono.GetTime() << '\n';
        cout << "total time in master moves: " << chain->GetModel()->GetMasterMoveTime() << '\n';
        cout << "mean total time in slave moves: " << chain->GetModel()->GetSlaveMoveTime() << '\n';
        cout << "mean total time in substitution mapping: " << chain->GetModel()->GetSlaveMapTime() << '\n';
    }

	MPI_Finalize();
}


