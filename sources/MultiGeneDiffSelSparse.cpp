#include <cmath>
#include <fstream>
#include "MultiGeneChain.hpp"
#include "MultiGeneDiffSelSparseModel.hpp"
#include "Chrono.hpp"

using namespace std;

MPI_Datatype Propagate_arg;

class MultiGeneDiffSelSparseChain : public MultiGeneChain  {

  private:
    // Chain parameters
    string modeltype, datafile, treefile;
    int ncond;
    int nlevel;
    int codonmodel;

  public:
    MultiGeneDiffSelSparseModel* GetModel() {
        return static_cast<MultiGeneDiffSelSparseModel*>(model);
    }

    string GetModelType() override { return modeltype; }

    MultiGeneDiffSelSparseChain(string indatafile, string intreefile, int inncond, int innlevel, int incodonmodel, int inevery, int inuntil, string inname, int force, int inmyid, int innprocs) : MultiGeneChain(inmyid,innprocs), modeltype("MULTIGENEDIFFSELSPARSE"), datafile(indatafile), treefile(intreefile), ncond(inncond), nlevel(innlevel), codonmodel(incodonmodel) {
        every = inevery;
        until = inuntil;
        name = inname;
        New(force);
    }

    MultiGeneDiffSelSparseChain(string filename, int inmyid, int innprocs) : MultiGeneChain(inmyid,innprocs) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new MultiGeneDiffSelSparseModel(datafile,treefile,ncond,nlevel,codonmodel,myid,nprocs);
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
            cerr << "-- initial ln prob = " << GetModel()->GetLogProb() << "\n";
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
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "MULTIGENEDIFFSELSPARSE") {
            model = new MultiGeneDiffSelSparseModel(datafile,treefile,ncond,nlevel,codonmodel,myid,nprocs);
        } else {
            cerr << "-- Error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
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
            param_os << 0 << '\n';
            param_os << every << '\t' << until << '\t' << size << '\n';
            GetModel()->MasterToStream(param_os);
        }
        else    {
            GetModel()->SlaveToStream();
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
    MultiGeneDiffSelSparseChain* chain = 0;

    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        chain = new MultiGeneDiffSelSparseChain(name,myid,nprocs);
    }

    // new chain
    else    {
        string datafile = "";
        string treefile = "";
        int ncond = 1;
        int nlevel = 1;
        int codonmodel = 1;
        int force = 1;
        int every = 1;
        int until = -1;

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
                else if (s == "-ncond")	{
                    i++;
                    ncond = atoi(argv[i]);
                }
                else if (s == "-nlevel")	{
                    i++;
                    nlevel = atoi(argv[i]);
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

        chain = new MultiGeneDiffSelSparseChain(datafile,treefile,ncond,nlevel,codonmodel,every,until,name,force,myid,nprocs);
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
        cout << "total time to run: " << chrono.GetTime() << '\n';
    }

	MPI_Finalize();
}


