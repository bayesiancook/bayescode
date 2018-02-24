#include <cmath>
#include <fstream>
#include "MultiGeneChain.hpp"
#include "MultiGeneBranchOmegaModel.hpp"

using namespace std;

MPI_Datatype Propagate_arg;

class MultiGeneBranchOmegaChain : public MultiGeneChain  {

  private:
    // Chain parameters
    string modeltype, datafile, treefile;

  public:
    MultiGeneBranchOmegaModel* GetModel() {
        return static_cast<MultiGeneBranchOmegaModel*>(model);
    }

    string GetModelType() override { return modeltype; }

    MultiGeneBranchOmegaChain(string indatafile, string intreefile, int inevery, int inuntil, string inname, int force, int inmyid, int innprocs) : MultiGeneChain(inmyid,innprocs), modeltype("MULTIGENEBRANCHOMEGA"), datafile(indatafile), treefile(intreefile) {
        every = inevery;
        until = inuntil;
        name = inname;
        New(force);
    }

    MultiGeneBranchOmegaChain(string filename, int inmyid, int innprocs) : MultiGeneChain(inmyid,innprocs) {
        name = filename;
        Open();
	Save();
    }

    void New(int force) override {
        model = new MultiGeneBranchOmegaModel(datafile,treefile,myid,nprocs);
        if (! myid) {
            cerr << "allocate\n";
        }
        GetModel()->Allocate();
        GetModel()->Update();
        Reset(force);

        if (! myid) {
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
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "MULTIGENEBRANCHOMEGA") {
            model = new MultiGeneBranchOmegaModel(datafile,treefile,myid,nprocs);
        } else {
            cerr << "error when opening file " << name
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
            param_os << 0 << '\n';
            param_os << every << '\t' << until << '\t' << size << '\n';
            GetModel()->MasterToStream(param_os);
        }
        else    {
            GetModel()->SlaveToStream();
        }
    }

    void MakeFiles(int force) override  {
        Chain::MakeFiles(force);
        ofstream gos((name + ".gene").c_str());
        ofstream bos((name + ".branch").c_str());
        ofstream bgos((name + ".branchgene").c_str());

        ofstream tos((name + ".branchindices").c_str());
        GetModel()->GetTree()->ToStreamWithBranchIndex(tos);

        ofstream nameos((name + ".genelist").c_str());
        GetModel()->PrintGeneList(nameos);
        nameos.close();
    }

    void SavePoint() override   {
        Chain::SavePoint();
        if (!myid)  {
            ofstream gos((name + ".gene").c_str(),ios_base::app);
            GetModel()->PrintGeneEffects(gos);
            ofstream bos((name + ".branch").c_str(),ios_base::app);
            GetModel()->PrintBranchEffects(bos);
            ofstream bgos((name + ".branchgene").c_str(),ios_base::app);
            GetModel()->PrintDeviations(bgos);
        }
    }
};

int main(int argc, char* argv[])	{

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

    if (nprocs <= 1)    {
        cerr << "error: should run the program with at least 2 cores\n";
        exit(1);
    }
        
    MultiGeneBranchOmegaChain* chain = 0;
    string name = "";

    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        chain = new MultiGeneBranchOmegaChain(name,myid,nprocs);
    }

    // new chain
    else    {
        string datafile = "";
        string treefile = "";
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
            cerr << "globom -d <alignment> -t <tree> <chainname> \n";
            cerr << '\n';
            exit(1);
        }

        chain = new MultiGeneBranchOmegaChain(datafile,treefile,every,until,name,force,myid,nprocs);
    }

    if (! myid) {
        cerr << "chain " << name << " started\n";
    }
    chain->Start();
    if (! myid) {
        cerr << "chain " << name << " stopped\n";
        cerr << chain->GetSize() << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
        chain->GetModel()->Trace(cerr);
    }

	MPI_Finalize();
}


