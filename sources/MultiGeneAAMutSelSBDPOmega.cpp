#include <cmath>
#include <fstream>
#include "MultiGeneChain.hpp"
#include "MultiGeneAAMutSelSBDPOmegaModel.hpp"

using namespace std;

MPI_Datatype Propagate_arg;

class MultiGeneAAMutSelSBDPOmegaChain : public MultiGeneChain  {

  private:
    // Chain parameters
    string modeltype, datafile, treefile;
    int Ncat;

  public:
    MultiGeneAAMutSelSBDPOmegaModel* GetModel() {
        return static_cast<MultiGeneAAMutSelSBDPOmegaModel*>(model);
    }

    string GetModelType() override { return modeltype; }

    MultiGeneAAMutSelSBDPOmegaChain(string indatafile, string intreefile, int inNcat, int inevery, int inuntil, string inname, int force, int inmyid, int innprocs) : MultiGeneChain(inmyid,innprocs), modeltype("MULTIGENEAAMUTSELSBDPOMEGA"), datafile(indatafile), treefile(intreefile), Ncat(inNcat) {
        every = inevery;
        until = inuntil;
        name = inname;
        New(force);
    }

    MultiGeneAAMutSelSBDPOmegaChain(string filename, int inmyid, int innprocs) : MultiGeneChain(inmyid,innprocs) {
        name = filename;
        Open();
        if (! myid) {
            Save();
        }
    }

    void New(int force) override {
        model = new MultiGeneAAMutSelSBDPOmegaModel(datafile,treefile,Ncat,myid,nprocs);
        if (! myid) {
            cerr << " -- master allocate\n";
        }
        GetModel()->Allocate();
        if (! myid) {
            cerr << " -- master unfold\n";
        }
        GetModel()->Unfold();

        if (! myid) {
            cerr << "-- Reset" << endl;
            Reset(force);
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
        is >> Ncat;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "MULTIGENEAAMUTSELSBDPOMEGA") {
            model = new MultiGeneAAMutSelSBDPOmegaModel(datafile,treefile,Ncat,myid,nprocs);
        } else {
            cerr << "-- Error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }
        GetModel()->Allocate();
        if (! myid) {
            model->FromStream(is);
            // broadcast parameter
        }
        else    {
            // receive parameter
        }
        model->Update();
        GetModel()->Unfold();
        if (! myid) {
            cerr << size << " points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
            model->Trace(cerr);
        }
    }

    void Save() override {
        if (myid)   {
            cerr << "error: slave in MultiGeneAAMutSelSBDPOmegaChain::Save\n";
            exit(1);
        }
        ofstream param_os((name + ".param").c_str());
        param_os << GetModelType() << '\n';
        param_os << datafile << '\t' << treefile << '\n';
        param_os << Ncat << '\n';
        param_os << 0 << '\n';
        param_os << every << '\t' << until << '\t' << size << '\n';
        model->ToStream(param_os);
    }

    void Monitor() override {
        Chain::Monitor();
        ofstream pos((name + ".geneom").c_str(), ios_base::app);
        GetModel()->TraceGeneOmegas(pos);
    }

    void MakeFiles(int force) override {
        Chain::MakeFiles(force);
        ofstream pos((name + ".geneom").c_str());
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

    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        string name = argv[1];
        MultiGeneAAMutSelSBDPOmegaChain* chain = new MultiGeneAAMutSelSBDPOmegaChain(name,myid,nprocs);
        if (!myid)  {
            cerr << "chain " << name << " started\n";
        }
        if (!myid)  {
            chain->Start();
            cerr << "chain " << name << " stopped\n";
            cerr << chain->GetSize() << " points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
            chain->GetModel()->Trace(cerr);
        }
    }

    // new chain
    else    {
        string datafile = "";
        string treefile = "";
        int Ncat = -1;
        string name = "";
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
                else if (s == "-ncat") {
                    i++;
                    Ncat = atoi(argv[i]);
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
            cerr << "multigeneaamutselhdp -d <list> -t <tree> -ncat <ncat> <chainname> \n";
            cerr << '\n';
            exit(1);
        }

        MultiGeneAAMutSelSBDPOmegaChain* chain = new MultiGeneAAMutSelSBDPOmegaChain(datafile,treefile,Ncat,every,until,name,force,myid,nprocs);
        if (! myid) {
            cerr << "chain " << name << " started\n";
        }
        chain->Start();
        if (! myid) {
            cerr << "chain " << name << " stopped\n";
            cerr << chain->GetSize() << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
            chain->GetModel()->Trace(cerr);
        }
    }

	MPI_Finalize();
}


