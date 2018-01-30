#include <cmath>
#include <fstream>
#include "MultiGeneChain.hpp"
#include "MultiGeneAAMutSelDSBDPOmegaModel.hpp"

using namespace std;

MPI_Datatype Propagate_arg;

class MultiGeneAAMutSelDSBDPOmegaChain : public MultiGeneChain  {

  private:
    // Chain parameters
    string modeltype, datafile, treefile;
    int Ncat;
    int baseNcat;
    int blmode, nucmode, basemode, omegamode;

  public:
    MultiGeneAAMutSelDSBDPOmegaModel* GetModel() {
        return static_cast<MultiGeneAAMutSelDSBDPOmegaModel*>(model);
    }

    string GetModelType() override { return modeltype; }

    MultiGeneAAMutSelDSBDPOmegaChain(string indatafile, string intreefile, int inNcat, int inbaseNcat, int inblmode, int innucmode, int inbasemode, int inomegamode, int inevery, int inuntil, string inname, int force, int inmyid, int innprocs) : MultiGeneChain(inmyid,innprocs), modeltype("MULTIGENEAAMUTSELDSBDPOMEGA"), datafile(indatafile), treefile(intreefile), Ncat(inNcat), baseNcat(inbaseNcat), blmode(inblmode), nucmode(innucmode), basemode(inbasemode), omegamode(inomegamode) {
        every = inevery;
        until = inuntil;
        name = inname;
        New(force);
    }

    MultiGeneAAMutSelDSBDPOmegaChain(string filename, int inmyid, int innprocs) : MultiGeneChain(inmyid,innprocs) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new MultiGeneAAMutSelDSBDPOmegaModel(datafile,treefile,Ncat,baseNcat,blmode,nucmode,basemode,omegamode,myid,nprocs);
        if (! myid) {
            cerr << " -- allocate\n";
        }
        GetModel()->Allocate();
        if (! myid) {
            cerr << " -- update\n";
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
        is >> Ncat >> baseNcat;
        is >> blmode >> nucmode >> basemode >> omegamode;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "MULTIGENEAAMUTSELDSBDPOMEGA") {
            model = new MultiGeneAAMutSelDSBDPOmegaModel(datafile,treefile,Ncat,baseNcat,blmode,nucmode,basemode,omegamode,myid,nprocs);
        } else {
            cerr << "-- Error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }

        GetModel()->Allocate();
        model->FromStream(is);
        GetModel()->Update();

        if (! myid) {
            cerr << size << " points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
            model->Trace(cerr);
        }
    }

    void Save() override {
        if (! myid) {
            ofstream param_os((name + ".param").c_str());
            param_os << GetModelType() << '\n';
            param_os << datafile << '\t' << treefile << '\n';
            param_os << Ncat << '\t' << baseNcat << '\n';
            param_os << blmode << '\t' << nucmode << '\t' << basemode << '\t' << omegamode << '\n';
            param_os << 0 << '\n';
            param_os << every << '\t' << until << '\t' << size << '\n';
            GetModel()->MasterToStream(param_os);
        }
        else    {
            GetModel()->SlaveToStream();
        }
    }

    void Monitor() override {
        if (myid)   {
            cerr << "error: in specialized monitor\n";
            exit(1);
        }
        Chain::Monitor();
        // ofstream trace_os((name + ".basemix").c_str(), ios_base::app);
        ofstream trace_os((name + ".basemix").c_str());
        GetModel()->TraceMixture(trace_os);
        ofstream logo_os((name + ".basemixlogo").c_str());
        GetModel()->PrintBaseMixtureLogo(logo_os);
        ofstream sample_os((name + ".basesamplelogo").c_str());
        GetModel()->PrintBaseSampleLogo(sample_os);
    }

    void MakeFiles(int force) override  {
        if (myid)   {
            cerr << "error: in specialized makefiles\n";
            exit(1);
        }
        Chain::MakeFiles(force);
        ofstream trace_os((name + ".basemix").c_str());
        ofstream logo_os((name + ".basemixlogo").c_str());
        ofstream samplelogo_os((name + ".basesamplelogo").c_str());
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
        MultiGeneAAMutSelDSBDPOmegaChain* chain = new MultiGeneAAMutSelDSBDPOmegaChain(name,myid,nprocs);
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
        int baseNcat = -1;
        string name = "";
        int force = 1;
        int every = 1;
        int until = -1;

        int blmode = 2;
        int nucmode = 0;
        int basemode = 2;
        int omegamode = 3;

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
                else if (s == "-basencat") {
                    i++;
                    baseNcat = atoi(argv[i]);
                }
                else if (s == "-basemix")    {
                    i++;
                    string tmp = argv[i];
                    if (tmp == "shared")    {
                        basemode = 2;
                    }
                    else if (tmp == "independent")  {
                        basemode = 0;
                    }
                    else    {
                        cerr << "error: does not recognize command after -basemix\n";
                        exit(1);
                    }
                }
                else if (s == "-fixomega")  {
                    omegamode = 3;
                }
                else if (s == "-freeomega") {
                    omegamode = 1;
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
            cerr << "multigeneaamutselddp -d <list> -t <tree> -ncat <ncat> <chainname> \n";
            cerr << '\n';
            exit(1);
        }

        MultiGeneAAMutSelDSBDPOmegaChain* chain = new MultiGeneAAMutSelDSBDPOmegaChain(datafile,treefile,Ncat,baseNcat,blmode,nucmode,basemode,omegamode,every,until,name,force,myid,nprocs);
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


