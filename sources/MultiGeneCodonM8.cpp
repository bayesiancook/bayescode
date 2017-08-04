#include <cmath>
#include <fstream>
#include "MultiGeneCodonM8Model.hpp"

using namespace std;

MPI_Datatype Propagate_arg;

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

	string datafile = "";
	string treefile = "";
	int ncat = 1;
    double pihypermean = 0.1;
    double pihyperinvconc = 0.2;

    int writegenedata = 0;

	string name = "";

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
            /*
            else if (s == "-f")	{
                force = 1;
            }
            */
            else if (s == "-ncat")  {
                i++;
                ncat = atoi(argv[i]);
            }
            else if (s == "-pi")    {
                i++;
                pihypermean = atof(argv[i]);
                i++;
                pihyperinvconc = atof(argv[i]);
            }
            else if (s == "-g")  {
                writegenedata = 1;
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
        cerr << "codonm8 -d <alignment> -t <tree> [-fixparam <paramfile> -fixhyper <hyperparamfile>] <chainname> \n";
        cerr << '\n';
        exit(1);
    }

	MultiGeneCodonM8Model* model = new MultiGeneCodonM8Model(datafile,treefile,ncat,pihypermean,pihyperinvconc,myid,nprocs);
    if (! myid) {
        cerr << " -- master unfold\n";
    }
    if (! myid) {
        cerr << " -- start\n";
    }
    model->Allocate();
    model->Unfold();
    if (! myid) {
        ofstream paramos((name + ".globalparam").c_str());
        ofstream hyperos((name + ".hyperparam").c_str());
        ofstream pos((name + ".posw").c_str());
        ofstream omos((name + ".posom").c_str());
        ofstream os((name + ".trace").c_str());
        model->TraceHeader(os);
        os.flush();
        while(1)	{
            model->MasterMove();
            model->MasterTraceGlobalParameters(paramos);
            model->MasterTraceHyperParameters(hyperos);
            model->MasterTrace(os);
            model->TracePosWeight(pos);
            model->TracePosOm(omos);
        }
    }
    else	{
        if (writegenedata)  {
            model->SlaveTracePostProbHeader(name);
        }
        while(1)	{
            model->SlaveMove();
            model->SlaveTrace();
            if (writegenedata)  {
                model->SlaveTracePostProb(name);
            }
        }
    }

	MPI_Finalize();
}


