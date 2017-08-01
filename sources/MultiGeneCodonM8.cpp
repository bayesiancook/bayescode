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

	string datafile = argv[1];
	string treefile = argv[2];
	int ncat = atoi(argv[3]);
	int withpos = atoi(argv[4]);
	string name = argv[5];

	MultiGeneCodonM8Model* model = new MultiGeneCodonM8Model(datafile,treefile,ncat,withpos,myid,nprocs);
    if (! myid) {
        cerr << " -- master unfold\n";
    }
    if (! myid) {
        cerr << " -- start\n";
    }
    model->Unfold();
    if (! myid) {
        ofstream pos((name + ".posw").c_str());
        ofstream os((name + ".trace").c_str());
        model->TraceHeader(os);
        os.flush();
        while(1)	{
            model->MasterMove();
            model->MasterTrace(os);
            model->TracePosWeight(pos);
        }
    }
    else	{
        while(1)	{
            model->SlaveMove();
            model->SlaveTrace();
        }
    }

	MPI_Finalize();
}


