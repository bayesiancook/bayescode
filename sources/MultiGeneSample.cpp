
#include "MultiGeneSample.hpp"

void MultiGeneSample::OpenChainFile()	{

	if (until == -1)	{
		until = chainsize;
	}
	if (until > chainsize)	{
		cerr << "number of points saved is less than " << until << '\n';
		until = chainsize;
	}
	if (burnin == -1)	{
		burnin = chainsize / 10;
	}
	size = (until-burnin) / every;
	if (size <= 0)	{
		cerr << "error : chain not long enough\n";
		exit(1);
	}
	currentpoint = 0;

    if (! myid) {
        chain_is = new ifstream((name + ".chain").c_str());
        if (! chain_is)	{
            cerr << "error: cannot find file " << name << ".chain\n";
            exit(1);
        }
        for (int i=0; i<burnin; i++)	{
            GetMultiGeneModel()->MasterFromStream(*chain_is);
        }
    }
    else    {
        for (int i=0; i<burnin; i++)	{
            GetMultiGeneModel()->SlaveFromStream();
        }
    }
}

void MultiGeneSample::GetNextPoint()	{

	if (currentpoint == size)	{
		cerr << "error in Sample::GetNextPoint: going past last points\n";
		exit(1);
	}
	if (currentpoint)	{
        if (!myid)  {
            for (int i=0; i<every-1; i++)	{
                GetMultiGeneModel()->MasterFromStream(*chain_is);
            }
        }
        else    {
            for (int i=0; i<every-1; i++)	{
                GetMultiGeneModel()->SlaveFromStream();
            }
        }
	}
	currentpoint++;
}

void MultiGeneSample::MasterPostPred() {

    cerr << size << " points to read\n";
    for (int i=0; i<size; i++)  {
        cerr << '.';
        GetNextPoint();
        ostringstream s;
        s << name << "_" << i << ".ali";
        model->PostPred(s.str());
        MPI_Barrier(MPI_COMM_WORLD);
    }
    cerr << '\n';
}

void MultiGeneSample::SlavePostPred() {

    for (int i=0; i<size; i++)  {
        GetNextPoint();
        ostringstream s;
        s << name << "_" << i << ".ali";
        model->PostPred(s.str());
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

