
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
