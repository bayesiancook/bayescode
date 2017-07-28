
#include "MultiGeneMPIModule.hpp"

#include <fstream>
#include "SequenceAlignment.hpp"

void MultiGeneMPIModule::AllocateAlignments(string datafile)	{

	ifstream is(datafile.c_str());
	is >> Ngene;
	genename.assign(Ngene,"NoName");
	genesize.assign(Ngene,0);
	genealloc.assign(Ngene,0);
	vector<int> geneweight(Ngene,0);

	for (int gene=0; gene<Ngene; gene++)	{
		is >> genename[gene];
		SequenceAlignment* tmpdata = new FileSequenceAlignment(genename[gene]);
        /*
		if (! gene)	{
			data = tmpdata;
		}
		else	{
			if (nstate != data->GetNstate())	{
				cerr << "error: all data files do not have the same alphabet\n";
				cerr << nstate << '\t' << data->GetNstate() << '\n';
				exit(1);
			}
		}
        */

		genesize[gene] = tmpdata->GetNsite();
		geneweight[gene] = tmpdata->GetNsite() * tmpdata->GetNtaxa();
        delete tmpdata;
        /*
		if (gene)	{
			delete tmpdata;
		}
        */
	}

	// sort alignments by decreasing size
    std::vector<int> permut(Ngene);
	for (int gene=0; gene<Ngene; gene++)	{
		permut[gene] = gene;
	}
	for (int i=0; i<Ngene; i++)	{
		for (int j=Ngene-1; j>i; j--)	{
			if (geneweight[permut[i]] < geneweight[permut[j]])	{
			// if (genesize[permut[i]] < genesize[permut[j]])	{
				int tmp = permut[i];
				permut[i] = permut[j];
				permut[j] = tmp;
			}
		}
	}

	int totsize[nprocs];
	for (int i=0; i<nprocs; i++)	{
		totsize[i] = 0;
	}

	for (int i=0; i<Ngene; i++)	{
		int gene = permut[i];
		int size = geneweight[gene];
		// int size = genesize[gene];

		int min = 0;
		int jmin = 0;
		for (int j=1; j<nprocs; j++)	{
			if ((j==1) || (min > totsize[j]))	{
				min = totsize[j];
				jmin = j;
			}
		}
		genealloc[gene] = jmin;
		totsize[jmin] += size;
	}

	if (totsize[0])	{
		cerr << "error in alloc\n";
		exit(1);
	}
	int total = 0;
	for (int i=1; i<nprocs; i++)	{
		int tot = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == i)	{
				tot += geneweight[gene];
				// tot += genesize[gene];
				total++;
			}
		}
		if (tot != totsize[i])	{
			cerr << "error in allocation\n";
			cerr << tot << '\t' << totsize[i] << '\n';
			exit(1);
		}
	}
	if (total != Ngene)	{
		cerr << "error in total allocation\n";
		exit(1);
	}

	for (int gene=0; gene<Ngene; gene++)	{
		if ((genealloc[gene] < 0) || (genealloc[gene] >= nprocs))	{
			cerr << "alloc : " << genealloc[gene] << '\t' << gene << '\n';
			exit(1);
		}
	}
}
