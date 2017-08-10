
#include "MultiGeneMPIModule.hpp"

#include <fstream>

void MultiGeneMPIModule::AllocateAlignments(string datafile)	{

	ifstream is(datafile.c_str());
	is >> Ngene;
    vector<string> genename(Ngene,"NoName");
	// genename.assign(Ngene,"NoName");
    vector<int> genesize(Ngene,0);
	// genesize.assign(Ngene,0);
    vector<int> genealloc(Ngene,0);
	// genealloc.assign(Ngene,0);
	vector<int> geneweight(Ngene,0);

	for (int gene=0; gene<Ngene; gene++)	{
		is >> genename[gene];
		SequenceAlignment* tmpdata = new FileSequenceAlignment(genename[gene]);

        if (! gene) {
            refdata = tmpdata;
        }

		genesize[gene] = tmpdata->GetNsite();
		geneweight[gene] = tmpdata->GetNsite() * tmpdata->GetNtaxa();

        if (gene)   {
            delete tmpdata;
        }
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

    SlaveNgene.assign(nprocs,0);
    for (int gene=0; gene<Ngene; gene++)	{
        SlaveNgene[genealloc[gene]]++;
    }
    
    if (! myid) {
        LocalNgene = Ngene;
        GeneName = genename;
        GeneAlloc = genealloc;
    }
    else    {
        GeneAlloc.assign(0,0);
        LocalNgene = SlaveNgene[myid];
        GeneName.assign(LocalNgene,"NoName");
        int i=0;
        for (int gene=0; gene<Ngene; gene++)	{
            if (genealloc[gene] == myid)    {
                GeneName[i] = genename[gene];
                i++;
            }
        }
    }
}
