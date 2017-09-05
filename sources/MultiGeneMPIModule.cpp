
#include "MultiGeneMPIModule.hpp"
// #include "Parallel.hpp"

#include <fstream>

void MultiGeneMPIModule::AllocateAlignments(string datafile)	{

	ifstream is(datafile.c_str());
	is >> Ngene;
    vector<string> genename(Ngene,"NoName");
    vector<int> genesize(Ngene,0);
    vector<int> genealloc(Ngene,0);
	vector<int> geneweight(Ngene,0);

	for (int gene=0; gene<Ngene; gene++)	{
		is >> genename[gene];
		SequenceAlignment* tmpdata = new FileSequenceAlignment(genename[gene]);

        if (! gene) {
            refdata = tmpdata;
        }

		genesize[gene] = tmpdata->GetNsite() / 3;
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
    SlaveTotNsite.assign(nprocs,0);
    for (int gene=0; gene<Ngene; gene++)	{
        if (! genealloc[gene])  {
            cerr << "error: gene allocated to master\n";
            exit(1);
        }
        SlaveNgene[genealloc[gene]]++;
        SlaveTotNsite[genealloc[gene]] += genesize[gene];
        SlaveTotNsite[0] += genesize[gene];
    }
    
    if (! myid) {
        LocalNgene = Ngene;
        GeneName.assign(Ngene,"noname");
        GeneNsite.assign(Ngene,0);
        GeneAlloc.assign(Ngene,0);
        int i = 0;
        cerr << '\n';
        cerr << "proc\tngene\ttotnsite\n";
        for (int proc=1; proc<nprocs; proc++)   {
            cerr << proc << '\t' << SlaveNgene[proc] << '\t' << SlaveTotNsite[proc] << '\n';
            for (int gene=0; gene<Ngene; gene++)    {
                if (genealloc[gene] == proc)    {
                    GeneAlloc[i] = proc;
                    GeneName[i] = genename[gene];
                    GeneNsite[i] = genesize[gene];
                    i++;
                }
            }
        }
        cerr << '\n';
        if (i != Ngene) {
            cerr << "error in mpimodule: non matching number of genes\n";
        }
    }
    else    {
        GeneAlloc.assign(0,0);
        LocalNgene = SlaveNgene[myid];
        GeneName.assign(LocalNgene,"NoName");
        GeneNsite.assign(LocalNgene,0);
        int i=0;
        for (int gene=0; gene<Ngene; gene++)	{
            if (genealloc[gene] == myid)    {
                GeneName[i] = genename[gene];
                GeneNsite[i] = genesize[gene];
                i++;
            }
        }
    }
}

void MultiGeneMPIModule::PrintGeneList(ostream& os) const {

    if (myid)   {
        cerr << "error: slave in MultiGeneMPIModule::PrintGeneList\n";
        exit(1);
    }
    os << Ngene << '\n';
    for (int gene=0; gene<Ngene; gene++)    {
        os << GeneName[gene] << '\t' << GeneNsite[gene] << '\n';
    }
}


/*
void MultiGeneMPIModule::MasterSendGlobalVariable(const MPIObject& object) {

    MPIBuffer buffer(object.GetMPISize());
    buffer << object;
    // object.MPIPut(buffer);
    MPI_Bcast(buffer.GetBuffer(),buffer.GetBufferSize(),MPI_DOUBLE,0,MPI_COMM_WORLD);
}


void MultiGeneMPIModule::SlaveReceiveGlobalVariable(MPIObject& object) {

    MPIBuffer buffer(object.GetMPISize());
    MPI_Bcast(buffer.GetBuffer(),buffer.GetBufferSize(),MPI_DOUBLE,0,MPI_COMM_WORLD);
    // object.MPIGet(buffer);
    buffer >> object;
}

void MultiGeneMPIModule::SlaveSendAdditive(const MPIAdditive& object) {

    MPIBuffer buffer(object.GetMPISize());
    buffer << object;
    // object.MPIPut(buffer);
    MPI_Send(buffer.GetBuffer(),buffer.GetBufferSize(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void MultiGeneMPIModule::MasterReceiveAdditive(MPIAdditive& object)   {

    object.Clear();
    MPI_Status stat;
    for (int proc=1; proc<GetNprocs(); proc++)  {
        MPIBuffer buffer(object.GetMPISize());
        MPI_Recv(buffer.GetBuffer(),buffer.GetBufferSize(),MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);
        // object.MPIAdd(buffer);
        object += buffer;
    }
}

template<class T> unsigned int MPIElementSize<T> = -1;

template<class T> void MultiGeneMPIModule::MasterSendGeneArray(const SimpleArray<T>&) {

}

template<class T> void MultiGeneMPIModule::SlaveReceiveGeneArray(SimpleArray<T>&) {

}

void MultiGeneMPIModule::MasterSendGeneArray(const MPIArray& array)   {

    unsigned int elemsize = array.GetMPIElementSize();
    int thusfar = 0;
    for (int proc=1; proc<GetNprocs(); proc++)  {
        unsigned int buffersize = elemsize * GetSlaveNgene(proc);
        MPIBuffer buffer(buffersize);
        array.MPIPutArray(buffer,thusfar,GetSlaveNgene(proc));
        thusfar += GetSlaveNgene(proc);
        MPI_Send(buffer.GetBuffer(),buffer.GetBufferSize(),MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD);
    }
}

        
void MultiGeneMPIModule::SlaveReceiveGeneArray(MPIArray& array)  {

    unsigned int elemsize = array.GetMPIElementSize();
    unsigned int buffersize = elemsize * GetLocalNgene();
    MPIBuffer buffer(buffersize);
    MPI_Status stat;
    MPI_Recv(buffer.GetBuffer(),buffer.GetBufferSize(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD,&stat);
    array.MPIGetArray(buffer,0,GetLocalNgene());
}

void MultiGeneMPIModule::SlaveSendGeneArray(const MPIArray& array)  {

    unsigned int elemsize = array.GetMPIElementSize();
    unsigned int buffersize = elemsize * GetLocalNgene();
    MPIBuffer buffer(buffersize);
    array.MPIPutArray(buffer,0,GetLocalNgene());
    MPI_Send(buffer.GetBuffer(),buffer.GetBufferSize(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void MultiGeneMPIModule::MasterReceiveGeneArray(MPIArray& array)    {

    unsigned int elemsize = array.GetMPIElementSize();
    MPI_Status stat;
    int thusfar = 0;
    for (int proc=1; proc<GetNprocs(); proc++)  {
        unsigned int buffersize = elemsize * GetSlaveNgene(proc);
        MPIBuffer buffer(buffersize);
        MPI_Recv(buffer.GetBuffer(),buffer.GetBufferSize(),MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);
        array.MPIGetArray(buffer,thusfar,GetSlaveNgene(proc));
        thusfar += GetSlaveNgene(proc);
    }
}
*/
