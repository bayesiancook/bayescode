#ifndef MULTIGENE_H
#define MULTIGENE_H

#include <vector>
#include <string>
using namespace std;
#include "SequenceAlignment.hpp"

class MultiGeneMPIModule    {

    public:

    MultiGeneMPIModule(int inmyid, int innprocs) : myid(inmyid), nprocs(innprocs) {}
    ~MultiGeneMPIModule() {}

	int GetMyid() const {
		return myid;
	}

	int GetNprocs() const {
		return nprocs;
	}

    int GetNgene() const {
        return Ngene;
    }

    int GetLocalNgene() const   {
        return LocalNgene;
    }

    int GetLocalTotNsite() const    {
        return SlaveTotNsite[myid];
    }

    int GetSlaveNgene(int proc) const   {
        if (myid)   {
            cerr << "error: slave in GetSlaveNgene\n";
            exit(1);
        }
        return SlaveNgene[proc];
    }

    string GetLocalGeneName(int gene) const {
        return GeneName[gene];
    }

    int GetLocalGeneNsite(int gene) const   {
        return GeneNsite[gene];
    }

    int GetSlaveTotNsite(int proc) const    {
        return SlaveTotNsite[proc];
    }

    void AllocateAlignments(string datafile);

    void PrintGeneList(ostream& os) const;

    protected:

    int myid;
    int nprocs;

	int Ngene;
    int LocalNgene;
    std::vector<int> SlaveNgene;
    std::vector<int> SlaveTotNsite;
    std::vector<int> GeneAlloc;
    std::vector<string> GeneName;
    std::vector<int> GeneNsite;

    SequenceAlignment* refdata;
};

#endif
