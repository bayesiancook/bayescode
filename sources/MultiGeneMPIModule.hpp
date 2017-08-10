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

    string GetLocalGeneName(int gene) const {
        return GeneName[gene];
    }

    void AllocateAlignments(string datafile);

    protected:

    int myid;
    int nprocs;

	int Ngene;
    int LocalNgene;
    std::vector<int> genealloc;
    // std::vector<int> genesize;
    std::vector<string> GeneName;

    SequenceAlignment* refdata;
};

#endif
