#pragma once

#include <string>
#include <vector>
#include "Array.hpp"
#include "SequenceAlignment.hpp"

class MultiGeneMPIModule {
  public:
    MultiGeneMPIModule(int inmyid, int innprocs) : myid(inmyid), nprocs(innprocs) {}
    ~MultiGeneMPIModule() {}

    int GetMyid() const { return myid; }

    int GetNprocs() const { return nprocs; }

    int GetNgene() const { return Ngene; }

    int GetLocalNgene() const { return LocalNgene; }

    int GetLocalTotNsite() const { return SlaveTotNsite[myid]; }

    int GetTotNsite() const { return SlaveTotNsite[0]; }

    int GetSlaveNgene(int proc) const {
        if (myid) {
            std::cerr << "error: slave in GetSlaveNgene\n";
            exit(1);
        }
        return SlaveNgene[proc];
    }

    std::string GetLocalGeneName(int gene) const { return GeneName[gene]; }

    int GetLocalGeneNsite(int gene) const { return GeneNsite[gene]; }

    int GetSlaveTotNsite(int proc) const { return SlaveTotNsite[proc]; }

    void AllocateAlignments(std::string datafile);

    void PrintGeneList(std::ostream &os) const;

    SequenceAlignment *refdata;  // FIXME: should not be public

  protected:
    int myid;
    int nprocs;

    int Ngene;
    int LocalNgene;
    std::vector<int> SlaveNgene;
    std::vector<int> SlaveTotNsite;
    std::vector<int> GeneAlloc;
    std::vector<std::string> GeneName;
    std::vector<int> GeneNsite;
};
