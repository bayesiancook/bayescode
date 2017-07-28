#ifndef MULTIGENE_H
#define MULTIGENE_H

class MultiGeneMPIModule    {

    public:

    MultiGeneMPIModule(int inmyid, int innprocs) : myid(inyid), nprocs(innprocs) {}
    ~MultiGeneMPIModule() {}

	int GetMyid() {
		return myid;
	}

	int GetNprocs() {
		return nprocs;
	}

    int GetNgene()  {
        return Ngene;
    }

    void AllocateAlignments(string datafile);

    protected:

    int myid;
    int nprocs;

	int Ngene;
    std::vector<int> genealloc;
    std::vector<int> genesize;
    std::vector<string> genename;
};

#endif
