#ifndef MULTISAMPLE_H
#define MULTISAMPLE_H

#include "Sample.hpp"
#include "MultiGeneProbModel.hpp"

class MultiGeneSample: public Sample    {

    public:

	MultiGeneSample(string filename, int burnin, int every, int until, int inmyid, int innprocs) : Sample(filename,burnin,every,until), myid(inmyid), nprocs(innprocs) {}

    MultiGeneProbModel *GetMultiGeneModel() { return static_cast<MultiGeneProbModel*>(model); }
    const MultiGeneProbModel *GetMultiGeneModel() const { return static_cast<const MultiGeneProbModel*>(model); }

    int GetNgene() const {
        return GetMultiGeneModel()->GetNgene();
    }

    void OpenChainFile() override;
    void GetNextPoint() override;

    protected:

    int myid;
    int nprocs;
};

#endif

