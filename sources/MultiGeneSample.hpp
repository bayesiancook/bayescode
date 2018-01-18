#ifndef MULTISAMPLE_H
#define MULTISAMPLE_H

#include "Sample.hpp"
#include "MultiGeneProbModel.hpp"

class MultiGeneSample: public Sample    {

    public:

	MultiGeneSample(string filename, int burnin, int every, int until, int inmyid, int innprocs) : Sample(filename,burnin,every,until), myid(inmyid), nprocs(innprocs) {}
    MultiGeneProbModel *GetMultiGeneModel() { return static_cast<MultiGeneProbModel*>(model); }

    protected:

    int myid;
    int nprocs;
};

#endif

