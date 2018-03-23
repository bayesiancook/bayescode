#ifndef MULTISAMPLE_H
#define MULTISAMPLE_H

#include "Sample.hpp"
#include "MultiGeneProbModel.hpp"

/**
 * \brief A generic interface for a multi-gene MCMC sample
 */

class MultiGeneSample: public Sample    {

    public:

    //! \brief Constructor
	MultiGeneSample(string filename, int burnin, int every, int until, int inmyid, int innprocs) : Sample(filename,burnin,every,until), myid(inmyid), nprocs(innprocs) {}

    //! \brief non-const access to current configuration of underlying multi-gene model
    MultiGeneProbModel *GetMultiGeneModel() { return static_cast<MultiGeneProbModel*>(model); }
    //! \brief const access to current configuration of underlying multi-gene model
    const MultiGeneProbModel *GetMultiGeneModel() const { return static_cast<const MultiGeneProbModel*>(model); }

    //! \brief return number of genes
    int GetNgene() const {
        return GetMultiGeneModel()->GetNgene();
    }

    void OpenChainFile() override;
    void GetNextPoint() override;
    void SlaveRead();

    virtual void PostPred() override;

    protected:

    int myid;
    int nprocs;
};

#endif

