#ifndef MULTICHAIN_H
#define MULTICHAIN_H

#include <string>
#include "Chain.hpp"
#include "MultiGeneProbModel.hpp"
#include "Chrono.hpp"

/**
 * \brief A generic interface for a multi-gene Monte Carlo Markov Chain
 *
 * MultiGeneChain overrides several functions of Chain and defines new methods,
 * implementing part of the behavior of the class relating to the multi-gene and
 * multi-process settings.
 */

class MultiGeneChain : public Chain {
  public:
    //! constructor parameterized by process id and number of processes (as in
    //! MultiGeneMPIModule)
    MultiGeneChain(int myid, int nprocs);

    virtual ~MultiGeneChain() throw() = default;

    virtual void SavePoint() override;

    virtual void Reset(int force) override;

    virtual void Move() override;

    virtual void Start() override;

    virtual void MakeFiles(int force) override;

    virtual int GetRunningStatus() override;

    //! \brief master broadcasts running status (1: run should continue / 0: run
    //! should now stop)
    void MasterSendRunningStatus(int status);
    //! \brief slave receive running status from master (and then stop if 0)
    int SlaveReceiveRunningStatus();

    virtual void Run() override;

    //! \brief returns pointer to multi-gene model (static cast of ProbModel*
    //! model of Chain)
    MultiGeneProbModel *GetMultiGeneModel() { return static_cast<MultiGeneProbModel *>(model); }

    void SetMaxTime(double inmaxtime)   {
        maxtime = inmaxtime;
    }

  protected:
    int myid;
    int nprocs;

    double maxtime;
    Chrono global_chrono;
};

#endif  // MULTICHAIN_H
