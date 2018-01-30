#ifndef MULTICHAIN_H
#define MULTICHAIN_H

#include <string>
#include "Chain.hpp"
#include "MultiGeneProbModel.hpp"

class MultiGeneChain : public Chain {
  public:

    MultiGeneChain(int myid, int nprocs);

    virtual ~MultiGeneChain() = default;

    virtual void SavePoint() override;

    virtual void Reset(int force) override;

    virtual void Move() override;
    // perform one cycle of Monte Carlo "moves" (updates)

    virtual void Start() override;
    // start Monte Carlo

    void MasterSendRunningStatus(int status);
    int SlaveReceiveRunningStatus();

    virtual void Run() override;
    // Move, Monitor amd Save while running status == 1

    MultiGeneProbModel *GetMultiGeneModel() { return static_cast<MultiGeneProbModel*>(model); }

  protected:
    int myid;
    int nprocs;

};

#endif  // MULTICHAIN_H
