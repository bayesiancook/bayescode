#ifndef CHAIN_H
#define CHAIN_H

#include <string>
#include "ProbModel.hpp"

/**
 * \brief A generic interface for a Monte Carlo Markov Chain
 *
 * Chain is responsible for creating a model,
 * by calling the model constructor with the relevant settings,
 * and then running the MCMC, regularly saving samples into a file.
 * A chain can be stopped, and then restarted from file.
 *
 * The files have the following extensions:
 * - <chainname>.param   : current state (detailed settings and complete model
 * configuration)
 * - <chainname>.chain   : list of all saved points since the beginning of the
 * MCMC (burnin included)
 * - <chainname>.trace   : trace file, each row corresponding to one point of
 * the points saved during the MCMC
 * - <chainname>.monitor : monitoring the success rate, time spent in each move,
 * numerical errors, etc
 * - <chainname>.run     : put 0 in this file to stop the chain
 */

class Chain {
  public:
    Chain();

    virtual ~Chain() = default;

    //! \brief return model type
    //!
    //! each derived class should define a unique string for each type of model
    //! (typically used to check that a chain restarted from file is from correct
    //! model), and then override this pure virtual function to return the type.
    virtual std::string GetModelType() = 0;

    //! make new chain (force == 1 : overwrite already existing files with same
    //! name)
    virtual void New(int force = 0) = 0;

    //! open chain from file
    virtual void Open() = 0;

    //! save chain to file
    virtual void Save() = 0;

    //! initialise model and make the files (force == 1: overwrite existing files
    //! with same name)
    virtual void Reset(int force = 0);

    //! create all files (chain, trace, monitor; called when creating a new chain)
    virtual void MakeFiles(int force = 0);

    //! start the MCMC
    virtual void Start();

    //! run the MCMC: cycle over Move, Monitor and Save while running status == 1
    virtual void Run();

    //! perform one cycle of Monte Carlo "moves" (updates)
    virtual void Move();

    //! save one point in the .chain file (called after each cycle)
    virtual void SavePoint();

    //! write current trace and monitoring statistics in the .trace and .monitor
    //! files (called after each cycle)
    virtual void Monitor();

    //! \brief returns running status (1: run should continue / 0: run should now
    //! stop)
    //!
    //! returns 1 (means continue) if one the following conditions holds true
    //! - <chainname>.run file contains a 1
    //! - size < until, or until == -1
    //!
    //! Thus, "echo 0 > <chainname>.run" is the proper way to stop a chain from a
    //! shell
    virtual int GetRunningStatus();

    //! return chain name: i.e. base name for all files corresponding to that
    //! chain
    std::string GetName() { return name; }

    //! return pointer to underlying model
    ProbModel *GetModel() { return model; }

    //! return current size (number of points saved to file thus far)
    int GetSize() { return size; }

  protected:
    //! saving frequency (i.e. number of move cycles performed between each point
    //! saved to file)
    int every;
    //! intended final size of the chain (until==-1 means no a priori specified
    //! upper limit)
    int until;
    //! current size (number of points saved to file)
    int size;
    //! pointer to the underlying model
    ProbModel *model;
    //! base name for all files corresponding to that chain
    std::string name;
    //! flag: if 1, then complete state is saved at each interation in .chain file
    int saveall;
};

#endif  // CHAIN_H
