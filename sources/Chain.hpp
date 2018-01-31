#ifndef CHAIN_H
#define CHAIN_H

#include <string>
#include "ProbModel.hpp"
#include "component_defs.hpp"

/**
 * \brief A generic interface for a Monte Carlo Markov Chain
 *
 * Chain is responsible for creating a model,
 * by calling the model constructor with the relevant settings,
 * and then running the MCMC, regularly saving samples into a file.
 * A chain can be stopped, and then restarted from file.
 *
 * The files have the following extensions:
 * - <chainname>.param   : current state (detailed settings and complete model configuration)
 * - <chainname>.chain   : list of all saved points since the beginning of the MCMC (burnin included)
 * - <chainname>.trace   : trace file, each row corresponding to one point of the points saved during the MCMC
 * - <chainname>.monitor : monitoring the success rate, time spent in each move, numerical errors, etc
 * - <chainname>.run     : put 0 in this file to stop the chain
 */

class Chain : public Go {
    ProbModel* model{nullptr};

    TraceFile* chainfile{nullptr};
    TraceFile* fitnessfile{nullptr};
    TraceFile* monitorfile{nullptr};
    TraceFile* paramfile{nullptr};
    TraceFile* tracefile{nullptr};

    RunToggle* run_toggle{nullptr};

  public:
    Chain() {
        port("model", &Chain::model);
        port("chainfile", &Chain::chainfile);
        port("fitnessfile", &Chain::fitnessfile);
        port("monitorfile", &Chain::monitorfile);
        port("paramfile", &Chain::paramfile);
        port("tracefile", &Chain::tracefile);
        port("runtoggle", &Chain::run_toggle);
    }

    ~Chain() = default;

    // //! make new chain (force == 1 : overwrite already existing files with same name)
    // void New(bool force = false) = 0;

    // //! open chain from file
    // void Open() = 0;

    // //! save chain to file
    // void Save() = 0;

    //! initialise model and make the files (force == 1: overwrite existing files with same name)
    void Reset(bool force = false);

    // //! create all files (chain, trace, monitor; called when creating a new chain)
    void MakeFiles(bool force = false);

    //! start the MCMC
    void go() override {}  // TODO

    //! start the MCMC
    void Start();

    //! run the MCMC: cycle over Move, Monitor and Save while running status == 1
    void Run();

    //! perform one cycle of Monte Carlo "moves" (updates)
    void Move();

    //! save one point in the .chain file (called after each cycle)
    void SavePoint();

    //! write current trace and monitoring statistics in the .trace and .monitor files (called after each cycle)
    void Monitor();

    //! \brief returns running status (1: run should continue / 0: run should now stop)
    //!
    //! returns 1 (means continue) if one the following conditions holds true
    //! - <chainname>.run file contains a 1
    //! - size < until, or until == -1
    //!
    //! Thus, "echo 0 > <chainname>.run" is the proper way to stop a chain from a shell
    bool IsRunning();

    //! return pointer to underlying model
    ProbModel* GetModel() { return model; }

    //! return current size (number of points saved to file thus far)
    int GetSize() { return size; }

  protected:
    //! saving frequency (i.e. number of move cycles performed between each point saved to file)
    int every{1};
    //! intended final size of the chain (until==-1 means no a priori specified upper limit)
    int until{-1};
    //! current size (number of points saved to file)
    int size{0};
};

#endif  // CHAIN_H
