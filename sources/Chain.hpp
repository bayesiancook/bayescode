#ifndef CHAIN_H
#define CHAIN_H

#include "component_defs.hpp"
class ProbModel;
/**
 * \brief A generic interface for a Monte Carlo Markov ChainDriver
 *
 * ChainDriver is responsible for creating a model,
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

class ChainDriver : public Go {
    ProbModel* model{nullptr};

    Lifecycle* lifecycle_handler;

    AbstractTraceFile* chainfile{nullptr};
    AbstractTraceFile* monitorfile{nullptr};
    AbstractTraceFile* paramfile{nullptr};
    AbstractTraceFile* tracefile{nullptr};

    RunToggle* run_toggle{nullptr};

  public:
    ChainDriver(int every, int until) : every(every), until(until) {
        port("lifecycle", &ChainDriver::lifecycle_handler);
        port("model", &ChainDriver::model);
        port("chainfile", &ChainDriver::chainfile);
        port("monitorfile", &ChainDriver::monitorfile);
        port("paramfile", &ChainDriver::paramfile);
        port("tracefile", &ChainDriver::tracefile);
        port("runtoggle", &ChainDriver::run_toggle);
    }

    //! start the MCMC
    void go() override;

    //! run the MCMC: cycle over Move, Monitor and Save while running status == 1
    void Run();

    //! perform one cycle of Monte Carlo "moves" (updates)
    void Move();

    //! \brief returns running status (1: run should continue / 0: run should now stop)
    //!
    //! returns 1 (means continue) if one the following conditions holds true
    //! - <chainname>.run file contains a 1
    //! - size < until, or until == -1
    //!
    //! Thus, "echo 0 > <chainname>.run" is the proper way to stop a chain from a shell
    bool IsRunning();

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
