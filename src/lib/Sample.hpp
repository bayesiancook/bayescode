#ifndef SAMPLE_H
#define SAMPLE_H

#include <cstdlib>
#include <fstream>

#include "ProbModel.hpp"

/**
 * \brief A generic interface for an MCMC sample
 *
 * In terms of external behavior and itended use,
 * Sample is in charge of reading chains from file,
 * and computing posterior quantities (averages, credible intervals,
 * distributions, posterior predictive samples. etc). As a simple example, see
 * SingleOmegaSample.
 */

class Sample {
  public:
    //! \brief Constructor, opening chain from file, with specified burn-in,
    //! thinning factor and upper limit.
    Sample(std::string filename, int in_burnin = 0, int in_every = 1, int in_until = -1);

    virtual ~Sample();

    //! \brief return the base name of all files (same as for the Chain class)
    std::string GetName() { return name; }

    //! \brief get the next point (automatically accounts for subsampling, as
    //! specified by the thinning parameter)
    //!
    //! the parameter configuration specified by this point of the MCMC is then
    //! accessible through GetModel()
    virtual void GetNextPoint();

    virtual void PostPred();

    //! \brief return a pointer to model configuration specified by current point
    //! (i.e. last point that was read from file)
    //!
    //! in general this function will be overriden by derived classes,
    //! since derived classes manipulate models that are derived from ProbModel
    virtual ProbModel *GetModel() { return model; }

    //! \brief return model type (see Chain)
    virtual std::string GetModelType() = 0;

    //! \brief create MCMC sample
    virtual void Open() = 0;

    //! \brief opens files, prepare data structures and sets the stream iterator
    virtual void OpenChainFile();

    int size;  // sample size (calculated from parameters above)

  protected:
    std::ifstream *chain_is;
    int chainevery;  // chain's saving frequency
    int chainuntil;  // chain's intended size of the run (number of saved points)
    int chainsize;   // chain's current size
    //! flag: if 1, then complete state is saved at each interation in .chain file
    int chainsaveall;
    int burnin;  // burnin
    int every;   // subsampling frequency
    int until;   // reading chain until this point
    int currentpoint;
    ProbModel *model;  // the model
    std::string name;  // the name of the chain in the filesystem
};

#endif  // SAMPLE_H
