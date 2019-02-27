#ifndef PROBMODEL_H
#define PROBMODEL_H

#include <iostream>
#include <set>
#include "global/Random.hpp"

/**
 * \brief A generic interface for MCMC probabilistic models
 *
 * ProbModel provides an interface for the fundamental methods that any model
 * should implement in BayesCode:
 * - Move: make a complete cycle of MCMC moves
 * - GetLogProb: return the log prob of the current model configuration
 * - Update: update the entire model (e.g. after reading an instance from file)
 * - FromStream / ToStream: intput output of entire model configuration
 * - Trace / TraceHeader / Monitor: tracing and monitoring the model
 * Currently, ProbModel does not declare those methods pure virtual (although
 * this could perhaps be enforced)
 *
 * In addition, ProbModel proposes three templates for sliding, scaling and
 * profile moves on model parameters.
 *
 */

using namespace std;

class ProbModel {
  public:
    ProbModel() {}
    ~ProbModel() {}

    //! make a complete cycle of MCMC moves -- in principle, should return average
    //! success rate (although rarely does so in practice)
    virtual double Move() { return 1; }

    //! update the entire model
    virtual void Update() {
        std::cerr << "error : in ProbModel::Update\n";
        exit(1);
    }

    //! post pred method
    virtual void PostPred(std::string name) {
        std::cerr << "error: in ProbModel::PostPred\n";
        exit(1);
    }

    //! return lof prob of the current model configuration
    virtual double GetLogProb() const { return 0; }

    //! save model configuration to stream
    virtual void ToStream(std::ostream &os) const {}
    //! get model configuration from stream
    virtual void FromStream(std::istream &is) {}

    //! write one line of trace of the current parameter configuration into trace
    //! file
    virtual void Trace(std::ostream & /*unused*/) const {}
    //! write one line of header for the trace file
    virtual void TraceHeader(std::ostream & /*unused*/) const {}
    //! output statistics monitoring the MCMC
    virtual void Monitor(std::ostream &os) const {}

    //! new type name for a const method of class C taking no argument and
    //! returning a double (intended: a log prob function)
    template <class C>
    using LogProbF = double (C::*)(void) const;
    //! new type name for a non-const method of class C taking no argument and
    //! with no return (intended: an update function)
    template <class C>
    using UpdateF = void (C::*)(void);

    //! \brief template for Metropolis Hastings sliding move on parameter x
    //!
    //! Parameterized by a tuning parameter, un number of iterations (nrep),
    //! and lower and upper constraint (min and max; if max<min, then no
    //! constraint is enforced). Should also give a pointer to a log prob and an
    //! update functions, as well as a pointer to the model itself. Returns
    //! success rate.
    template <class C>
    double SlidingMove(double &x, double tuning, int nrep, double min, double max,
        LogProbF<C> logprobf, UpdateF<C> updatef, C *This) {
        // C* This = dynamic_cast<C*>(this);

        double nacc = 0;
        double ntot = 0;
        for (int rep = 0; rep < nrep; rep++) {
            double bk = x;
            double deltalogprob = -(This->*logprobf)();
            double m = tuning * (Random::Uniform() - 0.5);
            x += m;
            if (max > min) {
                while ((x < min) || (x > max)) {
                    if (x < min) { x = 2 * min - x; }
                    if (x > max) { x = 2 * max - x; }
                }
            }
            (This->*updatef)();
            deltalogprob += (This->*logprobf)();
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                x = bk;
                (This->*updatef)();
            }
            ntot++;
        }
        return nacc / ntot;
    }

    //! \brief template for Metropolis Hastings scaling move on parameter x
    //!
    //! Parameterized by a tuning parameter, un number of iterations (nrep),
    //! and lower and upper constraint (min and max; if max<min, then no
    //! constraint is enforced). Should also give a pointer to a log prob and an
    //! update functions, as well as a pointer to the model itself. Returns
    //! success rate.
    template <class C>
    double ScalingMove(
        double &x, double tuning, int nrep, LogProbF<C> logprobf, UpdateF<C> updatef, C *This) {
        double nacc = 0;
        double ntot = 0;
        for (int rep = 0; rep < nrep; rep++) {
            double deltalogprob = -(This->*logprobf)();
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            x *= e;
            (This->*updatef)();
            deltalogprob += (This->*logprobf)();
            deltalogprob += m;
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                x /= e;
                (This->*updatef)();
            }
            ntot++;
        }
        return nacc / ntot;
    }

    //! \brief template for Metropolis Hastings move on a frequency vector
    //! (profile) x of the model
    //!
    //! Parameterized by a real and an integer tuning parameters (tuning and n), a
    //! number of iterations (nrep), and lower and upper constraint (min and max;
    //! if max<min, then no constraint is enforced). Should also give a pointer to
    //! a log prob and an update functions, as well as a pointer to the model
    //! itself. Returns success rate.
    template <class C>
    double ProfileMove(std::vector<double> &x, double tuning, int n, int nrep, LogProbF<C> logprobf,
        UpdateF<C> updatef, C *This) {
        double nacc = 0;
        double ntot = 0;
        std::vector<double> bk(x.size(), 0);
        for (int rep = 0; rep < nrep; rep++) {
            bk = x;
            double deltalogprob = -(This->*logprobf)();
            double loghastings = Random::ProfileProposeMove(x, x.size(), tuning, n);
            (This->*updatef)();
            deltalogprob += (This->*logprobf)();
            deltalogprob += loghastings;
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                x = bk;
                (This->*updatef)();
            }
            ntot++;
        }
        return nacc / ntot;
    }
};

#endif  // PROBMODEL_H
