#ifndef PROBMODEL_H
#define PROBMODEL_H

#include <iostream>
#include <set>
#include "Random.hpp"
using namespace std;

class ProbModel {
  public:
    ProbModel() {}
    ~ProbModel() {}

    virtual double Move() {return 1;}

    // void NoUpdate() {}
    virtual void Update() {}
    virtual double GetLogProb() const {return 0;}

    // save model configuration to stream
    virtual void ToStream(std::ostream &os) const {}
    // get model configuration from stream
    virtual void FromStream(std::istream &is) {}

    // monitoring the run
    virtual void Trace(std::ostream & /*unused*/) const {}
    virtual void TraceHeader(std::ostream & /*unused*/) const {}
    virtual void Monitor(std::ostream &os) const {}

    // templates for Metropolis Hastings Moves
    template<class C> using LogProbF = double (C::*)(void) const;
    template<class C> using UpdateF = void (C::*)(void);

	template<class C> double SlidingMove(double& x, double tuning, int nrep, double min, double max, LogProbF<C> logprobf, UpdateF<C> updatef, C* This) {
    
        // C* This = dynamic_cast<C*>(this);

        double nacc = 0;
        double ntot = 0;
        for (int rep=0; rep<nrep; rep++)	{
            double deltalogprob = -(This->*logprobf)();
            double m = tuning * (Random::Uniform() - 0.5);
            x += m;
            if (max > min)  {
                while ((x < min) || (x > max))  {
                    if (x < min)    {
                        x = 2*min - x;
                    }
                    if (x > max)    {
                        x = 2*max - x;
                    }
                }
            }
            (This->*updatef)();
            deltalogprob += (This->*logprobf)();
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted)	{
                nacc ++;
            }
            else	{
                x -= m;
                (This->*updatef)();
            }
            ntot++;
        }
        return nacc/ntot;
    }

	template<class C> double ScalingMove(double& x, double tuning, int nrep, LogProbF<C> logprobf, UpdateF<C> updatef, C* This) {
    
        double nacc = 0;
        double ntot = 0;
        for (int rep=0; rep<nrep; rep++)	{
            double deltalogprob = -(This->*logprobf)();
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            x *= e;
            (This->*updatef)();
            deltalogprob += (This->*logprobf)();
            deltalogprob += m;
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted)	{
                nacc ++;
            }
            else	{
                x /= e;
                (This->*updatef)();
            }
            ntot++;
        }
        return nacc/ntot;
    }

    template<class C> double ProfileMove(vector<double>& x, double tuning, int n, int nrep, LogProbF<C> logprobf, UpdateF<C> updatef, C* This)	{

        double nacc = 0;
        double ntot = 0;
        vector<double> bk(x.size(),0);
        for (int rep=0; rep<nrep; rep++)	{
            bk = x;
            double deltalogprob = -(This->*logprobf)();
            double loghastings = Random::ProfileProposeMove(x,x.size(),tuning,n);
            (This->*updatef)();
            deltalogprob += (This->*logprobf)();
            deltalogprob += loghastings;
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted)	{
                nacc ++;
            }
            else	{
                x = bk;
                (This->*updatef)();
            }
            ntot++;
        }
        return nacc/ntot;
    }
};

#endif  // PROBMODEL_H
