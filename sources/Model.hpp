#ifndef MODEL_H
#define MODEL_H

#include "Random.hpp"

using namespace std;

class Model {

    public:

    template<class C> using LogProbF = double (C::*)(void);
    template<class C> using UpdateF = void (C::*)(void);
    // template<class C> typedef double (C::*LogProbF)(void);
    // template<class C> typedef void (C::*UpdateF)(void);

	template<class C> double SlidingMove(double& x, double tuning, int nrep, double min, double max, LogProbF<C> logprobf, UpdateF<C> updatef, C* This) {
    
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
            }
            ntot++;
        }
        return nacc/ntot;
    }
};

#endif