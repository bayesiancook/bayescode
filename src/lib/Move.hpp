#pragma once

#include "global/Random.hpp"

class Move {
  public:
    //! new type name for a const method of class C taking no argument and
    //! returning a double (intended: a log prob function)
    template <class C>
    using LogProbF = double (C::*)(void) const;
    //! new type name for a non-const method of class C taking no argument and
    //! with no return (intended: an update function)
    template <class C>
    using UpdateF = void (C::*)(void);

    //! new type name for a const method of class C taking no argument and
    //! returning a double (intended: a log prob function)
    template <class C>
    using VectorLogProbF = double (C::*)(int) const;
    //! new type name for a non-const method of class C taking no argument and
    //! with no return (intended: an update function)
    template <class C>
    using VectorUpdateF = void (C::*)(int);

    //! \brief template for Metropolis Hastings sliding move on parameter x
    //!
    //! Parameterized by a tuning parameter, un number of iterations (nrep),
    //! and lower and upper constraint (min and max; if max<=min, then no
    //! constraint is enforced). Should also give a pointer to a log prob and an
    //! update functions, as well as a pointer to the model itself. Returns
    //! success rate.
    template <class C>
    static double Sliding(double &x, double tuning, int nrep, double min, double max,
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

    //! \brief template for Metropolis Hastings sliding move on parameter vector x
    //!
    //! Parameterized by a tuning parameter, un number of iterations (nrep),
    //! and lower and upper constraint (min and max; if max<=min, then no
    //! constraint is enforced). Should also give a pointer to a log prob and an
    //! update functions, as well as a pointer to the model itself. Returns
    //! success rate.
    template <class C>
    static double VectorSliding(std::vector<double> &x, double tuning, int nrep, double min,
        double max, VectorLogProbF<C> logprobf, VectorUpdateF<C> updatef, C *This) {
        double nacc = 0;
        double ntot = 0;
        for (int rep = 0; rep < nrep; rep++) {
            for (unsigned int i = 0; i < x.size(); i++) {
                double bk = x[i];
                double deltalogprob = -(This->*logprobf)(i);
                double m = tuning * (Random::Uniform() - 0.5);
                x[i] += m;
                if (max > min) {
                    while ((x[i] < min) || (x[i] > max)) {
                        if (x[i] < min) { x[i] = 2 * min - x[i]; }
                        if (x[i] > max) { x[i] = 2 * max - x[i]; }
                    }
                }
                (This->*updatef)(i);
                deltalogprob += (This->*logprobf)(i);
                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    x[i] = bk;
                    (This->*updatef)(i);
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    //! \brief template for Metropolis Hastings sliding move on parameter vector x
    //!
    //! Parameterized by a tuning parameter, un number of iterations (nrep),
    //! and lower and upper constraint (min and max; if max<=min, then no
    //! constraint is enforced). Should also give a pointer to a log prob and an
    //! update functions, as well as a pointer to the model itself. Returns
    //! success rate.
    template <class C>
    static double VectorSliding(std::vector<double> &x, double tuning, int nrep,
        const std::vector<double> &min, const std::vector<double> &max, VectorLogProbF<C> logprobf,
        VectorUpdateF<C> updatef, C *This) {
        assert((x.size() == min.size()) && (min.size() == max.size()));

        double nacc = 0;
        double ntot = 0;
        for (int rep = 0; rep < nrep; rep++) {
            for (unsigned int i = 0; i < x.size(); i++) {
                double bk = x[i];
                double deltalogprob = -(This->*logprobf)(i);
                double m = tuning * (Random::Uniform() - 0.5);
                x[i] += m;
                if (max[i] > min[i]) {
                    while ((x[i] < min[i]) || (x[i] > max[i])) {
                        if (x[i] < min[i]) { x[i] = 2 * min[i] - x[i]; }
                        if (x[i] > max[i]) { x[i] = 2 * max[i] - x[i]; }
                    }
                }
                (This->*updatef)(i);
                deltalogprob += (This->*logprobf)(i);
                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    x[i] = bk;
                    (This->*updatef)(i);
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    //! \brief template for Metropolis Hastings scaling move on parameter x
    //!
    //! Parameterized by a tuning parameter, un number of iterations (nrep),
    //! Should give a pointer to a log prob and an update function.
    //! Return success rate.
    template <class C>
    static double Scaling(
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

    //! \brief template for Metropolis Hastings scaling move on parameter vector x
    //!
    //! Parameterized by a tuning parameter, un number of iterations (nrep),
    //! Should give a pointer to a log prob and an update function.
    //! Return success rate.
    template <class C>
    static double VectorScaling(std::vector<double> &x, double tuning, int nrep,
        VectorLogProbF<C> logprobf, VectorUpdateF<C> updatef, C *This) {
        double nacc = 0;
        double ntot = 0;
        for (int rep = 0; rep < nrep; rep++) {
            for (unsigned int i = 0; i < x.size(); i++) {
                double deltalogprob = -(This->*logprobf)(i);
                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                x[i] *= e;
                (This->*updatef)(i);
                deltalogprob += (This->*logprobf)(i);
                deltalogprob += m;
                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    x[i] /= e;
                    (This->*updatef)(i);
                }
                ntot++;
            }
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
    static double Profile(std::vector<double> &x, double tuning, int n, int nrep,
        LogProbF<C> logprobf, UpdateF<C> updatef, C *This) {
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
