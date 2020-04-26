
#pragma once

#include "Array.hpp"
#include "MPIBuffer.hpp"
#include "PoissonSuffStat.hpp"
#include "Random.hpp"

/**
 * \brief An array of omega across sites under M9 model
 * mixture of 
 * - point mass at 0
 * - beta distribution between 0 and 1
 * - point mass at 1
 * - shifted gamma distribution for values > 1
 */

class M9SuffStat : public SuffStat {
  public:
    M9SuffStat() : count(4,0) {}
    ~M9SuffStat() {}

    //! set suff stats to 0
    void Clear() {
        for (int k=0; k<4; k++) {
            count[k] = 0;
        }
        betasumlog0 = 0;
        betasumlog1 = 0;
        gammasumlog = 0;
        gammasum = 0;
    }

    void AddSuffStat(double omega)    {
        if (omega == 0) {
            count[0] ++;
        }
        else if (omega == 1.0)  {
            count[2] ++;
        }
        else if (omega > 1.0)   {
            count[3]++;
            double dposom = omega -1;
            gammasum += dposom;
            gammasumlog += log(dposom);
        }
        else    {
            count[1]++;
            betasumlog0 += log(omega);
            betasumlog1 += log(1.0-omega);
        }
    }

    void AddSuffStat(const Selector<double>& omegaarray)	{
	    for (int i=0; i<omegaarray.GetSize(); i++)	{
		    AddSuffStat(omegaarray.GetVal(i));
	    }
    }

    //! (*this) += from
    void Add(const M9SuffStat &from) {
        for (int k=0; k<4; k++) {
            count[k] = from.count[k];
        }
        betasumlog0 += from.betasumlog0;
        betasumlog1 += from.betasumlog1;
        gammasumlog += from.gammasumlog;
        gammasum += from.gammasum;
    }

    //! (*this) += from, operator version
    M9SuffStat &operator+=(const M9SuffStat &from) {
        Add(from);
        return *this;
    }

    //! return object size, when put into an MPI buffer
    unsigned int GetMPISize() const { return 8; }

    //! put object into MPI buffer
    void MPIPut(MPIBuffer &buffer) const { 
        buffer << count << betasumlog0 << betasumlog1 << gammasumlog << gammasum;
    }

    //! read object from MPI buffer
    void MPIGet(const MPIBuffer &buffer) {
        buffer >> count >> betasumlog0 >> betasumlog1 >> gammasumlog >> gammasum;
    }

    //! read a M9SuffStat from MPI buffer and add it to this
    void Add(const MPIBuffer &buffer) {
        int tmp;
        for (int k=0; k<4; k++) {
            buffer >> tmp;
            count[k] += tmp;
        }
        double temp;
        buffer >> temp;
        betasumlog0 += temp;
        buffer >> temp;
        betasumlog1 += temp;
        buffer >> temp;
        gammasumlog += temp;
        buffer >> temp;
        gammasum += temp;
    }

    //! return log prob of suff stats, as a function of parameters of the mixture
    double GetLogProb(const vector<double>& weight, double posw, double betamean, double betainvconc, double gammamean, double gammainvshape) const {

        double alpha = betamean / betainvconc;
        double beta = (1.0 - betamean) / betainvconc;
        double shape = 1.0 / gammainvshape;
        double scale = shape / gammamean;

        double logprob = 0;

        logprob += count[0] * log((1-posw)*weight[0]);

        logprob += count[1] * (log((1-posw)*weight[1]) 
                + Random::logGamma(alpha + beta) 
                - Random::logGamma(alpha)
                - Random::logGamma(beta))
            + (alpha-1) * betasumlog0 + (beta-1) * betasumlog1;

        logprob += count[2] * log((1-posw)*weight[2]);

        logprob += count[3] * (log(posw) 
                + shape * log(scale)
                - Random::logGamma(shape))
            + (shape-1)*gammasumlog - scale*gammasum;

        return logprob;
    }

    int GetCount(int k) const { return count[k]; }

  private:
    vector<int> count;
    double betasumlog0;
    double betasumlog1;
    double gammasumlog;
    double gammasum;
};

class IIDM9 : public SimpleArray<double> {
  public:
    //! constructor specifies array size and initial shape and scale parameters
    IIDM9(int insize, const vector<double>& inweight, double inposw, double betamean, double betainvconc, double gammamean, double gammainvshape) :
        SimpleArray<double>(insize),
        weight(inweight),
        posw(inposw),
        alpha(betamean / betainvconc),
        beta((1-betamean) / betainvconc),
        shape(1.0 / gammainvshape),
        scale(1.0 / gammainvshape / gammamean)  {
        Sample();
    }

    ~IIDM9() {}

    void SetParameters(double inposw, double betamean, double betainvconc, double gammamean, double gammainvshape)  {
        posw = inposw;
        alpha = betamean / betainvconc;
        beta = (1.0 - betamean) / betainvconc;
        shape = 1.0 / gammainvshape;
        scale = shape / gammamean;
    }

    //! sample all entries, given current shape and scale params
    void Sample() {
        for (int i = 0; i < GetSize(); i++) {
            Sample(i);
        }
    }

    void Sample(int i)  {
        if (Random::Uniform() < posw)   {
            (*this)[i] = 1.0 + Random::GammaSample(shape, scale);
        }
        else    {
            int cat = Random::DrawFromDiscreteDistribution(weight);
            if (! cat)  {
                (*this)[i] = 0;
            }
            else if (cat == 2)  {
                (*this)[i] = 1.0;
            }
            else    {
                (*this)[i] = Random::BetaSample(alpha,beta);
            }
        }
    }

    //! return total log prob (summed over the array) given current shape and
    //! scale params
    double GetLogProb() const {
        double total = 0;
        for (int i = 0; i < GetSize(); i++) {
            total += GetLogProb(i);
        }
        return total;
    }

    //! return log prob of one specific entry
    double GetLogProb(int index) const {
        double ret = 0;
        if (GetVal(index) == 0) {
            ret = log((1-posw)*weight[0]);
        }
        else if (GetVal(index) == 1.0)  {
            ret = log((1-posw)*weight[2]);
        }
        else if (GetVal(index) > 1.0)   {
            ret = log(posw) + Random::logGammaDensity(GetVal(index) - 1.0, shape, scale);
        }
        else    {
            ret = log((1-posw)*weight[1]) + Random::logBetaDensity(GetVal(index), alpha, beta);
        }
        return ret;
    }

    // never used
    /*
    //! resample all entries, given current shape and scale parameters and given
    //! an array of Poisson sufficient statistics of same size
    double MHMove(double tuning, int nrep,  const Selector<PoissonSuffStat> &suffstatarray) {
        double nacc = 0;
        double ntot = 0;
        for (int rep=0; rep<nrep; rep++)    {
            for (int i = 0; i < GetSize(); i++) {
                double& v = (*this)[i];
                if ((v != 0) && (v != 1.0)) {
                    const PoissonSuffStat &suffstat = suffstatarray.GetVal(i);
                    double logratio = -GetLogProb(i) - suffstat.GetCount()*log(v) + suffstat.GetBeta()*v;
                    double m = tuning*(Random::Uniform() - 0.5);
                    double e = exp(m);
                    v *= e;
                    logratio += GetLogProb(i) + suffstat.GetCount()*log(v) - suffstat.GetBeta()*v;
                    logratio += m;
                    int acc = (log(Random::Uniform()) < logratio);
                    if (acc)    {
                        nacc++;
                    }
                    else    {
                        v /= e;
                    }
                    ntot++;
                }
            }
        }
        return nacc / ntot;
    }
    */

    double MultipleTryMove(int nsample, const Selector<PoissonSuffStat> &suffstatarray) {

        double nacc = 0;

        vector<double> purom(nsample,0);
        vector<double> purlogprob(nsample,0);
        vector<double> purprob(nsample,0);

        vector<double> posom(nsample,0);
        vector<double> poslogprob(nsample,0);
        vector<double> posprob(nsample,0);

        vector<double> logprob(4,0);
        vector<double> postprob(4,0);

        for (int i=0; i<GetSize(); i++) {

            double& omega = (*this)[i];
            const PoissonSuffStat &suffstat = suffstatarray.GetVal(i);

            // if omega == 0 and there are substitutions, then log likelihood is -infty
            logprob[0] = suffstat.GetCount() ? log(0) : 0;
            // for omega == 1
            logprob[2] = - suffstat.GetBeta();

            // for 0 < omega < 1
            double purmax = 0;
            for (int k=0; k<nsample; k++)   {
                double om = 0;
                if ((!k) && (omega>0) && (omega<1.0)) {
                    om = omega;
                }
                else    {
                    om = Random::BetaSample(alpha,beta);
                }
                purlogprob[k] = suffstat.GetCount()*log(om) - suffstat.GetBeta()*om;
                purom[k] = om;
                if ((!k) || (purmax < purlogprob[k]))   {
                    purmax = purlogprob[k];
                }
            }

            double purtot = 0;
            for (int k=0; k<nsample; k++)   {
                purprob[k] = exp(purlogprob[k] - purmax);
                purtot += purprob[k];
            }
            logprob[1] = log(purtot/nsample) + purmax;
            for (int k=0; k<nsample; k++)   {
                purprob[k] /= purtot;
            }

            double posmax = 0;
            for (int k=0; k<nsample; k++)   {
                double om = 0;
                if ((!k) && (omega>1)) {
                    om = omega;
                }
                else    {
                    om = 1.0 + Random::GammaSample(shape,scale);
                }
                poslogprob[k] = suffstat.GetCount()*log(om) - suffstat.GetBeta()*om;
                posom[k] = om;
                if ((!k) || (posmax < poslogprob[k]))   {
                    posmax = poslogprob[k];
                }
            }

            double postot = 0;
            for (int k=0; k<nsample; k++)   {
                posprob[k] = exp(poslogprob[k] - posmax);
                postot += posprob[k];
            }
            logprob[3] = log(postot/nsample) + posmax;
            for (int k=0; k<nsample; k++)   {
                posprob[k] /= postot;
            }

            double max = logprob[0];
            for (int k=1; k<4; k++) {
                if (max < logprob[k])   {
                    max = logprob[k];
                }
            }

            postprob[0] = (1-posw)*weight[0] * exp(logprob[0] - max);
            postprob[1] = (1-posw)*weight[1] * exp(logprob[1] - max);
            postprob[2] = (1-posw)*weight[2] * exp(logprob[2] - max);
            postprob[3] = posw * exp(logprob[3] - max);

            double tot = 0;
            for (int k=0; k<4; k++) {
                tot += postprob[k];
            }
            for (int k=0; k<4; k++) {
                postprob[k] /= tot;
            }
            double bkomega = omega;
            int choose = Random::DrawFromDiscreteDistribution(postprob);
            if (choose == 0)    {
                omega = 0;
            }
            else if (choose == 1)   {
                int c = Random::DrawFromDiscreteDistribution(purprob);
                omega = purom[c];
            }
            else if (choose == 2)   {
                omega = 1.0;
            }
            else    {
                int c = Random::DrawFromDiscreteDistribution(posprob);
                omega = posom[c];
            }

            if (omega != bkomega)    {
                nacc++;
            }
        }

        return nacc / GetSize();
    }

    //! get mean over the array
    double GetMean() const {
        double m1 = 0;
        for (int i = 0; i < GetSize(); i++) {
            m1 += GetVal(i);
        }
        m1 /= GetSize();
        return m1;
    }

    //! get variance over the array
    double GetVar() const {
        double m1 = 0;
        double m2 = 0;
        for (int i = 0; i < GetSize(); i++) {
            m1 += GetVal(i);
            m2 += GetVal(i) * GetVal(i);
        }
        m1 /= GetSize();
        m2 /= GetSize();
        m2 -= m1 * m1;
        return m2;
    }

  protected:
    const vector<double>& weight;
    double posw;
    double alpha;
    double beta;
    double shape;
    double scale;
};

