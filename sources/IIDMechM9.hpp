
#pragma once

#include "Array.hpp"
#include "MPIBuffer.hpp"
#include "PoissonSuffStat.hpp"
#include "Random.hpp"

/**
 * \brief An array of omega across sites under MechM9 model
 * mixture of 
 * - point mass at 0
 * - beta distribution between 0 and 1
 * - point mass at 1
 * - shifted gamma distribution for values > 1
 */

class MechM9SuffStat : public SuffStat {
  public:
    MechM9SuffStat() : count(2,0) {}
    ~MechM9SuffStat() {}

    //! set suff stats to 0
    void Clear() {
        for (int k=0; k<2; k++) {
            count[k] = 0;
        }
        gammasumlog = 0;
        gammasum = 0;
    }

    void AddSuffStat(double omega)    {
        if (omega == 1) {
            count[0] ++;
        }
        else    {
            count[1]++;
            double dposom = omega -1;
            gammasum += dposom;
            gammasumlog += log(dposom);
        }
    }

    void AddSuffStat(const Selector<double>& omegaarray)	{
	    for (int i=0; i<omegaarray.GetSize(); i++)	{
		    AddSuffStat(omegaarray.GetVal(i));
	    }
    }

    //! (*this) += from
    void Add(const MechM9SuffStat &from) {
        for (int k=0; k<2; k++) {
            count[k] = from.count[k];
        }
        gammasumlog += from.gammasumlog;
        gammasum += from.gammasum;
    }

    //! (*this) += from, operator version
    MechM9SuffStat &operator+=(const MechM9SuffStat &from) {
        Add(from);
        return *this;
    }

    //! return object size, when put into an MPI buffer
    unsigned int GetMPISize() const { return 4; }

    //! put object into MPI buffer
    void MPIPut(MPIBuffer &buffer) const { 
        buffer << count << gammasumlog << gammasum;
    }

    //! read object from MPI buffer
    void MPIGet(const MPIBuffer &buffer) {
        buffer >> count >> gammasumlog >> gammasum;
    }

    //! read a MechM9SuffStat from MPI buffer and add it to this
    void Add(const MPIBuffer &buffer) {
        int tmp;
        for (int k=0; k<2; k++) {
            buffer >> tmp;
            count[k] += tmp;
        }
        double temp;
        buffer >> temp;
        gammasumlog += temp;
        buffer >> temp;
        gammasum += temp;
    }

    //! return log prob of suff stats, as a function of parameters of the mixture
    double GetLogProb(double posw, double gammamean, double gammainvshape) const {

        double shape = 1.0 / gammainvshape;
        double scale = shape / gammamean;

        double logprob = 0;

        logprob += count[0] * log(1-posw);

        logprob += count[1] * (log(posw) 
                + shape * log(scale)
                - Random::logGamma(shape))
            + (shape-1)*gammasumlog - scale*gammasum;

        return logprob;
    }

    int GetCount(int k) const { return count[k]; }

  private:
    vector<int> count;
    double gammasumlog;
    double gammasum;
};

class IIDMechM9 : public SimpleArray<double> {
  public:
    //! constructor specifies array size and initial shape and scale parameters
    IIDMechM9(int insize, double inposw, double gammamean, double gammainvshape) :
        SimpleArray<double>(insize),
        posw(inposw),
        shape(1.0 / gammainvshape),
        scale(1.0 / gammainvshape / gammamean)  {
        Sample();
    }

    ~IIDMechM9() {}

    void SetParameters(double inposw, double gammamean, double gammainvshape)  {
        posw = inposw;
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
            (*this)[9] = 1.0;
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
        if (GetVal(index) == 1.0) {
            ret = log(1-posw);
        }
        else    {
            ret = log(posw) + Random::logGammaDensity(GetVal(index) - 1.0, shape, scale);
        }
        return ret;
    }

    double MultipleTryMove(int nsample, const Selector<PoissonSuffStat> &suffstatarray) {

        double nacc = 0;

        vector<double> posom(nsample,0);
        vector<double> poslogprob(nsample,0);
        vector<double> posprob(nsample,0);

        vector<double> logprob(2,0);
        vector<double> postprob(2,0);

        for (int i=0; i<GetSize(); i++) {

            double& omega = (*this)[i];
            const PoissonSuffStat &suffstat = suffstatarray.GetVal(i);

            logprob[0] = - suffstat.GetBeta();

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

            double max = (logprob[0] > logprob[1]) ? logprob[0] : logprob[1];
            postprob[0] = (1-posw) * exp(logprob[0] - max);
            postprob[1] = posw * exp(logprob[1] - max);
            double tot = postprob[0] + postprob[1];
            postprob[0] /= tot;
            postprob[1] /= tot;

            double bkomega = omega;
            int choose = Random::DrawFromDiscreteDistribution(postprob);
            if (choose == 0)    {
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
    double posw;
    double shape;
    double scale;
};

