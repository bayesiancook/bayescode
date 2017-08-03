
#ifndef BETA_H
#define BETA_H

#include "Array.hpp"
#include "Random.hpp"
#include "PoissonSuffStat.hpp"
#include "CodonSuffStat.hpp"

class BetaSuffStat : public SuffStat	{

	public:
	BetaSuffStat() {}
	~BetaSuffStat() {}

	void Clear()	{
		sum0 = 0;
		sum1 = 0;
        n = 0;
	}

	void AddSuffStat(double log0, double log1, int c = 1)	{
		sum0 += log0;
		sum1 += log1;
		n += c;
	}

    double GetMeanInvConcLogProb(double mean, double invconc)   {
        double alpha = mean / invconc;
        double beta = (1-mean) / invconc;
        return GetLogProb(alpha,beta);
    }

	double GetLogProb(double alpha, double beta) const {
        return n * (Random::logGamma(alpha + beta) - Random::logGamma(alpha) - Random::logGamma(beta)) + (alpha-1)*sum0 + (beta-1)*sum1;
	}
	
    int GetN() {return n;}

	private:

	double sum0;
	double sum1;
	int n;
};

class IIDBeta : public SimpleArray<double> {

    public:

    IIDBeta(int insize, double inalpha, double inbeta) : SimpleArray<double>(insize), alpha(inalpha), beta(inbeta)  {
        Sample();
    }

    ~IIDBeta() {}

    double GetAlpha() const {return alpha;}
    double GetBeta() const {return beta;}

    double GetMean() const {return alpha / (alpha + beta);}
    double GetInvConcentration() const {return 1.0 / (alpha + beta);}

    void SetMeanInvConc(double mean, double invconc)    {
        double alpha = mean / invconc;
        double beta = (1-mean) / invconc;
        SetAlpha(alpha);
        SetBeta(beta);
    }

    void SetAlpha(double inalpha) {alpha = inalpha;}
    void SetBeta(double inbeta) {beta = inbeta;}

	void Sample()	{
		for (int i=0; i<GetSize(); i++)	{
            double a = Random::sGamma(alpha);
            double b = Random::sGamma(beta);
            (*this)[i] = a / (a+b);
		}
	}

	double GetLogProb()	{
		double total = 0;
		for (int i=0; i<GetSize(); i++)	{
			total += GetLogProb(i);
		}
		return total;
	}

	double GetLogProb(int i)	{
        return Random::logGamma(alpha + beta) - Random::logGamma(alpha) - Random::logGamma(beta) + (alpha-1)*log(GetVal(i)) + (beta-1)*log(1.0 - GetVal(i));
	}

	void AddSuffStat(BetaSuffStat& suffstat)	{
		for (int i=0; i<GetSize(); i++)	{
            suffstat.AddSuffStat(log(GetVal(i)), log(1-GetVal(i)), 1);
		}
	}

    double GetPosMean()    {
        int tot = 0;
        double m1 = 0;
        for (int i=0; i<GetSize(); i++) {
            if (GetVal(i))  {
                m1 += GetVal(i);
                tot++;
            }
        }
        if (! tot)  {
            return 0;
        }
        m1 /= tot;
        return m1;
    }

    /*
    double GetPosVar() {
        int tot = 0;
        double m1 = 0;
        double m2 = 0;
        for (int i=0; i<GetSize(); i++) {
            if (GetVal(i))  {
                m1 += GetVal(i);
                m2 += GetVal(i) * GetVal(i);
                tot++;
            }
        }
        if (! tot)  {
            return 0;
        }
        m1 /= tot;
        m2 /= tot;
        m2 -= m1*m1;
        return m2;
    }
    */

	protected:

    double alpha;
    double beta;
};


#endif


