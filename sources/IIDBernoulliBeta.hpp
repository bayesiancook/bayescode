
#ifndef BERNBETA_H
#define BERNBETA_H

#include "Array.hpp"
#include "Random.hpp"
#include "PoissonSuffStat.hpp"
#include "CodonSuffStat.hpp"

class BernoulliBetaSuffStat : public SuffStat	{

	public:
	BernoulliBetaSuffStat() {}
	~BernoulliBetaSuffStat() {}

	void Clear()	{
		sum0 = 0;
		sum1 = 0;
        posn = 0;
		nulln = 0;
	}

    void AddNullSuffStat(int c = 1) {
        nulln += c;
    }

	void AddPosSuffStat(double log0, double log1, int c = 1)	{
		sum0 += log0;
		sum1 += log1;
		posn += c;
	}

	double GetLogProb(double pi, double alpha, double beta) const {
        double logbern = nulln * log(1-pi) + posn * log(pi);
        double logbeta = posn * (Random::logGamma(alpha + beta) - Random::logGamma(alpha) - Random::logGamma(beta)) + (alpha-1)*sum0 + (beta-1)*sum1;
        return logbern + logbeta;
	}
	
    int GetN0() const {return nulln;}
    int GetN1() const {return posn;}

	private:

	double sum0;
	double sum1;
	int nulln;
	int posn;
};

class IIDBernoulliBeta : public SimpleArray<double> {

    public:

    IIDBernoulliBeta(int insize, double inpi, double inalpha, double inbeta) : SimpleArray<double>(insize), pi(inpi), alpha(inalpha), beta(inbeta)  {
        Sample();
    }

    ~IIDBernoulliBeta() {}

    double GetPi() const {return pi;}
    double GetAlpha() const {return alpha;}
    double GetBeta() const {return beta;}

    void SetPi(double inpi) {pi = inpi;}
    void SetAlpha(double inalpha) {alpha = inalpha;}
    void SetBeta(double inbeta) {beta = inbeta;}

	void Sample()	{
		for (int i=0; i<GetSize(); i++)	{
            if (Random::Uniform() < pi) {
                double a = Random::sGamma(alpha);
                double b = Random::sGamma(beta);
                (*this)[i] = a / (a+b);
            }
            else    {
                (*this)[i] = 0;
            }
		}
	}

    int GetNullSet() {
        int tot = 0;
        for (int i=0; i<GetSize(); i++) {
            if (!GetVal(i)) {
                tot++;
            }
        }
        return tot;
    }

	double GetLogProb()	{
		double total = 0;
		for (int i=0; i<GetSize(); i++)	{
			total += GetLogProb(i);
		}
		return total;
	}

	double GetLogProb(int i)	{
        double ret = 0;
        if (! GetVal(i))    {
            ret += log(1-pi);
        }
        else    {
            ret += log(pi);
            ret += Random::logGamma(alpha + beta) - Random::logGamma(alpha) - Random::logGamma(beta) + (alpha-1)*log(GetVal(i)) + (beta-1)*log(1.0 - GetVal(i));
        }
        return ret;
	}

	void AddSuffStat(BernoulliBetaSuffStat& suffstat)	{
		for (int i=0; i<GetSize(); i++)	{
            if (! GetVal(i))    {
                suffstat.AddNullSuffStat(1);
            }
            else    {
                suffstat.AddPosSuffStat(log(GetVal(i)), log(1-GetVal(i)), 1);
            }
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

    double pi;
    double alpha;
    double beta;
};


#endif


