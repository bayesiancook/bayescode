
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
		sumlog0 = 0;
		sumlog1 = 0;
        n0 = 0;
		n1 = 0;
	}

    void AddNullSuffStat(int c = 1) {
        n0 += c;
    }

	void AddPosSuffStat(double log0, double log1, int c = 1)	{
		sumlog0 += log0;
		sumlog1 += log1;
		n1 += c;
	}

	double GetLogProb(double pi, double alpha, double beta) const {
        double logbern = n0 * log(1-pi) + n1 * log(pi);
        double logbeta = n1 * (Random::logGamma(alpha + beta) - Random::logGamma(alpha) - Random::logGamma(beta)) + (alpha-1)*sumlog0 + (beta-1)*sumlog1;
        return logbern + logbeta;
	}
	
    int GetN0() const {return n0;}
    int GetN1() const {return n1;}
    int GetSumLog0() const {return sumlog0;}
    int GetSumLog1() const {return sumlog1;}

	private:

	double sumlog0;
	double sumlog1;
	int n0;
	int n1;
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


