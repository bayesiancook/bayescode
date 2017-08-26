
#ifndef IIDDIR_H
#define IIDDIR_H

#include "Array.hpp"
#include "Random.hpp"
#include "SuffStat.hpp"

class DirichletSuffStat : public SuffStat	{

	public:
    DirichletSuffStat(int indim) : sumlog(indim,0), n(0) {}
    ~DirichletSuffStat() {}

	void Clear()	{
        for (unsigned int i=0; i<sumlog.size(); i++)    {
            sumlog[i] = 0;
        }
        n = 0;
	}

    void AddSuffStat(const vector<double>& pi)  {
        for (unsigned int i=0; i<sumlog.size(); i++)    {
            if (pi[i] <= 0) {
                cerr << "error: negative pi in DirichletSuffStat: " << pi[i] << '\n';
                exit(1);
            }
            sumlog[i] += log(pi[i]);
        }
        n++;
    }

    void AddSuffStat(const double* insumlog, int d)  {
        for (unsigned int i=0; i<sumlog.size(); i++)    {
            sumlog[i] += insumlog[i];
        }
        n += d;
    }

    double GetSumLog(int i) const   {
        return sumlog[i];
    }

    int GetN() const    {
        return n;
    }

	double GetLogProb(const vector<double>& center, double concentration) const    {
        
        double tot = n * Random::logGamma(concentration);
        for (unsigned int i=0; i<sumlog.size(); i++)    {
            tot += - n * Random::logGamma(concentration*center[i]) + (concentration*center[i]-1)*sumlog[i];
        }
        return tot;
    }

	private:

    vector<double> sumlog;
    int n;
};

class IIDDirichlet: public SimpleArray<vector<double> >	{

	public: 

	IIDDirichlet(int insize, const vector<double>& incenter, double inconcentration) : SimpleArray<vector<double> >(insize), center(incenter), concentration(inconcentration) {
        for (int i=0; i<GetSize(); i++) {
            (*this)[i].assign(center.size(),0);
        }
		Sample();
	}

	~IIDDirichlet() {}

    void SetCenter(const vector<double>& incenter)    {
        center = incenter;
    }
    
    void SetConcentration(double inconcentration)   {
        concentration = inconcentration;
    }

    int GetDim() const {
        return center.size();
    }

	void Sample()	{
		for (int i=0; i<GetSize(); i++)	{
            double total = 0;
            for (int k=0; k<GetDim(); k++)  {
                (*this)[i][k] = Random::sGamma(concentration*center[k]);
                total += (*this)[i][k];
            }
            for (int k=0; k<GetDim(); k++)  {
                (*this)[i][k] /= total;
            }
		}
	}

    // to be rewritten: multinomial suffstat given as argument
    /*
	void GibbsResample(const ConstArray<MultinomialSuffStat>& suffstatarray)	{
		for (int i=0; i<GetSize(); i++)	{
			const PoissonSuffStat& suffstat = suffstatarray.GetVal(i);
			(*this)[i] = Random::Gamma(shape + suffstat.GetCount(), scale + suffstat.GetBeta());
		}
	}
    */

	double GetLogProb()	const {
		double total = 0;
		for (int i=0; i<GetSize(); i++)	{
			total += GetLogProb(i);
		}
		return total;
	}

	double GetLogProb(int i) const {
        double tot = Random::logGamma(concentration);
        for (int k=0; k<GetDim(); k++)  {
            tot += -Random::logGamma(concentration*center[k]) + (concentration*center[k]-1)*log(GetVal(i)[k]);
        }
        return tot;
	}

	void AddSuffStat(DirichletSuffStat& suffstat) const {
		for (int i=0; i<GetSize(); i++)	{
			suffstat.AddSuffStat(GetVal(i));
		}
	}

	double GetEntropy(const std::vector<double>& profile) const {
		double tot = 0;
		for (unsigned int i=0; i<profile.size(); i++)	{
			tot -= (profile[i] < 1e-6) ? 0 : profile[i]*log(profile[i]);
		}
		return tot;
	}

    double GetMeanEntropy() const   {

        double mean = 0;
        for (int i=0; i<GetSize(); i++) {
            mean += GetEntropy(GetVal(i));
        }
        mean /= GetSize();
        return mean;
    }

    double GetMean(int k) const {
        double m1 = 0;
        for (int i=0; i<GetSize(); i++) {
            m1 += GetVal(i)[k];
        }
        m1 /= GetSize();
        return m1;
    }

    double GetVar(int k) const {
        double m1 = 0;
        double m2 = 0;
        for (int i=0; i<GetSize(); i++) {
            m1 += GetVal(i)[k];
            m2 += GetVal(i)[k] * GetVal(i)[k];
        }
        m1 /= GetSize();
        m2 /= GetSize();
        m2 -= m1*m1;
        return m2;
    }

	protected:
    vector<double> center;
    double concentration;
};

#endif
