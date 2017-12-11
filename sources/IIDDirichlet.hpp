
#ifndef IIDDIR_H
#define IIDDIR_H

#include "Array.hpp"
#include "Random.hpp"
#include "SuffStat.hpp"
#include "MPIBuffer.hpp"

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

    int GetDim()    {
        return (int) sumlog.size();
    }

    void AddSuffStat(const vector<double>& pi)  {
        for (unsigned int i=0; i<sumlog.size(); i++)    {
            if (pi[i] <= 0) {
                cerr << "error: negative pi in DirichletSuffStat: " << pi[i] << '\n';
                for (unsigned int j=0; j<sumlog.size(); j++)    {
                    cerr << pi[j] << '\t';
                }
                cerr << '\n';
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

    void Add(const DirichletSuffStat& from) {
        for (unsigned int i=0; i<sumlog.size(); i++)    {
            sumlog[i] += from.GetSumLog(i);
        }
        n += from.GetN();
    }

    DirichletSuffStat& operator+=(const DirichletSuffStat& from)    {
        Add(from);
        return *this;
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

    unsigned int GetMPISize() const {return sumlog.size() + 1;}

    void MPIPut(MPIBuffer& buffer) const {
        for (unsigned int i=0; i<sumlog.size(); i++)    {
            buffer << sumlog[i];
        }
        buffer << n;
    }

    void MPIGet(const MPIBuffer& buffer)    {
        for (unsigned int i=0; i<sumlog.size(); i++)    {
            buffer >> sumlog[i];
        }
        buffer >> n;
    }

    void Add(const MPIBuffer& buffer)   {
        double tmp;
        for (unsigned int i=0; i<sumlog.size(); i++)    {
            buffer >> tmp;
            sumlog[i] += tmp;
        }
        int temp;
        buffer >> temp;
        n += temp;
    }

	private:

    vector<double> sumlog;
    int n;
};


class DirichletSuffStatArray : public SimpleArray<DirichletSuffStat>    {

    public:
    DirichletSuffStatArray(int insize, int indim) : SimpleArray<DirichletSuffStat>(insize,DirichletSuffStat(indim)), dim(indim) {}
    ~DirichletSuffStatArray() {}

    int GetDim() const {
        return dim;
    }

    void Clear()    {
        for (int i=0; i<GetSize(); i++) {
            (*this)[i].Clear();
        }
    }

    void Add(const DirichletSuffStatArray& from)    {
		for (int i=0; i<GetSize(); i++)	{
            (*this)[i].Add(from.GetVal(i));
        }
    }

    DirichletSuffStatArray& operator+=(const DirichletSuffStatArray& from)  {
        Add(from);
        return *this;
    }

    unsigned int GetMPISize() const {return GetSize() * GetVal(0).GetMPISize();}

    void MPIPut(MPIBuffer& buffer) const    {
		for (int i=0; i<GetSize(); i++)	{
            buffer << GetVal(i);
        }
    }

    void MPIGet(const MPIBuffer& buffer)    {
		for (int i=0; i<GetSize(); i++)	{
            buffer >> (*this)[i];
        }
    }

    void Add(const MPIBuffer& buffer)   {
		for (int i=0; i<GetSize(); i++)	{
            (*this)[i] += buffer;
        }
    }

    private:
    int dim;
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

    void SetUniform()   {
        int dim = GetDim();
        for (int i=0; i<GetSize(); i++) {
            for (int k=0; k<dim; k++)   {
                (*this)[i][k] = 1.0/dim;
            }
        }
    }

    int GetDim() const {
        return center.size();
    }

	void Sample()	{
		for (int i=0; i<GetSize(); i++)	{
            Random::DirichletSample((*this)[i],center,concentration);
		}
	}

	double GetLogProb()	const {
		double total = 0;
		for (int i=0; i<GetSize(); i++)	{
			total += GetLogProb(i);
		}
		return total;
	}

	double GetLogProb(int i) const {
        return Random::logDirichletDensity(GetVal(i),center,concentration);
	}

	void AddSuffStat(DirichletSuffStat& suffstat) const {
		for (int i=0; i<GetSize(); i++)	{
			suffstat.AddSuffStat(GetVal(i));
		}
	}

	void AddSuffStat(DirichletSuffStat& suffstat, const Selector<int>& occupancy) const   {
		for (int i=0; i<GetSize(); i++)	{
            if (occupancy.GetVal(i))   {
                suffstat.AddSuffStat(GetVal(i));
            }
		}
	}

    void PriorResample(const Selector<int>& occupancy)    {
		for (int i=0; i<GetSize(); i++)	{
            if (! occupancy.GetVal(i)) {
                Random::DirichletSample((*this)[i],center,concentration);
            }
		}
    }

    double GetMeanEntropy() const   {

        double mean = 0;
        for (int i=0; i<GetSize(); i++) {
            mean += Random::GetEntropy(GetVal(i));
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

class MultiDirichlet: public SimpleArray<vector<double> >	{

	public: 

	MultiDirichlet(const Selector<vector<double> >* incenterarray, const Selector<double>* inconcentrationarray) : SimpleArray<vector<double> >(incenterarray->GetSize()), dim(incenterarray->GetVal(0).size()), centerarray(incenterarray), concentrationarray(inconcentrationarray) {
        if (centerarray->GetSize() != concentrationarray->GetSize())    {
            cerr << "error in multi dirichlet: center and concentration arrays should have same size\n";
            exit(1);
        }

        for (int i=0; i<GetSize(); i++) {
            (*this)[i].assign(dim,0);
        }
		Sample();
	}

	~MultiDirichlet() {}

    int GetDim() const {
        return dim;
    }

	void Sample()	{
		for (int i=0; i<GetSize(); i++)	{
            Random::DirichletSample((*this)[i],centerarray->GetVal(i),concentrationarray->GetVal(i));
		}
	}

    bool CheckPositivity()  {
        int allpos = 1;
		for (int i=0; i<GetSize(); i++)	{
            int pos = 1;
            for (int j=0; j<GetDim(); j++)  {
                if (! GetVal(i)[j]) {
                    pos = 0;
                }
            }
            if (! pos)  {
                allpos = 0;
                for (int j=0; j<GetDim(); j++)  {
                    cerr << GetVal(i)[j] << '\t';
                }
                cerr << '\n';
                cerr << "hyperconcentration: " << concentrationarray->GetVal(i) << '\n';
            }
        }
        return allpos;
    }

	double GetLogProb()	const {
		double total = 0;
		for (int i=0; i<GetSize(); i++)	{
			total += GetLogProb(i);
		}
		return total;
	}

	double GetLogProb(int i) const {
        return Random::logDirichletDensity(GetVal(i),centerarray->GetVal(i),concentrationarray->GetVal(i));
	}

    void AddSuffStat(Array<DirichletSuffStat>& suffstatarray, const Selector<int>& alloc) {
		for (int i=0; i<GetSize(); i++)	{
			suffstatarray[alloc.GetVal(i)].AddSuffStat(GetVal(i));
		}
	}

    void PriorResample(const Selector<int>& occupancy)    {
		for (int i=0; i<GetSize(); i++)	{
            if (! occupancy.GetVal(i)) {
                Random::DirichletSample((*this)[i],centerarray->GetVal(i),concentrationarray->GetVal(i));
                /*
                int pos = 1;
                for (int j=0; j<GetDim(); j++)  {
                    if (! GetVal(i)[j]) {
                        pos = 0;
                    }
                }
                if (! pos)  {
                    cerr << "in MultiDirichlet::PriorResample\n";
                    for (int j=0; j<GetDim(); j++)  {
                        cerr << GetVal(i)[j] << '\t';
                    }
                    cerr << '\n';
                    cerr << "hyperconcentration: " << concentrationarray->GetVal(i) << '\n';
                    exit(1);
                }
                */
            }
		}
    }


    double GetMeanEntropy() const   {

        double mean = 0;
        for (int i=0; i<GetSize(); i++) {
            mean += Random::GetEntropy(GetVal(i));
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
    int dim;
    const Selector<vector<double> >* centerarray;
    const Selector<double>* concentrationarray;
};

#endif
