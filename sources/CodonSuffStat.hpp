
#ifndef CODONSUFFSTAT_H
#define CODONSUFFSTAT_H

#include "PathSuffStat.hpp"
#include "CodonSubMatrixArray.hpp"
#include "PoissonSuffStat.hpp"
#include <typeinfo>

class NucPathSuffStat : public SuffStat	{

	public:

	NucPathSuffStat() :
        rootcount(4,0),
        paircount(4,vector<int>(4,0)),
        pairbeta(4,vector<double>(4,0)) {}

	~NucPathSuffStat() {}

	void Clear()	{
		for (int i=0; i<Nnuc; i++)	{
			rootcount[i] = 0;
			for (int j=0; j<Nnuc; j++)	{
				paircount[i][j] = 0;
				pairbeta[i][j] = 0;
			}
		}
	}

	// assumes pathsuffstat is 61x61
	// collect the 4x4 path suff stat out of codonpathsuffstat
	void AddSuffStat(const MGOmegaCodonSubMatrix& codonmatrix, const PathSuffStat& codonpathsuffstat)    {

		const CodonStateSpace* cod = codonmatrix.GetCodonStateSpace();
        const SubMatrix* nucmatrix = codonmatrix.GetNucMatrix();

		// root part
        const std::map<int,int>& codonrootcount = codonpathsuffstat.GetRootCountMap();
        for (std::map<int,int>::const_iterator i = codonrootcount.begin(); i!= codonrootcount.end(); i++)	{
            int codon = i->first;
            rootcount[cod->GetCodonPosition(0,codon)] += i->second;
            rootcount[cod->GetCodonPosition(1,codon)] += i->second;
            rootcount[cod->GetCodonPosition(2,codon)] += i->second;
        }

        const std::map<int,double>& waitingtime = codonpathsuffstat.GetWaitingTimeMap();
        for (std::map<int,double>::const_iterator i = waitingtime.begin(); i!= waitingtime.end(); i++)	{
            int codon = i->first;
            for (int c2 = 0; c2 < cod->GetNstate(); c2++)	{
                if (c2 != codon)	{
                    int pos = cod->GetDifferingPosition(codon,c2);
                    if (pos < 3)	{
                        int n1 = cod->GetCodonPosition(pos,codon);
                        int n2 = cod->GetCodonPosition(pos,c2);
                        pairbeta[n1][n2] += i->second * codonmatrix(codon,c2) / (*nucmatrix)(n1,n2);
                    }
                }
            }
        }

        const std::map<pair<int,int>,int>& codonpaircount = codonpathsuffstat.GetPairCountMap();
        for (std::map<pair<int,int>, int>::const_iterator i = codonpaircount.begin(); i!= codonpaircount.end(); i++)	{
            int cod1 = i->first.first;
            int cod2 = i->first.second;
            int pos = cod->GetDifferingPosition(cod1,cod2);
            if (pos == 3)	{
                cerr << "error in codon conj path suffstat\n";
                exit(1);
            }
            int n1 = cod->GetCodonPosition(pos,cod1);
            int n2 = cod->GetCodonPosition(pos,cod2);
            paircount[n1][n2] += i->second;
        }
	}

	void AddSuffStat(const MGOmegaCodonSubMatrixArray& codonmatrixarray, const PathSuffStatArray& codonpathsuffstatarray)    {
		for (int i=0; i<codonmatrixarray.GetSize(); i++)	{
            AddSuffStat(codonmatrixarray.GetVal(i),codonpathsuffstatarray.GetVal(i));
        }
    }

	double GetLogProb(const SubMatrix& mat, const CodonStateSpace& cod) const {

        double total = 0;
		// root part
        int nroot = 0;
        const double* rootstat = mat.GetStationary();
        for (int i=0; i<Nnuc; i++)	{
            total += rootcount[i] * log(rootstat[i]);
            nroot += rootcount[i];
        }
        total -= nroot / 3 * log(cod.GetNormStat(rootstat));

		// non root part
        for (int i=0; i<Nnuc; i++)	{
            for (int j=0; j<Nnuc; j++)	{
                if (i != j)	{
                    total += paircount[i][j] * log(mat(i,j));
                    total -= pairbeta[i][j] * mat(i,j);
                }
            }
        }

		return total;
	}

    void Add(const NucPathSuffStat& from)   {

        for (int i=0; i<Nnuc; i++)  {
            rootcount[i] += from.rootcount[i];
        }
        for (int i=0; i<Nnuc; i++)  {
            for (int j=0; j<Nnuc; j++)  {
                paircount[i][j] += from.paircount[i][j];
                pairbeta[i][j] += from.pairbeta[i][j];
            }
        }
    }

    void Push(int* incount, double* inbeta) const {
        int index = 0;
        for (int i=0; i<Nnuc; i++)  {
            incount[index++] = rootcount[i];
        }
        for (int i=0; i<Nnuc; i++)  {
            for (int j=0; j<Nnuc; j++)  {
                incount[index++] = paircount[i][j];
            }
        }
        index = 0;
        for (int i=0; i<Nnuc; i++)  {
            for (int j=0; j<Nnuc; j++)  {
                inbeta[index++] = pairbeta[i][j];
            }
        }
    }

    void Add(const int* incount, const double* inbeta)  {
        int index = 0;
        for (int i=0; i<Nnuc; i++)  {
            rootcount[i] += incount[index++];
        }
        for (int i=0; i<Nnuc; i++)  {
            for (int j=0; j<Nnuc; j++)  {
                paircount[i][j] += incount[index++];
            }
        }
        index = 0;
        for (int i=0; i<Nnuc; i++)  {
            for (int j=0; j<Nnuc; j++)  {
                pairbeta[i][j] += inbeta[index++];
            }
        }
    }

    private:

    std::vector<int> rootcount;
    std::vector<vector<int> > paircount;
    std::vector<vector<double> > pairbeta;
};

class OmegaSuffStat : public PoissonSuffStat {

	public:

	OmegaSuffStat() {}
	~OmegaSuffStat() {}

	// assumes pathsuffstat is 61x61
	// tease out syn and non-syn substitutions and sum up count and beta stats  
	void AddSuffStat(const MGOmegaCodonSubMatrix& codonsubmatrix, const PathSuffStat& pathsuffstat)	{

        int ncodon = codonsubmatrix.GetNstate();
        const CodonStateSpace* statespace = codonsubmatrix.GetCodonStateSpace();

        const std::map<pair<int,int>,int>& paircount = pathsuffstat.GetPairCountMap();
        const std::map<int,double>& waitingtime = pathsuffstat.GetWaitingTimeMap();

        double tmpbeta = 0;
        for (std::map<int,double>::const_iterator i = waitingtime.begin(); i!= waitingtime.end(); i++)	{
            double totnonsynrate = 0;
            int a = i->first;
            for (int b=0; b<ncodon; b++)	{
                if (b != a)	{
                    if (codonsubmatrix(a,b) != 0)	{
                        if (!statespace->Synonymous(a,b))	{
                            totnonsynrate += codonsubmatrix(a,b);
                        }
                    }
                }
            }
            tmpbeta += i->second * totnonsynrate;
        }
        tmpbeta /= codonsubmatrix.GetOmega();

        int tmpcount = 0;
        for (std::map<pair<int,int>, int>::const_iterator i = paircount.begin(); i!= paircount.end(); i++)	{
            if (! statespace->Synonymous(i->first.first,i->first.second))	{
                tmpcount += i->second;
            }
        }

        PoissonSuffStat::AddSuffStat(tmpcount,tmpbeta);
    }
};

class OmegaSuffStatArray : public SimpleArray<OmegaSuffStat>, public Array<PoissonSuffStat>    {

	public:

	OmegaSuffStatArray(int insize) : SimpleArray<OmegaSuffStat>(insize) {}
	~OmegaSuffStatArray() {}

    int GetSize() const override {return array.size();}
    const OmegaSuffStat& GetVal(int i) const override {return array[i];}
    OmegaSuffStat& operator[](int i) override {return array[i];}

	void Clear()	{
		for (int i=0; i<GetSize(); i++)	{
			(*this)[i].Clear();
		}
	}

	void AddSuffStat(const ConstArray<MGOmegaCodonSubMatrix>& codonsubmatrixarray, const ConstArray<PathSuffStat>& pathsuffstatarray)	{
		for (int i=0; i<GetSize(); i++)	{
            (*this)[i].AddSuffStat(codonsubmatrixarray.GetVal(i),pathsuffstatarray.GetVal(i));
		}
	}

	double GetLogProb(const Array<double>* omegaarray) const{
		double total = 0;
		for (int i=0; i<GetSize(); i++)	{
			total += GetVal(i).GetLogProb(omegaarray->GetVal(i));
		}
		return total;
	}

	double GetMarginalLogProb(double shape, double scale)	const {
		double total = 0;
		// factoring out prior factor
		for (int i=0; i<GetSize(); i++)	{
			int count = GetVal(i).GetCount();
			double beta = GetVal(i).GetBeta();
			total += -(shape+count)*log(scale+beta) + Random::logGamma(shape+count);
		}
		total += GetSize() * (shape*log(scale) - Random::logGamma(shape));
		return total;
	}
};

#endif
