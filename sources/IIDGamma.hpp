
#ifndef IIDGAMMA_H
#define IIDGAMMA_H

#include "Array.hpp"
#include "BranchArray.hpp"
#include "Random.hpp"
#include "PoissonSuffStat.hpp"
#include "MPIBuffer.hpp"

/**
 * \brief A sufficient statistic for a collection of gamma variates, as a function of the shape and scale parameters
 *
 * Suppose you have x = (x_i)_i=1..N, iid Gamma(shape,scale).
 * Then, p(X | shape,scale) can be expressed as a function of compact sufficient statistics: sum x_i's, sum log(x_i)'s and N.
 * GammaSuffStat implements this idea, by providing methods for collecting these suff stats and returning the log prob for a given
 * value for the shape and scale parameters.
 */

class GammaSuffStat : public SuffStat	{

	public:
	GammaSuffStat() {}
	~GammaSuffStat() {}

    //! set suff stats to 0
	void Clear()	{
		sum = 0;
		sumlog = 0;
		n = 0;
	}

    //! add the contribution of one gamma variate (x) to this suffstat
	void AddSuffStat(double x, double logx, int c = 1)	{
		sum += x;
		sumlog += logx;
		n += c;
	}

    //! (*this) += from
    void Add(const GammaSuffStat& from) {
        sum += from.GetSum();
        sumlog += from.GetSumLog();
        n += from.GetN();
    }

    //! (*this) += from, operator version
    GammaSuffStat& operator+=(const GammaSuffStat& from)    {
        Add(from);
        return *this;
    }

    //! return object size, when put into an MPI buffer
    unsigned int GetMPISize() const {return 3;}

    //! put object into MPI buffer
    void MPIPut(MPIBuffer& buffer) const    {
        buffer << sum << sumlog << n;
    }

    //! read object from MPI buffer
    void MPIGet(const MPIBuffer& buffer)    {
        buffer >> sum >> sumlog >> n;
    }

    //! read a GammaSuffStat from MPI buffer and add it to this
    void Add(const MPIBuffer& buffer)   {
        double temp;
        buffer >> temp;
        sum += temp;
        buffer >> temp;
        sumlog += temp;

        int tmp;
        buffer >> tmp;
        n += tmp;
    }

    //! return log prob, as a function of the given shape and scale parameters
	double GetLogProb(double shape, double scale) const {
		return n*(shape*log(scale) - Random::logGamma(shape)) + (shape-1)*sumlog - scale*sum;
	}
	
    //! return sum x_i's
    double GetSum() const {return sum;}
    //! return sum log x_i's
    double GetSumLog() const {return sumlog;}
    //! return N, total number of gamma variates contributing to the suff stat
    int GetN() const {return n;}

	private:

	double sum;
	double sumlog;
	int n;
};

/**
 * \brief A tree-structured branch-wise array of gamma sufficient statistics
 *
 * Useful for gene-specific branch lengths (see for instance MultiGeneCodonM2aModel).
 * In this model,
 * for a given branch of the tree,
 * genes have differing lengths, which are iid Gamma, of shape and scale that are branch-specific.
 * Thus, it is useful to collect suff stats across genes, and do this branchwise (storing the result into a branch-wise array). 
 * This array of suffstats can then be used to do fast MCMC moves on the hyperparameters of the distribution of 
 * branch lengths across genes.
 */

class GammaSuffStatBranchArray : public SimpleBranchArray<GammaSuffStat>    {

    public:
	GammaSuffStatBranchArray(const Tree& intree) : SimpleBranchArray<GammaSuffStat>(intree) {}
    ~GammaSuffStatBranchArray() {}

    //! member-wise addition between the two arrays ((*this) += from)
    void Add(const GammaSuffStatBranchArray& from)  {
        for (int i=0; i<GetNbranch(); i++)  {
            (*this)[i].Add(from.GetVal(i));
        }
    }

    //! member-wise addition between the two arrays, operator version
    GammaSuffStatBranchArray& operator+=(const GammaSuffStatBranchArray& from)  {
        Add(from);
        return *this;
    }

    //! return array size, when put into an MPI buffer
    unsigned int GetMPISize() const {return 3*GetNbranch();}

    //! put array into MPI buffer
    void MPIPut(MPIBuffer& buffer) const    {
        for (int i=0; i<GetNbranch(); i++)  {
            buffer << GetVal(i);
        }
    }

    //! read array from MPI buffer
    void MPIGet(const MPIBuffer& buffer)    {
        for (int i=0; i<GetNbranch(); i++)  {
            buffer >> (*this)[i];
        }
    }

    //! read from MPI buffer and add to current array
    void Add(const MPIBuffer& buffer)    {
        for (int i=0; i<GetNbranch(); i++)  {
            (*this)[i].Add(buffer);
        }
    }

    //! set all suff stats to 0
    void Clear()    {
        for (int i=0; i<GetNbranch(); i++)  {
            (*this)[i].Clear();
        }
    }

    //! get total log prob
    double GetLogProb(const BranchSelector<double>& blmean, double invshape) const {

        double total = 0;
        for (int i=0; i<GetNbranch(); i++)  {
            double shape = 1.0 / invshape;
            double scale = 1.0 / blmean.GetVal(i);
            total += GetVal(i).GetLogProb(shape,scale);
        }
        return total;
    }
};

/**
 * \brief An array of IID gamma random variables
 *
 * One should be careful about the fact that the shape and scale parameters are given by copy (not by ref) to the array.
 * Thus, each time the shape and scale parameters are modified during the MCMC,
 * the new values should be given to the array (using the SetShape and SetScale methods).
 */

class IIDGamma: public SimpleArray<double>	{

	public: 

    //! constructor specifies array size and initial shape and scale parameters
	IIDGamma(int insize, double inshape, double inscale) : SimpleArray<double>(insize), shape(inshape), scale(inscale)	{
		Sample();
	}

	~IIDGamma() {}

    //! return shape parameter
	double GetShape() const {return shape;}

    //! return scale parameter
	double GetScale() const {return scale;}

    //! set shape parameter to a new value
	void SetShape(double inshape)	{
		shape = inshape;
	}

    //! set scale parameter to a new value
	void SetScale(double inscale)	{
		scale = inscale;
	}

    //! sample all entries, given current shape and scale params
	void Sample()	{
		for (int i=0; i<GetSize(); i++)	{
			(*this)[i] = Random::GammaSample(shape,scale);
		}
	}

    //! resample all entries, given current shape and scale parameters and given an array of Poisson sufficient statistics of same size
	void GibbsResample(const Selector<PoissonSuffStat>& suffstatarray)	{
		for (int i=0; i<GetSize(); i++)	{
			const PoissonSuffStat& suffstat = suffstatarray.GetVal(i);
			(*this)[i] = Random::GammaSample(shape + suffstat.GetCount(), scale + suffstat.GetBeta());
		}
	}

    //! return total log prob (summed over the array) given current shape and scale params
	double GetLogProb()	const {
		double total = 0;
		for (int i=0; i<GetSize(); i++)	{
			total += GetLogProb(i);
		}
		return total;
	}

    //! return log prob of one specific entry
	double GetLogProb(int index) const {
        return Random::logGammaDensity(GetVal(index),shape,scale);
	}

    //! \brief given a Poisson suffstat S, calculates postprob[i] propto weight[i] * p(S | (*this)[i]), for i=0..GetSize()-1.
    //!
    //! interpretation: postprob[i] = posterior probability that the data summarized by S have been produced by a process of rate (*this)[i], for i=0..GetSize()-1
    void GetAllocPostProb(const PoissonSuffStat& suffstat, const vector<double>& weight, vector<double>& postprob) const {

        double max = 0;
        for (int i=0; i<GetSize(); i++) {
            double tmp = suffstat.GetLogProb(GetVal(i));
            postprob[i] = tmp;
            if ((!i) || (max < tmp))    {
                max = tmp;
            }
        }

        double total = 0;
        for (int i=0; i<GetSize(); i++) {
            postprob[i] = weight[i] * exp(postprob[i] - max);
            total += postprob[i];
        }

        for (int i=0; i<GetSize(); i++) {
            postprob[i] /= total;
        }

    }

    //! add all entries to a GammaSuffStat
	void AddSuffStat(GammaSuffStat& suffstat) const {
		for (int i=0; i<GetSize(); i++)	{
			suffstat.AddSuffStat(GetVal(i),log(GetVal(i)));
		}
	}

    //! add all entries for which occupancy[i] != 0 to a GammaSuffStat
	void AddSuffStat(GammaSuffStat& suffstat, const Selector<int>& occupancy) const	{
		for (int i=0; i<GetSize(); i++)	{
            if (occupancy.GetVal(i))   {
                suffstat.AddSuffStat(GetVal(i),log(GetVal(i)));
            }
		}
	}

    //! resample all entries for which occupancy[i] == 0 from the prior (from a Gamma(shape,scale))
    void PriorResample(const Selector<int>& occupancy)	{
		for (int i=0; i<GetSize(); i++)	{
            if (! occupancy.GetVal(i)) {
                (*this)[i] = Random::GammaSample(shape,scale);
            }
		}
    }

    //! add all entries for which poswarray[i] > 0 to a GammaSuffStat
	void AddSuffStat(GammaSuffStat& suffstat, const Selector<double>& poswarray) const {
		for (int i=0; i<GetSize(); i++)	{
            if (poswarray.GetVal(i)) {
                suffstat.AddSuffStat(GetVal(i),log(GetVal(i)));
            }
		}
	}

    //! resample all entries for which poswarray[i] == 0 from the prior (from a Gamma(shape,scale))
    void PriorResample(const Selector<double>& poswarray) {
        cerr << "in prior resample\n";
        cerr << "check that\n";
        exit(1);
		for (int i=0; i<GetSize(); i++)	{
            if (poswarray.GetVal(i)) {
                (*this)[i] = Random::GammaSample(shape,scale);
            }
		}
    }

    //! get mean over the array
    double GetMean() const {
        double m1 = 0;
        for (int i=0; i<GetSize(); i++) {
            m1 += GetVal(i);
        }
        m1 /= GetSize();
        return m1;
    }

    //! get variance over the array
    double GetVar() const {
        double m1 = 0;
        double m2 = 0;
        for (int i=0; i<GetSize(); i++) {
            m1 += GetVal(i);
            m2 += GetVal(i) * GetVal(i);
        }
        m1 /= GetSize();
        m2 /= GetSize();
        m2 -= m1*m1;
        return m2;
    }

	protected:
	double shape;
	double scale;
};

/**
 * \brief A tree-structured branch-wise array of iid Gamma variables (tree-stuctured version of IIDGamma)
 *
 * One should be careful about the fact that the shape and scale parameters are given by copy (not by ref) to the array.
 * Thus, each time the shape and scale parameters are modified during the MCMC,
 * the new values should be given to the array (using the SetShape and SetScale methods).
 */

class BranchIIDGamma: public SimpleBranchArray<double>	{

	public: 

	BranchIIDGamma(const Tree& intree, double inshape, double inscale) : SimpleBranchArray<double>(intree), shape(inshape), scale(inscale)	{
		Sample();
	}

	~BranchIIDGamma() {}

	double GetShape() const {return shape;}
	double GetScale() const {return scale;}

	void SetShape(double inshape)	{
		shape = inshape;
	}

	void SetScale(double inscale)	{
		scale = inscale;
	}

    //! set all entries equal to inval
	void SetAllBranches(double inval)	{
		for (int i=0; i<GetNbranch(); i++)	{
			(*this)[i] = inval;
		}
	}

    //! sample all entries from prior
	void Sample()	{
		for (int i=0; i<GetNbranch(); i++)	{
			(*this)[i] = Random::GammaSample(shape,scale);
		}
	}

    //! resample all entries from posterior, conditional on BranchArray of PoissonSuffStat
	void GibbsResample(const PoissonSuffStatBranchArray& suffstatarray)	{
		for (int i=0; i<GetNbranch(); i++)	{
			const PoissonSuffStat& suffstat = suffstatarray.GetVal(i);
			(*this)[i] = Random::GammaSample(shape + suffstat.GetCount(), scale + suffstat.GetBeta());
		}
	}

    //! get total log prob summed over all branches
	double GetLogProb()	{
		double total = 0;
		for (int i=0; i<GetNbranch(); i++)	{
			total += GetLogProb(i);
		}
		return total;
	}

    //! get log prob for a given branch
	double GetLogProb(int index)	{
        return Random::logGammaDensity(GetVal(index),shape,scale);
	}

    //! get sum over all entries (name is rather specialized... could change..)
    double GetTotalLength() const {
        double m1 = 0;
        for (int i=0; i<GetNbranch(); i++) {
            m1 += GetVal(i);
        }
        return m1;
    }

    //! add all entries to a GammaSuffStat
	void AddSuffStat(GammaSuffStat& suffstat)	{
		for (int i=0; i<GetNbranch(); i++)	{
			suffstat.AddSuffStat((*this)[i],log((*this)[i]));
		}
	}

    //! get mean over the array
    double GetMean() const {
        double m1 = 0;
        for (int i=0; i<GetNbranch(); i++) {
            m1 += GetVal(i);
        }
        m1 /= GetNbranch();
        return m1;
    }

    //! get variance over the array
    double GetVar() const {
        double m1 = 0;
        double m2 = 0;
        for (int i=0; i<GetNbranch(); i++) {
            m1 += GetVal(i);
            m2 += GetVal(i) * GetVal(i);
        }
        m1 /= GetNbranch();
        m2 /= GetNbranch();
        m2 -= m1*m1;
        return m2;
    }

	protected:
	double shape;
	double scale;
};

/**
 * \brief A tree-structured branch-wise array of Gamma variables, with branch-specific means but same shape parameter
 *
 * One should be careful about the fact that the shape parameter is given by copy (not by ref) to the array.
 * Thus, each time the shape parameter is modified during the MCMC,
 * the new value should be given to the array (using the SetShape method).
 */
class GammaWhiteNoise: public SimpleBranchArray<double>	{

	public: 

	GammaWhiteNoise(const Tree& intree, const BranchSelector<double>& inblmean, double inshape) : SimpleBranchArray<double>(intree), blmean(inblmean), shape(inshape)	{
		Sample();
	}

	~GammaWhiteNoise() {}

	double GetShape() const {return shape;}

	void SetShape(double inshape)	{
		shape = inshape;
	}

    //! sample all entries from prior
	void Sample()	{
		for (int i=0; i<GetNbranch(); i++)	{
            double scale = shape / blmean.GetVal(i);
			(*this)[i] = Random::GammaSample(shape,scale);
		}
	}

    //! resample entries based on a BranchArray of PoissonSuffStat
	void GibbsResample(const PoissonSuffStatBranchArray& suffstatarray)	{
		for (int i=0; i<GetNbranch(); i++)	{
			const PoissonSuffStat& suffstat = suffstatarray.GetVal(i);
            double scale = shape / blmean.GetVal(i);
			(*this)[i] = Random::GammaSample(shape + suffstat.GetCount(), scale + suffstat.GetBeta());
		}
	}

    //! return total log prob summed over all entries
	double GetLogProb()	{
		double total = 0;
		for (int i=0; i<GetNbranch(); i++)	{
			total += GetLogProb(i);
		}
		return total;
	}

    //! return log prob for one entry
	double GetLogProb(int index) const {
        double scale = shape / blmean.GetVal(index);
        return Random::logGammaDensity(GetVal(index),shape,scale);
	}

    //! add entries to a BranchArray of GammaSuffStat (!! valid only if all blmean[i]'s are the same across all branches)
	void AddSuffStat(GammaSuffStat& suffstat)	{
		for (int i=0; i<GetNbranch(); i++)	{
			suffstat.AddSuffStat((*this)[i],log((*this)[i]));
		}
	}

    //! add entries to a BranchArray of GammaSuffStat
    void AddSuffStat(GammaSuffStatBranchArray& suffstatarray)   {
		for (int i=0; i<GetNbranch(); i++)	{
			suffstatarray[i].AddSuffStat((*this)[i],log((*this)[i]));
		}
    }

    double GetTotalLength() const {
        double m1 = 0;
        for (int i=0; i<GetNbranch(); i++) {
            m1 += GetVal(i);
        }
        return m1;
    }

	protected:
    const BranchSelector<double>& blmean;
	double shape;
};

/**
 * \brief An array of GammaWhiteNoise
 *
 * This is an Array of BranchArray's. It is used for implementing modulations of branch-lengths across genes,
 * while implementing some basic shrinkage.
 * Specifically, for a given branch j,
 * the branch-lengths across genes are gamma of mean blmean[j] and shape parameter uniform across all branches.
 */

class GammaWhiteNoiseArray : public Array<GammaWhiteNoise>    {

    public:

    //! constructor: parameterized by the number of genes, the tree, the means over branches and the shape parameter
    GammaWhiteNoiseArray(int inNgene, const Tree& intree, const BranchSelector<double>& inblmean, double inshape) : Ngene(inNgene), tree(intree), blmean(inblmean), shape(inshape), blarray(Ngene, (GammaWhiteNoise*) 0)    {
        
        for (int gene=0; gene<Ngene; gene++)    {
            blarray[gene] = new GammaWhiteNoise(tree,blmean,shape);
        }
    }

    ~GammaWhiteNoiseArray()  {
        for (int gene=0; gene<Ngene; gene++)    {
            delete blarray[gene];
        }
    }

    //! set the shape parameter (should be called whenever the shape parameter has changed during the MCMC)
    void SetShape(double inshape)   {
        shape = inshape;
        for (int gene=0; gene<Ngene; gene++)    {
            blarray[gene]->SetShape(shape);
        }
    }

    //! return total number of entries (number of genes)
    int GetSize() const {
        return Ngene;
    }

    //! return total number of genes
    int GetNgene() const {
        return Ngene;
    }

    //! return total number of branches of the underlying tree
    int GetNbranch() const {
        return blarray[0]->GetNbranch();
    }

    //! const access to the GammaWhiteNoise BranchArray for the given gene
    const GammaWhiteNoise& GetVal(int gene) const {
        return *blarray[gene];
    }

    //! non-const access to the GammaWhiteNoise BranchArray for the given gene
    GammaWhiteNoise& operator[](int gene)  {
        return *blarray[gene];
    }

    //! return mean tree length over genes
    double GetMeanLength() const {
        double tot = 0;
        for (int j=0; j<GetNbranch(); j++)   {
            double mean = 0;
            for (int g=0; g<GetNgene(); g++) {
                double tmp = blarray[g]->GetVal(j);
                mean += tmp;
            }
            mean /= GetNgene();
            tot += mean;
        }
        return tot;
    }

    //! return variance of tree length across genes
    double GetVarLength() const {
        double tot = 0;
        for (int j=0; j<GetNbranch(); j++)   {
            double mean = 0;
            double var = 0;
            for (int g=0; g<GetNgene(); g++) {
                double tmp = blarray[g]->GetVal(j);
                mean += tmp;
                var += tmp*tmp;
            }
            mean /= GetNgene();
            var /= GetNgene();
            var -= mean*mean;
            tot += var;
        }
        tot /= GetNbranch();
        return tot;
    }

    //! return total log prob (over all genes and over all branches)
    double GetLogProb() const {
        double total = 0;
        for (int gene=0; gene<GetNgene(); gene++)  {
            total += blarray[gene]->GetLogProb();
        }
        return total;
    }

    //! add all gamma variates to a BranchArray of GammaSuffStat (gene-wise sum, separately for each branch)
    void AddSuffStat(GammaSuffStatBranchArray& suffstatarray) const    {
        for (int gene=0; gene<GetNgene(); gene++)  {
            blarray[gene]->AddSuffStat(suffstatarray);
        }
    }

    private:

    int Ngene;
    const Tree& tree;
    const BranchSelector<double>& blmean;
	double shape;
    vector<GammaWhiteNoise*> blarray;
};

#endif
