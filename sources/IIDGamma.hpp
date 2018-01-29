
#ifndef IIDGAMMA_H
#define IIDGAMMA_H

#include "Array.hpp"
#include "BranchArray.hpp"
#include "Random.hpp"
#include "PoissonSuffStat.hpp"
#include "MPIBuffer.hpp"

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

    //! resample all entries for which occupancy[i] == 0 from the prior (from a Gamma(shape,scale))
    void PriorResample(const Selector<int>& occupancy)	{
		for (int i=0; i<GetSize(); i++)	{
            if (! occupancy.GetVal(i)) {
                (*this)[i] = Random::GammaSample(shape,scale);
            }
		}
    }

    //! resample all entries for which poswarray[i] == 0 from the prior (from a Gamma(shape,scale))
    void PriorResample(const Selector<double>& poswarray) {
		for (int i=0; i<GetSize(); i++)	{
            if (! poswarray.GetVal(i)) {
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
	// void GibbsResample(const PoissonSuffStatBranchArray& suffstatarray)	{
	void GibbsResample(const BranchArray<PoissonSuffStat>& suffstatarray)	{
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

    private:

    int Ngene;
    const Tree& tree;
    const BranchSelector<double>& blmean;
	double shape;
    vector<GammaWhiteNoise*> blarray;
};

#endif
