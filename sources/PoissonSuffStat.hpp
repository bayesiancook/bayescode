
#ifndef POISSONSUFFSTAT_H
#define POISSONSUFFSTAT_H

#include "SuffStat.hpp"
#include "Array.hpp"
#include "BranchArray.hpp"
#include "PhyloProcess.hpp"
#include <cmath>
#include "MPIBuffer.hpp"

/**
 * \brief A Poisson-like sufficient statistic
 *
 * This sufficient statistic deals with all cases where the probability of the variable(s) of interest (say X),
 * as a function of some rate parameter r (positive real number),
 * can be written as P(X | r) \propto r^count exp(-beta*r),
 * for some suff stats count (integer) and beta (positive real number).
 * When several independent variables are entailed by the generic X mentioned above, the suff stats just have to be 
 * summed over all items (separately for count and beta).
 *
 * As an example, the probability of the substitution histories over all sites, for a given branch of the tree,
 * as a function of the branch length l, can be written as follows: p(S | l) = K l^count exp(-beta*l),
 * where count is the total number of substitution events along the branch (across all sites),
 * beta is the average rate away from the current state (averaged over all paths)
 * and K is some normalization constant (that may depend on other parameters than l).
 * Thus, conditional on the current substitution histories, 
 * the sufficient statistics count and beta can be first computed, 
 * and then MCMC moves on branch lengths can be done based on the knowledge of these two numbers.
 * Other examples include the probability of substitution histories as a function of omega = dN/dS
 * (in that case, the count suff stat is the number of non-synonymous substitutions only).
 * In these two cases, 
 * the suffstats have to be calculated based on the substitution histories (see BranchSitePath::AddLengthSuffStat)
 * and then summed over branches, sites, or any other more subtle pattern,
 * depending on the exact structure of the model.
 *
 * PoissonSuffStat implements this general idea:
 * collecting count and beta suff stats across all relevant items
 * and providing methods for calculating the probability, as a function of the rate parameter.
 */

class PoissonSuffStat : public SuffStat	{

	public:

	PoissonSuffStat() {}
	~PoissonSuffStat() {}

    //! set count and beta to 0
	void Clear()	{
		count = beta = 0;
	}

	void IncrementCount()	{
		count++;
	}

	void AddCount(int in)	{
		count += in;
	}

	void AddBeta(double in)	{
		beta += in;
	}

	void AddSuffStat(int incount, double inbeta)	{
		count += incount;
		beta += inbeta;
	}

    void Add(const PoissonSuffStat& from)   {
        count += from.GetCount();
        beta += from.GetBeta();
    }

    PoissonSuffStat& operator+=(const PoissonSuffStat& from)    {
        Add(from);
        return *this;
    }

    //! return size when put into an MPI buffer
    unsigned int GetMPISize() const {return 2;}

    //! put current value of count and beta into an MPI buffer
    void MPIPut(MPIBuffer& buffer) const    {
        buffer << beta << count;
    }

    //! get value from MPI buffer
    void MPIGet(const MPIBuffer& buffer)    {
        buffer >> beta >> count;
    }

    //! get a PoissonSuffStat from MPI buffer and then add it to this object
    void Add(const MPIBuffer& buffer)   {
        double temp;
        buffer >> temp;
        beta += temp;

        int tmp;
        buffer >> tmp;
        count += tmp;
    }

	int GetCount() const {
		return count;
	}

	double GetBeta() const {
		return beta;
	}

    //! return the log probability as a function of the rate: essentially count*log(rate) - beta*rate
	double GetLogProb(double rate) const {
		return count*log(rate) - beta*rate;
	}

    //! return the log of the marginal probability when the rate is from a gamma distribution
	double GetMarginalLogProb(double shape, double scale) const {
		return shape*log(scale) - Random::logGamma(shape) - (shape + count)*log(scale + beta) + Random::logGamma(shape + count);
	}

	protected:

	int count;
	double beta;
};

/**
 * \brief An array of Poisson sufficient statistics
 *
 * This class provides an interface for dealing with cases where
 * each item has a different rate parameter (e.g. rates r_i's or omega_i's across sites).
 * In that case, we need a distinct PoissonSuffStat for each item.
 * Each suff stat might still be summing over other subcases (e.g. across all branches for a given site)
 */

class PoissonSuffStatArray : public SimpleArray<PoissonSuffStat>	{

	public:

    //! constructor with array size
	PoissonSuffStatArray(int insize) : SimpleArray<PoissonSuffStat>(insize) {}
	~PoissonSuffStatArray() {}

    //! set suff stats to 0
	void Clear()	{
		for (int i=0; i<GetSize(); i++)	{
			(*this)[i].Clear();
		}
	}

    //! add path sufficient statistics for resampling branch lengths from PhyloProcess
    void AddRatePathSuffStat(const PhyloProcess& process)   {
        process.AddRateSuffStat(*this);
    }

    //! member-wise addition: essentially (*this)[i] += from[i], for i=0..GetSize()-1
    void Add(const PoissonSuffStatArray& from)    {
		for (int i=0; i<GetSize(); i++)	{
            (*this)[i].Add(from.GetVal(i));
        }
    }

    //! member-wise addition, operator version
    PoissonSuffStatArray& operator+=(const PoissonSuffStatArray& from)  {
        Add(from);
        return *this;
    }

    //! return size when put into an MPI buffer
    unsigned int GetMPISize() const {return 2 * GetSize();}

    //! put array into MPI buffer
    void MPIPut(MPIBuffer& buffer) const    {
		for (int i=0; i<GetSize(); i++)	{
            buffer << GetVal(i);
        }
    }

    //! get array from MPI buffer
    void MPIGet(const MPIBuffer& buffer)    {
		for (int i=0; i<GetSize(); i++)	{
            buffer >> (*this)[i];
        }
    }

    //! get array from MPI buffer and add it to this array (member-wise addition)
    void Add(const MPIBuffer& buffer)   {
		for (int i=0; i<GetSize(); i++)	{
            (*this)[i] += buffer;
        }
    }

    //! \brief get logprob, based on an array of rates (of same size)
    //!
    //! specifically, total log prob such as returned by the function is sum_i log p((*this)[i] | ratearray[i])
	double GetLogProb(const Array<double>& ratearray) const{
		double total = 0;
		for (int i=0; i<GetSize(); i++)	{
			total += GetVal(i).GetLogProb(ratearray.GetVal(i));
		}
		return total;
	}

    //! get marginal log prob, based on an array of rates that are iid from a gamma of given shape and scale parameters
	double GetMarginalLogProb(double shape, double scale)	const {
		double total = 0;
		/*
		for (int i=0; i<GetSize(); i++)	{
			total += GetVal(i).GetMarginalLogProb(shape,scale);
		}
		*/
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


/**
 * \brief A tree-structured branch-wise array of Poisson sufficient statistics
 *
 * This class provides an interface for dealing with cases where
 * each branch has a different rate (or, for that matter, length) parameter 
 * In that case, we need a distinct PoissonSuffStat for each item (for each branch).
 * Each suff stat might still be summing over subcases (e.g. across all sites for a given branch)
 */

class PoissonSuffStatBranchArray : public SimpleBranchArray<PoissonSuffStat>	{

	public:

    //! constructor parameterized by underlying treee
	PoissonSuffStatBranchArray(const Tree& intree) : SimpleBranchArray<PoissonSuffStat>(intree) {}
	~PoissonSuffStatBranchArray() {}

    //! set all suff stats to 0
	void Clear()	{
		for (int i=0; i<GetNbranch(); i++)	{
			(*this)[i].Clear();
		}
	}

    //! get total log prob, based on an array of branch-specific rates (or lengths)
	double GetLogProb(const BranchSelector<double>& ratearray) const{
		double total = 0;
		for (int i=0; i<GetNbranch(); i++)	{
			total += GetVal(i).GetLogProb(ratearray.GetVal(i));
		}
		return total;
	}

    //! member-wise addition (*this)[i] = from[i] for all i=0..GetNbranch()-1
    void Add(const PoissonSuffStatBranchArray& from)    {
		for (int i=0; i<GetNbranch(); i++)	{
            (*this)[i].Add(from.GetVal(i));
        }
    }

    //! member-wise addition: operator version
    PoissonSuffStatBranchArray& operator+=(const PoissonSuffStatBranchArray& from)  {
        Add(from);
        return *this;
    }

    //! add path sufficient statistics for resampling branch lengths from PhyloProcess
    void AddLengthPathSuffStat(const PhyloProcess& process)   {
        process.AddLengthSuffStat(*this);
    }

    //! return array size when put into an MPI buffer
    unsigned int GetMPISize() const {return 2 * GetNbranch();}

    //! put array into MPI buffer
    void MPIPut(MPIBuffer& buffer) const    {
		for (int i=0; i<GetNbranch(); i++)	{
            buffer << GetVal(i);
        }
    }

    //! get array from MPI buffer
    void MPIGet(const MPIBuffer& buffer)    {
		for (int i=0; i<GetNbranch(); i++)	{
            buffer >> (*this)[i];
        }
    }

    //! get an array from MPI buffer and then add it to this array
    void Add(const MPIBuffer& buffer)   {
		for (int i=0; i<GetNbranch(); i++)	{
            (*this)[i] += buffer;
        }
    }

    //! get total (summed) marginal log prob integrated over branch-specific rates (or lengths) iid from a gamma(shape,scale)
	double GetMarginalLogProb(double shape, double scale)	const {
		double total = 0;
		for (int i=0; i<GetNbranch(); i++)	{
			int count = GetVal(i).GetCount();
			double beta = GetVal(i).GetBeta();
			total += -(shape+count)*log(scale+beta) + Random::logGamma(shape+count);
		}
		total += GetNbranch() * (shape*log(scale) - Random::logGamma(shape));
		return total;
	}

    //! push two C-style arrays of counts and betas into this array (erase current values)
    void Push(int* count, double* beta) const   {
        for (int i=0; i<GetNbranch(); i++)  {
            count[i] = GetVal(i).GetCount();
            beta[i] = GetVal(i).GetBeta();
        }
    }

    //! push-and-add two C-style arrays of counts and betas into this array (member-wise addition)
    void Add(const int* count, const double* beta)  {
        for (int i=0; i<GetNbranch(); i++)  {
            (*this)[i].AddCount(count[i]);
            (*this)[i].AddBeta(beta[i]);
        }
    }

};

#endif
