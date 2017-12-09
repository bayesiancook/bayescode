
#ifndef PHYLOPROCESS_H
#define PHYLOPROCESS_H

#include <map>
#include "BranchSitePath.hpp"
#include "SequenceAlignment.hpp"
#include "Chrono.hpp"
#include "SubMatrix.hpp"
#include "Tree.hpp"
#include "BranchSiteSelector.hpp"
#include "BidimArray.hpp"
#include "BranchAllocationSystem.hpp"
 
// PhyloProcess is a dispatcher:
// its responsibility is to create a random branch/site path
// for each branch/site pair

/**
 * \brief The core class of phylogenetic likelihood calculation and stochastic mapping of substitution histories
 *
 * PhyloProcess takes as an input a tree, a sequence alignment (data), a set of branch lengths, of site-specific rates and a selector of substitution matrices across sites and branches.
 * It is then responsible for organizing all likelihood calculations by pruning, as stochastic mapping of substitution histories.
 */

class PhyloProcess	{

public:

	//! \brief generic constructor
	PhyloProcess(const Tree* intree, const SequenceAlignment* indata, const BranchSelector<double>* inbranchlength, const Selector<double>* insiterate, const BranchSiteSelector<SubMatrix>* insubmatrixarray, const Selector<SubMatrix>* inrootsubmatrixarray);

	//! \brief special (short-cut) constructor for branch-homogeneous and site-homogeneous model
	PhyloProcess(const Tree* intree, const SequenceAlignment* indata, const BranchSelector<double>* inbranchlength, const Selector<double>* insiterate, const SubMatrix* insubmatrix);

	//! \brief special (short-cut) constructor for branch-homogeneous and site-heterogeneous model
	PhyloProcess(const Tree* intree, const SequenceAlignment* indata, const BranchSelector<double>* inbranchlength, const Selector<double>* insiterate, const Selector<SubMatrix>* insubmatrixarray);

	~PhyloProcess();

    //! return log likelihood (integrated over all substitution histories, by pruning, Felsenstein 1981)
	double GetLogProb() const;              

	double GetFastLogProb() const;

	double SiteLogLikelihood(int site) const;
	double FastSiteLogLikelihood(int site) const;

    //! perform stochastic mapping of substitution events under current parameter configuration
	void ResampleSub();

    //! perform stochastic mapping of substitution events for a given (random) fraction of all sites
	double Move(double fraction);

    //! unfold all data structures, do the likelihood computation and stochastic mapping of substitution histories
	void Unfold();

    //! delete data structures
	void Cleanup();

	//! homogeneous across sites and branches
	void AddPathSuffStat(PathSuffStat& suffstat);
	//! heterogeneeous across sites, homogeneous across branches
	void AddPathSuffStat(Array<PathSuffStat>& suffstatarray);
	//! heterogeneeous across sites, branches partitioned into conditions
	void AddPathSuffStat(BidimArray<PathSuffStat>& suffstatarray, const BranchAllocationSystem& branchalloc);
	//! homogeneous across sites
	void AddLengthSuffStat(BranchArray<PoissonSuffStat>& branchlengthsuffstatarray);

	// homogeneous across branches
	// void AddPoissonSuffStat(Array<PoissonSuffStat>& poissonsuffstatarray);
	// homogeneous across sites, heterogeneous across branches
	// void AddSuffStat(BranchArray<PathSuffStat>& branchsuffstatarray, PathSuffStat& rootsuffstat);
	// heterogeneous across sites and branches
	// void AddSuffStat(BranchSiteArray<PathSuffStat>& branchsitesuffstatarray);


    private:

    //! return branch length for given branch
	double GetBranchLength(int branch) const {
		return branchlength->GetVal(branch);
	}

    //! return site rate for given site (if no rates-across-sites array was given to phyloprocess, returns 1)
	double GetSiteRate(int site) const {
		if (! siterate)	{
			return 1.0;
		}
		return siterate->GetVal(site);
	}

    //! return matrix that should be used 
	const SubMatrix& GetSubMatrix(int branch, int site) const {
		return submatrixarray->GetVal(branch,site);
	}

	const EVector& GetRootFreq(int site) const {
		return rootsubmatrixarray->GetVal(site).GetStationary();
	}

	const StateSpace* GetStateSpace() const { return data->GetStateSpace(); }
	const TaxonSet* GetTaxonSet() const { return data->GetTaxonSet(); }

	int GetNsite() const { return data->GetNsite(); }
	int GetNtaxa() const { return data->GetNtaxa(); }

	int GetNstate() const { return Nstate;}

	const SequenceAlignment *GetData() const { return data; }
	int GetData(int taxon, int site) const { return data->GetState(taxon, site); }

	const Tree* GetTree() const { return tree;}
	const Link* GetRoot() const { return GetTree()->GetRoot(); }

	int GetMaxTrial() const { return maxtrial; }
	void SetMaxTrial(int i) { maxtrial = i; }

	void SetData(const SequenceAlignment *indata);
	void ClampData() { clampdata = true; }
	void UnclampData() { clampdata = false; }

	int& GetState(const Node *node, int site) { return statemap[node][site]; }

	void GetLeafData(SequenceAlignment *data);
	void RecursiveGetLeafData(const Link *from, SequenceAlignment *data);

	bool isDataCompatible(int taxon, int site, int state) const {
		return GetStateSpace()->isCompatible(GetData(taxon, site), state);
	}

	void DrawSites(double fraction); // draw a fraction of sites which will be resampled
	void ResampleSub(int site);

	void RecursiveAddPathSuffStat(const Link* from, PathSuffStat& suffstat);
	void LocalAddPathSuffStat(const Link* from, PathSuffStat& suffstat);

	void RecursiveAddPathSuffStat(const Link* from, Array<PathSuffStat>& suffstatarray);
	void LocalAddPathSuffStat(const Link* from, Array<PathSuffStat>& suffstatarray);

	void RecursiveAddPathSuffStat(const Link* from, BidimArray<PathSuffStat>& suffstatarray, const BranchAllocationSystem& branchalloc);
	void LocalAddPathSuffStat(const Link* from, BidimArray<PathSuffStat>& suffstatarray, int cond);

	void RecursiveAddLengthSuffStat(const Link* from, BranchArray<PoissonSuffStat>& branchlengthsuffstatarray);
	void LocalAddLengthSuffStat(const Link* from, PoissonSuffStat& branchlengthsuffstat);

	void PostPredSample(bool rootprior = false);  // unclamped Nielsen
	void PostPredSample(int site, bool rootprior = false);
	// rootprior == true : root state drawn from stationary probability of the
	// process
	// rootprior == false: root state drawn from posterior distribution

	// various accessors

	bool isMissing(const Node *node, int site) const { 
        return false; 
    }

	bool isMissing(const Link *link, int site) const {
		return false;
		// return (missingmap[link->GetNode()][site] || missingmap[link->Out()->GetNode()][site]);
	}

	void CreateMissingMap();
	void DeleteMissingMap();
	void RecursiveCreateMissingMap(const Link *from);
	void FillMissingMap();
	void BackwardFillMissingMap(const Link *from);
	void ForwardFillMissingMap(const Link *from, const Link* up);

	double* GetCondLikelihood(const Link *from) const {
        map<const Link*, double*>::const_iterator i = condlmap.find(from); 
        if (i == condlmap.end())    {
            cerr << "error in PhyloProcess::GetCondLikelihood\n";
            exit(1);
        }
        return i->second;
    }

	double GetPruningTime() const { return pruningchrono.GetTime(); }
	double GetResampleTime() const { return resamplechrono.GetTime(); }

	void RecursiveCreate(const Link *from);
	void RecursiveDelete(const Link *from);

	void RecursiveCreateTBL(const Link *from);
	void RecursiveDeleteTBL(const Link *from);

	void Pruning(const Link *from, int site) const ;
	void ResampleSub(const Link *from, int site);
	void ResampleState();
	void ResampleState(int site);
	void PruningAncestral(const Link *from, int site);
	void PriorSample(const Link *from, int site, bool rootprior);
	void PriorSample();
	void RootPosteriorDraw(int site);

	// borrowed from phylobayes
	// where should that be?
	BranchSitePath* SamplePath(int stateup, int statedown, double time, double rate, const SubMatrix& matrix);
	BranchSitePath* SampleRootPath(int rootstate);
	BranchSitePath* ResampleAcceptReject(int maxtrial, int stateup, int statedown, double rate, double totaltime, const SubMatrix& matrix);
	BranchSitePath* ResampleUniformized(int stateup, int statedown, double rate, double totaltime, const SubMatrix& matrix);

	const Tree *tree;
	const SequenceAlignment *data;
	const BranchSelector<double>* branchlength;
	const Selector<double>* siterate;
	const BranchSiteSelector<SubMatrix>* submatrixarray; 
	const Selector<SubMatrix>* rootsubmatrixarray;
	bool allocsubmatrixarray;
	bool allocrootsubmatrixarray;

	int *sitearray;
	mutable double *sitelnL;

	int Nstate;

	bool clampdata;

	BranchSitePath* GetPath(const Node* node, int site)	{
	    if (pathmap[node][site] == nullptr) {
		std::cerr << "error in phyloprocess::getpath: null path\n";
		exit(1);
	    }
	    return pathmap[node][site];
	}

	mutable std::map<const Link *, double *> condlmap;
	std::map<const Node*, BranchSitePath **> pathmap;
	std::map<const Node *, int *> statemap;
	// std::map<const Node *, int> totmissingmap;

    int** missingmap;

	int maxtrial;
	static const int unknown = -1;

	static const int DEFAULTMAXTRIAL = 100;

	mutable Chrono pruningchrono;
	mutable Chrono resamplechrono;
};

#endif  // PHYLOPROCESS_H
