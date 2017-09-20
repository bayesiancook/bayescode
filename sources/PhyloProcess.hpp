
#ifndef PHYLOPROCESS_H
#define PHYLOPROCESS_H

#include <map>
#include "BranchSitePath.hpp"
#include "SequenceAlignment.hpp"
#include "Chrono.hpp"
#include "SubMatrix.hpp"
#include "Tree.hpp"
#include "BranchSiteArray.hpp"
#include "BidimArray.hpp"
#include "BranchAllocationSystem.hpp"
 
// PhyloProcess is a dispatcher:
// its responsibility is to create a random branch/site path
// for each branch/site pair

class PhyloProcess	{

public:

	// generic constructor
	PhyloProcess(const Tree* intree, const SequenceAlignment* indata, const ConstBranchArray<double>* inbranchlength, const ConstArray<double>* insiterate, const ConstBranchSiteArray<SubMatrix>* insubmatrixarray, const ConstArray<SubMatrix>* inrootsubmatrixarray);

	// branch and site homogeneous
	PhyloProcess(const Tree* intree, const SequenceAlignment* indata, const ConstBranchArray<double>* inbranchlength, const ConstArray<double>* insiterate, const SubMatrix* insubmatrix);

	// branch homogeneous, site heterogeneous
	PhyloProcess(const Tree* intree, const SequenceAlignment* indata, const ConstBranchArray<double>* inbranchlength, const ConstArray<double>* insiterate, const ConstArray<SubMatrix>* insubmatrixarray);

	~PhyloProcess();

	// accessors

	double GetBranchLength(int branch) const {
		return branchlength->GetVal(branch);
	}

	double GetSiteRate(int site) const {
		if (! siterate)	{
			return 1.0;
		}
		return siterate->GetVal(site);
	}

	const SubMatrix& GetSubMatrix(int branch, int site) const {
		return submatrixarray->GetVal(branch,site);
	}

	const EVector& GetRootFreq(int site) const {
		return rootsubmatrixarray->GetVal(site).GetStationary();
	}

	double SiteLogLikelihood(int site) const;
	double FastSiteLogLikelihood(int site) const;

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

	// probability, pruning, sampling
	double GetLogProb() const;                        // likelihood Felsenstein 1981
	double GetFastLogProb() const;                            // likelihood Felsenstein 1981

	double Move(double fraction);
	void DrawSites(double fraction); // draw a fraction of sites which will be resampled
	void ResampleSub();  // clamped Nielsen
	void ResampleSub(int site);

	// basic building blocks for compiling suffstats...

    /*
	void AddRootSuffStat(int site, PathSuffStat& suffstat);
	void AddPathSuffStat(const Link* link, int site, PathSuffStat& suffstat);
	void AddLengthSuffStat(const Link* link, int site, PoissonSuffStat& suffstat);
    */
	// void AddPoissonSuffStat(const Link* link, int site, PoissonSuffStat& suffstat);

	// ... which are then used for looping over branches and sites
    //
    // void PrintMissing(const Link* from, int site);

	// homogeneous across sites and branches
	void AddPathSuffStat(PathSuffStat& suffstat);
	void RecursiveAddPathSuffStat(const Link* from, PathSuffStat& suffstat);
	void LocalAddPathSuffStat(const Link* from, PathSuffStat& suffstat);

	// heterogeneeous across sites, homogeneous across branches
	void AddPathSuffStat(Array<PathSuffStat>& suffstatarray);
	void RecursiveAddPathSuffStat(const Link* from, Array<PathSuffStat>& suffstatarray);
	void LocalAddPathSuffStat(const Link* from, Array<PathSuffStat>& suffstatarray);

	// heterogeneeous across sites, branches partitioned into conditions
	void AddPathSuffStat(BidimArray<PathSuffStat>& suffstatarray, const BranchAllocationSystem& branchalloc);
	void RecursiveAddPathSuffStat(const Link* from, BidimArray<PathSuffStat>& suffstatarray, const BranchAllocationSystem& branchalloc);
	void LocalAddPathSuffStat(const Link* from, BidimArray<PathSuffStat>& suffstatarray, int cond);

	// homogeneous across sites, heterogeneous across branches
	// void AddSuffStat(BranchArray<PathSuffStat>& branchsuffstatarray, PathSuffStat& rootsuffstat);

	// heterogeneous across sites and branches
	// void AddSuffStat(BranchSiteArray<PathSuffStat>& branchsitesuffstatarray);

	// homogeneous across sites
	void AddLengthSuffStat(BranchArray<PoissonSuffStat>& branchlengthsuffstatarray);
	void RecursiveAddLengthSuffStat(const Link* from, BranchArray<PoissonSuffStat>& branchlengthsuffstatarray);
	void LocalAddLengthSuffStat(const Link* from, PoissonSuffStat& branchlengthsuffstat);

	// homogeneous across branches
	// void AddPoissonSuffStat(Array<PoissonSuffStat>& poissonsuffstatarray);

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

	void Unfold();
	void Cleanup();

	protected:

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
	const ConstBranchArray<double>* branchlength;
	const ConstArray<double>* siterate;
	const ConstBranchSiteArray<SubMatrix>* submatrixarray; 
	const ConstArray<SubMatrix>* rootsubmatrixarray;
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

	private:

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
