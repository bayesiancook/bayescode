
#ifndef PHYLOPROCESS_H
#define PHYLOPROCESS_H

#include <map>
#include "BidimArray.hpp"
#include "BranchAllocationSystem.hpp"
#include "BranchSitePath.hpp"
#include "BranchSiteSelector.hpp"
#include "Chrono.hpp"
#include "SequenceAlignment.hpp"
#include "SubMatrix.hpp"
#include "Tree.hpp"

// PhyloProcess is a dispatcher:
// its responsibility is to create a random branch/site path
// for each branch/site pair

/**
 * \brief The core class of phylogenetic likelihood calculation and stochastic mapping of substitution histories
 *
 * PhyloProcess takes as an input a tree, a sequence alignment (data), a set of branch lengths, of site-specific rates and a
 * selector of substitution matrices across sites and branches.
 * It is then responsible for organizing all likelihood calculations by pruning, as stochastic mapping of substitution
 * histories.
 */

class PhyloProcess : public tc::Component {
  public:
    friend class PathSuffStat;
    friend class PathSuffStatArray;
    friend class PathSuffStatBidimArray;
    friend class PoissonSuffStatBranchArray;

    //! \brief generic constructor
    //!
    //! Constructor takes as parameters (pointers):
    //! - tree
    //! - sequence alignment
    //! - a BranchSelector of branch lengths
    //! - a site Selector of site rates (if pointer is null, then rates across sites are all equal to 1)
    //! - a BranchSiteSelector specifying which substitution matrix should be used for each branch site pair
    //! - a site Selector of substitution matrices, specifying which matrix should be used for getting the equilibrium
    //! frequencies at each site, from which to draw the root state
    PhyloProcess(const Tree* intree, const SequenceAlignment* indata, BranchSelector<double>* inbranchlength,
                 const Selector<double>* insiterate, const BranchSiteSelector<SubMatrix>* insubmatrixarray,
                 const Selector<SubMatrix>* inrootsubmatrixarray);

    //! \brief special (short-cut) constructor for branch-homogeneous and site-homogeneous model
    //!
    //! Compared to the generic constructor, this constructor takes a pointer to a single substitution matrix.
    //! The branch/site and site selectors of substitution matrices (submatrixarray and rootsubmatrixarray)
    //! are then internally allocated by PhyloProcess based on this matrix.
    //! If insiterate pointer is null, then rates across sites are all equal to 1.
    PhyloProcess(const Tree* intree, const SequenceAlignment* indata, BranchSelector<double>* inbranchlength,
                 const Selector<double>* insiterate, const SubMatrix* insubmatrix);

    // Component-compatible version (separated in two)
    PhyloProcess(const Tree* intree, const SequenceAlignment* indata)
        : tree(intree), data(indata), Nstate(data->GetNstate()), maxtrial(DEFAULTMAXTRIAL) {
        port("branchlength", &PhyloProcess::branchlength);
        port("siterate", &PhyloProcess::siterate);
        port("submatrix", &PhyloProcess::set_submatrix);
    }

    void set_submatrix(SubMatrix* insubmatrix) {
        submatrixarray = new BranchHomogeneousSiteHomogeneousSelector<SubMatrix>(*tree, GetNsite(), *insubmatrix);
        allocsubmatrixarray = true;
        rootsubmatrixarray = new HomogeneousSelector<SubMatrix>(GetNsite(), *insubmatrix);
        allocrootsubmatrixarray = true;
    }

    //! \brief special (short-cut) constructor for branch-homogeneous and site-heterogeneous model
    //!
    //! Compared to the generic constructor, this constructor takes a pointer to a (site) Selector<SubMatrix>*
    //! insubmatrixarray.
    //! The branch/site selector of substitution matrices (submatrixarray)
    //! is then internally allocated by PhyloProcess based on this matrix based on this array,
    //! while rootmatrixarray is set to insubmatrixarray.
    //! If insiterate pointer is null, then rates across sites are all equal to 1.
    PhyloProcess(const Tree* intree, const SequenceAlignment* indata, BranchSelector<double>* inbranchlength,
                 const Selector<double>* insiterate, const Selector<SubMatrix>* insubmatrixarray);

    ~PhyloProcess();

    //! return log likelihood (computed using the pruning algorithm, Felsenstein 1981)
    double GetLogLikelihood() const;

    //! return log likelihood for given site
    double SiteLogLikelihood(int site) const;

    //! stochastic sampling of substitution history under current parameter configuration
    void ResampleSub();

    //! stochastic sampling of substitution history for given random fraction of sites
    double Move(double fraction);

    //! create all data structures necessary for computation
    void Unfold();

    //! delete data structures
    void Cleanup();

  private:
    //! \brief const access to substitution history (BranchSitePath) for given node and given site
    //!
    //! Substitution histories are indexed by node (not by branch);
    //! root node also has a substitution history (starting state).
    const BranchSitePath* GetPath(const Node* node, int site) const {
        map<const Node*, BranchSitePath**>::const_iterator it = pathmap.find(node);
        if (it == pathmap.end()) {
            std::cerr << "error in phyloprocess::did not find entry for node in branch site path\n";
            exit(1);
        }
        const BranchSitePath* path = it->second[site];
        if (path == nullptr) {
            std::cerr << "error in phyloprocess::getpath: null path\n";
            exit(1);
        }
        return path;
    }

    double GetFastLogProb() const;
    double FastSiteLogLikelihood(int site) const;

    //! return branch length for given branch
    double GetBranchLength(int branch) const { return branchlength->GetVal(branch); }

    //! return site rate for given site (if no rates-across-sites array was given to phyloprocess, returns 1)
    double GetSiteRate(int site) const {
        if (!siterate) {
            return 1.0;
        }
        return siterate->GetVal(site);
    }

    //! return matrix that should be used
    const SubMatrix& GetSubMatrix(int branch, int site) const { return submatrixarray->GetVal(branch, site); }

    const EVector& GetRootFreq(int site) const { return rootsubmatrixarray->GetVal(site).GetStationary(); }

    const StateSpace* GetStateSpace() const { return data->GetStateSpace(); }
    const TaxonSet* GetTaxonSet() const { return data->GetTaxonSet(); }

    int GetNsite() const { return data->GetNsite(); }
    int GetNtaxa() const { return data->GetNtaxa(); }

    int GetNstate() const { return Nstate; }

    const SequenceAlignment* GetData() const { return data; }
    int GetData(int taxon, int site) const { return data->GetState(taxon, site); }

    const Tree* GetTree() const { return tree; }
    const Link* GetRoot() const { return GetTree()->GetRoot(); }

    int GetMaxTrial() const { return maxtrial; }
    void SetMaxTrial(int i) { maxtrial = i; }

    void SetData(const SequenceAlignment* indata);
    void ClampData() { clampdata = true; }
    void UnclampData() { clampdata = false; }

    int& GetState(const Node* node, int site) { return statemap[node][site]; }
    const int& GetState(const Node* node, int site) const { return statemap[node][site]; }

    void GetLeafData(SequenceAlignment* data);
    void RecursiveGetLeafData(const Link* from, SequenceAlignment* data);

    bool isDataCompatible(int taxon, int site, int state) const {
        return GetStateSpace()->isCompatible(GetData(taxon, site), state);
    }

    void DrawSites(double fraction);  // draw a fraction of sites which will be resampled
    void ResampleSub(int site);

    //! compute path sufficient statistics across all sites and branches and add them to suffstat (site-branch-homogeneous
    //! model)
    void AddPathSuffStat(PathSuffStat& suffstat) const;
    //! compute path sufficient statistics across all sites and branches and add them to suffstatarray (site-heterogeneous
    //! branch-homogeneous model)
    void AddPathSuffStat(Array<PathSuffStat>& suffstatarray) const;
    //! compute path sufficient statistics across all sites and branches and add them to bidim suffstatarray (branches
    //! partitioned into conditions, see DiffSelModel)
    //! heterogeneeous across sites, branches partitioned into conditions
    void AddPathSuffStat(BidimArray<PathSuffStat>& suffstatarray, const BranchAllocationSystem& branchalloc) const;

    //! compute path sufficient statistics for resampling branch lengths add them to branchlengthpathsuffstatarray
    void AddLengthSuffStat(BranchArray<PoissonSuffStat>& branchlengthpathsuffstatarray) const;

    void RecursiveAddPathSuffStat(const Link* from, PathSuffStat& suffstat) const;
    void LocalAddPathSuffStat(const Link* from, PathSuffStat& suffstat) const;

    void RecursiveAddPathSuffStat(const Link* from, Array<PathSuffStat>& suffstatarray) const;
    void LocalAddPathSuffStat(const Link* from, Array<PathSuffStat>& suffstatarray) const;

    void RecursiveAddPathSuffStat(const Link* from, BidimArray<PathSuffStat>& suffstatarray,
                                  const BranchAllocationSystem& branchalloc) const;
    void LocalAddPathSuffStat(const Link* from, BidimArray<PathSuffStat>& suffstatarray, int cond) const;

    void RecursiveAddLengthSuffStat(const Link* from, BranchArray<PoissonSuffStat>& branchlengthpathsuffstatarray) const;
    void LocalAddLengthSuffStat(const Link* from, PoissonSuffStat& branchlengthsuffstat) const;

    void PostPredSample(bool rootprior = false);  // unclamped Nielsen
    void PostPredSample(int site, bool rootprior = false);
    // rootprior == true : root state drawn from stationary probability of the
    // process
    // rootprior == false: root state drawn from posterior distribution

    // various accessors

    bool isMissing(const Node* node, int site) const { return false; }

    bool isMissing(const Link* link, int site) const {
        return false;
        // return (missingmap[link->GetNode()][site] || missingmap[link->Out()->GetNode()][site]);
    }

    void CreateMissingMap();
    void DeleteMissingMap();
    void RecursiveCreateMissingMap(const Link* from);
    void FillMissingMap();
    void BackwardFillMissingMap(const Link* from);
    void ForwardFillMissingMap(const Link* from, const Link* up);

    double* GetCondLikelihood(const Link* from) const {
        map<const Link*, double*>::const_iterator i = condlmap.find(from);
        if (i == condlmap.end()) {
            cerr << "error in PhyloProcess::GetCondLikelihood\n";
            exit(1);
        }
        return i->second;
    }

    double GetPruningTime() const { return pruningchrono.GetTime(); }
    double GetResampleTime() const { return resamplechrono.GetTime(); }

    void RecursiveCreate(const Link* from);
    void RecursiveDelete(const Link* from);

    void RecursiveCreateTBL(const Link* from);
    void RecursiveDeleteTBL(const Link* from);

    void Pruning(const Link* from, int site) const;
    void ResampleSub(const Link* from, int site);
    void ResampleState();
    void ResampleState(int site);
    void PruningAncestral(const Link* from, int site);
    void PriorSample(const Link* from, int site, bool rootprior);
    void PriorSample();
    void RootPosteriorDraw(int site);

    // borrowed from phylobayes
    // where should that be?
    BranchSitePath* SamplePath(int stateup, int statedown, double time, double rate, const SubMatrix& matrix);
    BranchSitePath* SampleRootPath(int rootstate);
    BranchSitePath* ResampleAcceptReject(int maxtrial, int stateup, int statedown, double rate, double totaltime,
                                         const SubMatrix& matrix);
    BranchSitePath* ResampleUniformized(int stateup, int statedown, double rate, double totaltime, const SubMatrix& matrix);

    const Tree* tree;
    const SequenceAlignment* data;
    BranchSelector<double>* branchlength;
    const Selector<double>* siterate;
    const BranchSiteSelector<SubMatrix>* submatrixarray;
    const Selector<SubMatrix>* rootsubmatrixarray;
    bool allocsubmatrixarray;
    bool allocrootsubmatrixarray;

    int* sitearray;
    mutable double* sitelnL;

    int Nstate;

    bool clampdata;

    BranchSitePath* GetPath(const Node* node, int site) {
        if (pathmap[node][site] == nullptr) {
            std::cerr << "error in phyloprocess::getpath: null path\n";
            exit(1);
        }
        return pathmap[node][site];
    }

    mutable std::map<const Link*, double*> condlmap;
    mutable std::map<const Node*, BranchSitePath**> pathmap;
    mutable std::map<const Node*, int*> statemap;
    // std::map<const Node *, int> totmissingmap;

    int** missingmap;

    int maxtrial;
    static const int unknown = -1;

    static const int DEFAULTMAXTRIAL = 100;

    mutable Chrono pruningchrono;
    mutable Chrono resamplechrono;
};

#endif  // PHYLOPROCESS_H
