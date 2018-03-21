/*Copyright or © or Copr. Centre National de la Recherche Scientifique (CNRS) (2017-06-14).
Contributors:
* Nicolas LARTILLOT - nicolas.lartillot@univ-lyon1.fr

This software is a computer program whose purpose is to detect convergent evolution using Bayesian
phylogenetic codon models.

This software is governed by the CeCILL-C license under French law and abiding by the rules of
distribution of free software. You can use, modify and/ or redistribute the software under the terms
of the CeCILL-C license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and rights to copy, modify and redistribute
granted by the license, users are provided only with a limited warranty and the software's author,
the holder of the economic rights, and the successive licensors have only limited liability.

In this respect, the user's attention is drawn to the risks associated with loading, using,
modifying and/or developing or reproducing the software by the user in light of its specific status
of free software, that may mean that it is complicated to manipulate, and that also therefore means
that it is reserved for developers and experienced professionals having in-depth computer knowledge.
Users are therefore encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or data to be ensured and,
more generally, to use and operate it in the same conditions as regards security.

The fact that you are presently reading this means that you have had knowledge of the CeCILL-C
license and that you accept its terms.*/


#include "CodonSequenceAlignment.hpp"
#include "GTRSubMatrix.hpp"
#include "PhyloProcess.hpp"
#include "ProbModel.hpp"
#include "Tree.hpp"
#include "IIDMultiGamma.hpp"
#include "IIDMultiBernoulli.hpp"
#include "MultiGammaSuffStat.hpp"
#include "DiffSelSparseFitnessArray.hpp"
#include "BranchAllocationSystem.hpp"
#include "AADiffSelCodonMatrixBidimArray.hpp"
#include "SubMatrixSelector.hpp"
#include "IIDGamma.hpp"
#include "GammaSuffStat.hpp"
#include "IIDDirichlet.hpp"
#include "PathSuffStat.hpp"
#include "IIDProfileMask.hpp"

/**
 * \brief A doubly-sparse version of the differential selection model (see DiffSelModel)
 *
 * This is a codon model, based on the mutation-selection formalism,
 * such that the fitness landscape (defined by site-specific amimo-acid fitness profiles)
 * is modulated across a set of alternative 'conditions', or environments, specified over the tree.
 * The model is meant to identify positions displaying significant patterns of directional convergent selection
 * (i.e. sites that tend to substitute their amino-acid state in a consistent manner, upon repeated transitions into a specific ecological condition).
 *
 */

class DiffSelDoublySparseModel : public ProbModel {

    // -----
    // model selectors
    // -----

    int codonmodel;

    // 0: free wo shrinkage
    // 1: free with shrinkage
    // 2: shared across genes
    // 3: fixed

    int blmode;
    int nucmode;
    int fitnessshapemode;
    int fitnesscentermode;
    int maskmode;
    int maskepsilonmode;

    int withtoggle;

    // -----
    // external parameters
    // -----

    Tree* tree;
    FileSequenceAlignment* data;
    const TaxonSet* taxonset;
    CodonSequenceAlignment* codondata;

    // number of sites
    int Nsite;
    int Ntaxa;
    int Nbranch;

    // number of diff sel categories
    int Ncond;

    // number of levels of the model
    // with 2 levels, structure of the model is as follows:
    // baseline (condition 0)
    // baseline  || fitness1 (for condition 1)
    // baseline || fitness1  || fitnessk  (for condition k=2..Ncond)
    int Nlevel;

    // which branch is under which condition
    BranchAllocationSystem* branchalloc;

    // -----
    //  model structure
    // -----

    // branch lengths 
	double lambda;
	BranchIIDGamma* blhypermean;
    double blhyperinvshape;
    GammaWhiteNoise* branchlength;
	PoissonSuffStatBranchArray* lengthpathsuffstatarray;
	GammaSuffStat hyperlengthsuffstat;

    // nucleotide exchange rates and equilibrium frequencies (stationary probabilities)
    // hyperparameters
    vector<double> nucrelratehypercenter;
    double nucrelratehyperinvconc;
    vector<double> nucstathypercenter;
    double nucstathyperinvconc;
    // parameters
	std::vector<double> nucrelrate;
	std::vector<double> nucstat;
	GTRSubMatrix* nucmatrix;

    double fitnessshape;
    vector<double> fitnesscenter;
    BidimIIDMultiGamma* fitness;

    double maskprob;
    double maskepsilon;
    IIDProfileMask* sitemaskarray;

    // shiftprob (across conditions):
    // either Beta(shiftprobhypermean,shiftprobhyperinvconc), estimated across genes
    // or mixture 1-pi * 0 + pi * Beta: and this, for each condition separately
    vector<double> pi;
    vector<double> shiftprobhypermean;
    vector<double> shiftprobhyperinvconc;
    vector<double> shiftprob;

    BidimIIDMultiBernoulli* toggle;

    // fitness profiles (combinations of baseline and delta)
    // across conditions and across sites
    DiffSelDoublySparseFitnessArray* fitnessprofile;

    // codon substitution matrices
    // across conditions and sites
    AADiffSelCodonMatrixBidimArray* condsubmatrixarray;

    // branch- and site-substitution matrices (for phyloprocess)
    SubMatrixSelector* submatrixarray;
    // and for root (condition 0)
    RootSubMatrixSelector* rootsubmatrixarray;

    // phyloprocess
    PhyloProcess* phyloprocess;

    // suff stats

    // path suff stats across conditions and sites
    PathSuffStatBidimArray* suffstatarray;

    MultiGammaSuffStat hyperfitnesssuffstat;

    int gammanullcount;

  public:

    //! \brief constructor
    //!
    //! parameters:
    //! - datafile: name of file containing codon sequence alignment
    //! - treefile: name of file containing tree topology (and branch conditions, such as specified by branch names)
    //! - Ncond: number of conditions (K)
    //! - Nlevel: number of levels (if Nlevel == 1: each condition is defined w.r.t. condition 0; if Nlevel == 2, condition 1 is defined w.r.t. condition 0, and condition 2..K-1 all defined w.r.t. condition 1)
    //! - codonmodel: type of codon substitution model (1: canonical mutation-selection model, 0: square-root model, see Parto and Lartillot, 2017)
    DiffSelDoublySparseModel(const std::string& datafile, const std::string& treefile, int inNcond, int inNlevel, int incodonmodel, double inepsilon, double inshape) : hyperfitnesssuffstat(Naa) {

        withtoggle = 0;
        fitnesscentermode = 3;
        fitnessshapemode = 3;
        fitnessshape = 20.0;

        codonmodel = incodonmodel;

        blmode = 0;
        nucmode = 0;

        if (inshape > 0)    {
            fitnessshapemode = 3;
            fitnessshape = inshape;
        }
        else    {
            fitnessshapemode = 0;
            fitnessshape = 20.0;
        }

        if (inepsilon == 1)   {
            maskepsilon = 1;
            maskmode = 3;
            maskepsilonmode = 3;
        }
        else if (inepsilon >= 0) {
            maskepsilon = inepsilon;
            maskepsilonmode = 3;
            maskmode = 0;
        }
        else    {
            maskepsilonmode = 0;
            maskmode = 0;
        }

        Ncond = inNcond;
        Nlevel = inNlevel;
        if (Ncond <= 2) {
            Nlevel = 1;
        }

        ReadFiles(datafile, treefile);

        // specifies which condition for which branch
        branchalloc = new BranchAllocationSystem(*tree,Ncond);
    }

    // DiffSelDoublySparseModel(const DiffSelDoublySparseModel&) = delete;

    ~DiffSelDoublySparseModel() {}

    //! read files (and read out the distribution of conditions across branches, based on the tree read from treefile)
    void ReadFiles(string datafile, string treefile) {
        // nucleotide sequence alignment
        data = new FileSequenceAlignment(datafile);

        // translated into codon sequence alignment
        codondata = new CodonSequenceAlignment(data, true);

        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();

        std::cerr << "-- Number of sites: " << Nsite << std::endl;

        taxonset = codondata->GetTaxonSet();

        // get tree from file (newick format)
        tree = new Tree(treefile);

        // check whether tree and data fits together
        tree->RegisterWith(taxonset);

        // traversal of the tree, so as to number links, branches and nodes
        // convention is: branches start at 1 (branch number 0 is the null branch behind the root)
        // nodes start at 0 (for the root), and nodes 1..Ntaxa are tip nodes (corresponding to taxa
        // in sequence alignment)
        tree->SetIndices();
        Nbranch = tree->GetNbranch();
    }

    //! allocate the model (data structures)
    void Allocate() {

        // ----------
        // construction of the model
        // ----------

        // allocating data structures and sampling initial configuration

        // branch lengths
		lambda = 10;
        blhypermean = new BranchIIDGamma(*tree,1.0,lambda);
        blhypermean->SetAllBranches(1.0 / lambda);
        blhyperinvshape = 1.0;
        branchlength = new GammaWhiteNoise(*tree,*blhypermean,1.0/blhyperinvshape);
        lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);

        nucrelratehypercenter.assign(Nrr,1.0/Nrr);
        nucrelratehyperinvconc = 1.0 / Nrr;

        nucstathypercenter.assign(Nnuc,1.0/Nnuc);
        nucstathyperinvconc = 1.0 / Nnuc;

        // nucleotide mutation matrix
		nucrelrate.assign(Nrr,0);
        Random::DirichletSample(nucrelrate,vector<double>(Nrr,1.0/Nrr),((double) Nrr));
		nucstat.assign(Nnuc,0);
        Random::DirichletSample(nucstat,vector<double>(Nnuc,1.0/Nnuc),((double) Nnuc));
		nucmatrix = new GTRSubMatrix(Nnuc,nucrelrate,nucstat,true);

        // fitness parameters: IID Gamma, across all conditions, sites, and amino-acids
        // those are not the final fitness values (depends on the system of masks and toggles, specified below)
        fitnesscenter.assign(Naa,1.0/Naa);
        fitness = new BidimIIDMultiGamma(Ncond,Nsite,Naa,fitnessshape,fitnesscenter);

        // profiles across sites are masked:
        // each site has a 20-dim mask, iid bernoulli of prob maskprob, conditional on at least one entry being 1
        // if mask[a] == 0 for amino-acid a, then its fitness is equal to maskepsilon, across all conditions, for that site
        maskprob = 0.1;
        sitemaskarray = new IIDProfileMask(Nsite,Naa,maskprob);

        // hyperparameters for the system of toggles
        // shiftprob is a vector of Ncond-1 probabilities;
        // for each non-baseline condition, k=1..Ncond:
        //     - with probability 1-pi[k-1], shiftprob[k-1] = 0 
        //     - with probability pi[k-1]  , shiftprob[k-1] ~ Beta(shiftprobhypermean, shiftprobhyperinvconc)
        pi.assign(Ncond-1,0.1);
        shiftprobhypermean.assign(Ncond-1,0.1);
        shiftprobhyperinvconc.assign(Ncond-1,0.5);
        shiftprob.assign(Ncond-1,0.1);

        // toggles specifying the sites and amino-acids displaying fitness modulations across conditions
        // for each k=1..Ncond;
        // all toggles across all sites and amino-acids are iid Bernoulli of parameter shiftprob[k-1]
        toggle = new BidimIIDMultiBernoulli(Ncond-1,Nsite,Naa,shiftprob);
        if (!withtoggle)    {
            toggle->Reset();
        }

        // final amino-acid fitness profiles across sites
        // (*fitnessprofile)(k,i)[a]: fitness of amino-acid a for site i under condition k
        // deterministic functions of the gamma fitness parameters, the system of masks across sites and toggles across sites and conditions
        // for site i and amino-acid a:
        // if sitemaskarray[i][a] == 0 : fitnessprofile(k,i)[a] = maskepsilon across all conditions
        // if sitemaskarray[i][a] == 1 : fitnessprofile(k,i)[a] determined by the system of fitness(0:Ncond,i)[a] and toggle(1:Ncond,i)[a]
        // as in the simple DiffSelSparseModel
        fitnessprofile = new DiffSelDoublySparseFitnessArray(*fitness,*sitemaskarray,*toggle,Nlevel,maskepsilon);
        
        // codon matrices
        // per condition and per site
        condsubmatrixarray = new AADiffSelCodonMatrixBidimArray(*fitnessprofile,*GetCodonStateSpace(),*nucmatrix);

        // sub matrices per branch and per site
        submatrixarray = new SubMatrixSelector(*condsubmatrixarray,*branchalloc);
        // sub matrices for root, across sites
        rootsubmatrixarray = new RootSubMatrixSelector(*condsubmatrixarray);

        // create phyloprocess
        phyloprocess = new PhyloProcess(tree, codondata, branchlength, 0, submatrixarray, rootsubmatrixarray);
        phyloprocess->Unfold();

        // create suffstat arrays
        suffstatarray = new PathSuffStatBidimArray(Ncond,Nsite);
    }

    //! \brief set toggle status: 0: toggles all to 0, no move on them, 1:random toggles, MCMC move on them
    void SetWithToggles(int in)   {
        withtoggle = in;
    }

    //! \brief set estimation method for branch lengths
    //!
    //! Used in a multigene context.
    //! - mode == 2: global
    //! - mode == 1: gene specific, with hyperparameters estimated across genes
    //! - mode == 0: gene-specific, with fixed hyperparameters
    void SetBLMode(int in)   {
        blmode = in;
    }

    //! \brief set estimation method for nuc rates
    //!
    //! Used in a multigene context.
    //! - mode == 2: global
    //! - mode == 1: gene specific, with hyperparameters estimated across genes
    //! - mode == 0: gene-specific, with fixed hyperparameters
    void SetNucMode(int in) {
        nucmode = in;
    }

    //! \brief set estimation method for fitness hyperparameter (center of Dirichlet distribution)
    //!
    //! - mode == 3: fixed (uniform)
    //! - mode == 2: shared across genes, estimated 
    //! - mode == 1: gene specific, with hyperparameters estimated across genes
    //! - mode == 0: gene-specific, with fixed hyperparameters
    void SetFitnessCenterMode(int in)    {
        fitnesscentermode = in;
    }

    //! \brief set estimation method for fitness hyperparameter (shape of Dirichlet distribution)
    //!
    //! - mode == 3: fixed
    //! - mode == 2: shared across genes, estimated 
    //! - mode == 1: gene specific, with hyperparameters estimated across genes
    //! - mode == 0: gene-specific, with fixed hyperparameters
    void SetFitnessShapeMode(int in, double inshape = 0)    {
        fitnessshapemode = in;
        if (inshape > 0)    {
            fitnessshape = inshape;
        }
    }

    //! \brief set estimation method for site profile masks
    //!
    //! Used in a multigene context.
    //! - mode == 3: no mask
    //! - mode == 2: parameter (maskprob) shared across genes
    //! - mode == 1: gene-specific parameter (maskprob), hyperparameters estimated across genes
    //! - mode == 0: gene-specific parameter (maskprob) with fixed hyperparameters
    void SetMaskMode(int in)    {
        maskmode = in;
    }

    //! \brief set estimation method for background fitness (maskepsilon)
    void SetMaskEpsilonMode(int in)    {
        maskepsilonmode = in;
    }

    // ------------------
    // Update system
    // ------------------

    //! \brief set branch lengths to a new value
    //! 
    //! Used in a multigene context.
    void SetBranchLengths(const BranchSelector<double>& inbranchlength)    {
        branchlength->Copy(inbranchlength);
    }

    //! get a copy of branch lengths into array given as argument
    void GetBranchLengths(BranchArray<double>& inbranchlength) const    {
        inbranchlength.Copy(*branchlength);
    }

    //! set branch lengths hyperparameters to a new value (multi-gene analyses)
    void SetBranchLengthsHyperParameters(const BranchSelector<double>& inblmean, double inblinvshape)   {
        blhypermean->Copy(inblmean);
        blhyperinvshape = inblinvshape;
        branchlength->SetShape(1.0 / blhyperinvshape);
    }

    //! set nucleotide rates hyperparameters to a new value (multi-gene analyses)
    void SetNucRatesHyperParameters(const std::vector<double>& innucrelratehypercenter, double innucrelratehyperinvconc, const std::vector<double>& innucstathypercenter, double innucstathyperinvconc) {
        nucrelratehypercenter = innucrelratehypercenter;
        nucrelratehyperinvconc = innucrelratehyperinvconc;
        nucstathypercenter = innucstathypercenter;
        nucstathyperinvconc = innucstathyperinvconc;
    }

    //! set nucleotide rates to a new value (multi-gene analyses)
    void SetNucRates(const std::vector<double>& innucrelrate, const std::vector<double>& innucstat) {
        nucrelrate = innucrelrate;
        nucstat = innucstat;
        CorruptMatrices();
    }

    //! copy nucleotide rates into vectors given as arguments (multi-gene analyses)
    void GetNucRates(std::vector<double>& innucrelrate, std::vector<double>& innucstat) const {
        innucrelrate = nucrelrate;
        innucstat = nucstat;
    }

    //! \brief set value of background fitness of low-fitness amino-acids
    void SetMaskEpsilon(double in)  {
        maskepsilon = in;
        fitnessprofile->SetEpsilon(maskepsilon);
    }

    //! set shift prob hyperparameters (pi, shiftprobhypermean and hyperinvconc) to specified values (used in multi-gene context)
    void SetShiftProbHyperParameters(const vector<double>& inpi, const vector<double>& inshiftprobhypermean, const vector<double>& inshiftprobhyperinvconc)    {
        pi = inpi;
        shiftprobhypermean = inshiftprobhypermean;
        shiftprobhyperinvconc = inshiftprobhyperinvconc;
    }

    //! const access to shift prob vector
    const vector<double>& GetShiftProbVector() const {
        return shiftprob;
    }


    //! get a copy of fitness array (for all sites and amino-acids) for condition k
    void GetFitnessArray(int k, double* array) const  {
        int j = 0;
        for (int i=0; i<GetNsite(); i++)    {
            for (int a=0; a<Naa; a++)   {
                array[j++] = fitness->GetVal(k,i)[a];
            }
        }
    }

    //! get a copy of toggle values (for all sites and amino-acids) for condition k
    void GetShiftToggleArray(int k, int* array) const  {
        int j = 0;
        for (int i=0; i<GetNsite(); i++)    {
            for (int a=0; a<Naa; a++)   {
                array[j++] = toggle->GetVal(k-1,i)[a];
            }
        }
    }

    //! const ref access to toggles for condition k=1..Ncond
    const vector<vector<int> > & GetCondToggleArray(int k) const {
        return toggle->GetSubArray(k-1);
    }

    //! const ref access to masks across sites
    const vector<vector<int> >& GetMaskArray() const    {
        return sitemaskarray->GetArray();
    }

    double GetMaskEpsilon() const   {
        return maskepsilon;
    }

    void Update() override {
        if (blmode == 0)    {
            blhypermean->SetAllBranches(1.0/lambda);
        }
        UpdateMask();
		fitness->SetShape(fitnessshape);
        UpdateAll();
        ResampleSub(1.0);
    }

    //! update mask array
    void UpdateMask()   {
        sitemaskarray->SetPi(maskprob);
    }

    //! \brief dummy function that does not do anything.
    //! 
    //! Used for the templates of ScalingMove, SlidingMove and ProfileMove (defined in ProbModel),
    //! all of which require a void (*f)(void) function pointer to be called after changing the value of the focal parameter.
    void NoUpdate() {}

    //! \brief tell the nucleotide and the codon matrices that their parameters have changed and that it should be updated
    //!
    //! The matrices are not directly updated at that step. Instead, corruption is notified,
    //! such that the matrices know that they will have to recalculate whichever component is requested later on upon demand.
    void CorruptMatrices()  {
        CorruptNucMatrix();
        condsubmatrixarray->Corrupt();
    }

    //! \brief tell the nucleotide matrix that its parameters have changed and that it should be updated
    //!
    //! The matrix is not directly updated at that step. Instead, corruption is notified,
    //! such that the matrix knows that it will have to recalculate whichever component is requested later on upon demand.
    void CorruptNucMatrix() {
        nucmatrix->CopyStationary(nucstat);
        nucmatrix->CorruptMatrix();
    }

    //! update fitness profiles and matrices across all sites and conditions
    void UpdateAll() {
        fitnessprofile->SetEpsilon(maskepsilon);
        fitnessprofile->Update();
        CorruptMatrices();
    }

    //! update fitness profiles and matrices across all conditions for site i
    void UpdateSite(int i) {
        fitnessprofile->UpdateColumn(i);
        condsubmatrixarray->CorruptColumn(i);
    }

    // ---------------
    // log priors
    // ---------------

    //! \brief return total log prior
    //!
    //! Note: up to some multiplicative constant
    double GetLogPrior() const {
        double total = 0;
        if (blmode < 2) {
            total += BranchLengthsLogPrior();
        }
        if (nucmode < 2)    {
            total += NucRatesLogPrior();
        }
        if ((fitnessshapemode < 2) || (fitnesscentermode < 2))   {
            total += FitnessHyperLogPrior();
        }
        // not updated at all times 
        // total += FitnessLogPrior();
        if (maskmode < 2)   {
            total += MaskHyperLogPrior();
            total += MaskLogPrior();
        }
        total += ToggleHyperLogPrior();
        // not updated at all times
        // total += ToggleLogPrior();
        return total;
    }

    //! \brief log prior over hyperparameter of prior over branch lengths (here, lambda ~ exponential of rate 10)
	double BranchLengthsHyperLogPrior() const {
		return -log(10.0) - lambda / 10;
	}

    //! log prior over branch lengths (iid exponential of rate lambda)
	double BranchLengthsLogPrior() const {
		double ret = branchlength->GetLogProb();
        if (blmode == 0)    {
            ret += BranchLengthsHyperLogPrior();
        }
        return ret;
	}

    //! log prior over nuc rates rho and pi (uniform)
    double NucRatesLogPrior() const {
        double total = 0;
        total += Random::logDirichletDensity(nucrelrate,nucrelratehypercenter,1.0/nucrelratehyperinvconc);
        total += Random::logDirichletDensity(nucstat,nucstathypercenter,1.0/nucstathyperinvconc);
        return total;
    }

    //! log prior over fitness hyperparameters
    double FitnessHyperLogPrior() const {
        double ret = 0;
        if (fitnessshapemode < 2)   {
            ret += FitnessShapeLogPrior();
        }
        if (fitnesscentermode < 2)  {
            ret += FitnessCenterLogPrior();
        }
        return ret;
    }

    //! log prior over fitness center hyperparameter
    double FitnessCenterLogPrior() const {
        return 0;
    }

    //! log prior over fitness shape hyperparameter
    double FitnessShapeLogPrior() const {
        // exponential on shape
        return -fitnessshape;
    }

    //! log prior over input fitness parameters
    double FitnessLogPrior() const  {
        return fitness->GetLogProb();
    }

    //! log prior over mask array hyperparameter (maskprob: uniform between 0 and 1 -- could be hyperparameterized)
    double MaskHyperLogPrior() const  {
        if (maskepsilon > 1)    {
            return Random::INFPROB;
        }
        return 0;
    }

    //! log prior over mask arrays across sites (iid Bernoulli conditional on at least one entry being 1)
    double MaskLogPrior() const   {
        return sitemaskarray->GetLogProb();
    }

    //! log prior over mask array for site i
    double MaskLogPrior(int i) const   {
        return sitemaskarray->GetLogProb(i);
    }

    //! log prior over toggle array hyperparameters (shiftprob vector)
    double ToggleHyperLogPrior() const  {
        double total = 0;
        for (int k=1; k<Ncond; k++) {
            double alpha = shiftprobhypermean[k-1] / shiftprobhyperinvconc[k-1];
            double beta = (1-shiftprobhypermean[k-1]) / shiftprobhyperinvconc[k-1];
            if (shiftprob[k-1] != 0)    {
                total += log(pi[k-1]) + Random::logBetaDensity(shiftprob[k-1],alpha,beta);
            }
            else    {
                if (pi[k-1] == 1.0)   {
                    cerr << "error in ToggleHyperLogPrior: inf\n";
                    exit(1);
                }
                total += log(1 - pi[k-1]);
            }
        }
        return total;
    }

    //! log prior over toggle array (IID bernoulli)
    double ToggleLogPrior() const   {
        return toggle->GetLogProb();
    }

    //! return log likelihood
    double GetLogLikelihood() const { 
        return phyloprocess->GetLogLikelihood();
    }

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const {
        return GetLogPrior() + GetLogLikelihood();
    }

    // ---------------
    // collecting suff stats
    // ---------------

    //! \brief const access to array of length-pathsuffstats across branches
    const PoissonSuffStatBranchArray* GetLengthPathSuffStatArray() const {
        return lengthpathsuffstatarray;
    }

    //! collect generic sufficient statistics from substitution mappings
    void CollectPathSuffStat() {
        suffstatarray->Clear();
        suffstatarray->AddSuffStat(*phyloprocess,*branchalloc);
    }

    //! collect sufficient statistics for moving branch lengths (directly from the substitution mappings)
    void CollectLengthSuffStat()    {
		lengthpathsuffstatarray->Clear();
        lengthpathsuffstatarray->AddLengthPathSuffStat(*phyloprocess);
    }

    //! \brief return log prob of the current substitution mapping, as a function of the current codon substitution process
    double SuffStatLogProb() const   {
        return suffstatarray->GetLogProb(*condsubmatrixarray);
    }

    //! \brief return log prob of the current substitution mapping, as a function of the current codon substitution process, at site i
    double SiteSuffStatLogProb(int site) const   {
        return suffstatarray->GetLogProb(site,*condsubmatrixarray);
    }

    //! \brief return log prob of current branch lengths, as a function of branch lengths hyperparameter lambda
	double BranchLengthsHyperSuffStatLogProb()	const {
		return hyperlengthsuffstat.GetLogProb(1.0,lambda);
	}

    //! return log prob of current fitness parameters, conditional on their hyperparameters
	double FitnessHyperSuffStatLogProb()	const {
		return hyperfitnesssuffstat.GetLogProb(fitnessshape,fitnesscenter);
	}

    // ---------------
    // log probs for MH moves
    // ---------------

    //! \brief log prob factor to be recomputed when moving branch lengths hyperparameters (here, lambda)
    double BranchLengthsHyperLogProb() const {
        return BranchLengthsHyperLogPrior() + BranchLengthsHyperSuffStatLogProb();
    }

    //! \brief log prob factor to be recomputed when moving nucleotide mutation rate parameters (nucrelrate and nucstat)
    double NucRatesLogProb() const {
        return NucRatesLogPrior() + SuffStatLogProb();
    }

    //! \brief log prob factor to be recomputed when moving fitness hyperparameters
    double FitnessHyperLogProb() const  {
        return FitnessHyperLogPrior() + FitnessHyperSuffStatLogProb();
    }

    //! \brief log prob factor to be recomputed when moving mask hyperparameter pi
    double MaskLogProb() const  {
        return MaskHyperLogPrior() + MaskLogPrior();
    }

    //! \brief log prob factor to be recomputed when moving maskepsilon
    double MaskEpsilonLogProb() const  {
        return MaskHyperLogPrior() + SuffStatLogProb();
    }

    // ---------------
    // Moves
    // ---------------

    //! \brief complete MCMC move schedule
	double Move() override {
        gammanullcount = 0;
        ResampleSub(1.0);
        MoveParameters(3,20);
        return 1.0;
	}

    //! complete series of MCMC moves on all parameters (repeated nrep times)
    void MoveParameters(int nrep0, int nrep) {

        for (int rep0 = 0; rep0 < nrep0; rep0++) {

            if (blmode < 2)    {
                MoveBranchLengths();
            }

            CollectPathSuffStat();
            UpdateAll();

            for (int rep = 0; rep < nrep; rep++) {
                MoveBaselineFitness();
                CompMoveFitness();
                if (maskmode < 3)   {
                    MoveMasks(10);
                }
                if (maskmode < 2)   {
                    MoveMaskHyperParameters();
                }
                if (withtoggle) {
                    MoveFitnessShifts(10);
                    MoveShiftToggles(10);
                }
                if ((fitnessshapemode < 2) || (fitnesscentermode < 2))   {
                    MoveFitnessHyperParameters();
                }
                if (maskepsilonmode < 2)    {
                    MoveMaskEpsilon();
                }
            }

            if (nucmode < 2)    {
                MoveNucRates();
            }
        }

        UpdateAll();
    }

    //! Gibbs resample substitution mappings conditional on current parameter configuration
    void ResampleSub(double frac)   {
        CorruptMatrices();
		phyloprocess->Move(frac);
    }

    //! Gibbs resample branch lengths (based on sufficient statistics and current value of lambda)
	void ResampleBranchLengths()	{
        CollectLengthSuffStat();
		branchlength->GibbsResample(*lengthpathsuffstatarray);
	}


    //! MCMC move schedule on branch lengths 
    void MoveBranchLengths()    {
        ResampleBranchLengths();
        if (blmode == 0)    {
            MoveLambda();
        }
    }

    //! MH move on branch lengths hyperparameters (here, scaling move on lambda, based on suffstats for branch lengths)
	void MoveLambda()	{
		hyperlengthsuffstat.Clear();
		hyperlengthsuffstat.AddSuffStat(*branchlength);
        ScalingMove(lambda,1.0,10,&DiffSelDoublySparseModel::BranchLengthsHyperLogProb,&DiffSelDoublySparseModel::NoUpdate,this);
        ScalingMove(lambda,0.3,10,&DiffSelDoublySparseModel::BranchLengthsHyperLogProb,&DiffSelDoublySparseModel::NoUpdate,this);
        blhypermean->SetAllBranches(1.0/lambda);
	}

    //! MH moves on nucleotide rate parameters (nucrelrate and nucstat: using ProfileMove)
	void MoveNucRates()	{

        CorruptMatrices();

        ProfileMove(nucrelrate,0.1,1,10,&DiffSelDoublySparseModel::NucRatesLogProb,&DiffSelDoublySparseModel::CorruptMatrices,this);
        ProfileMove(nucrelrate,0.03,3,10,&DiffSelDoublySparseModel::NucRatesLogProb,&DiffSelDoublySparseModel::CorruptMatrices,this);
        ProfileMove(nucrelrate,0.01,3,10,&DiffSelDoublySparseModel::NucRatesLogProb,&DiffSelDoublySparseModel::CorruptMatrices,this);

        ProfileMove(nucstat,0.1,1,10,&DiffSelDoublySparseModel::NucRatesLogProb,&DiffSelDoublySparseModel::CorruptMatrices,this);
        ProfileMove(nucstat,0.01,1,10,&DiffSelDoublySparseModel::NucRatesLogProb,&DiffSelDoublySparseModel::CorruptMatrices,this);

        CorruptMatrices();
	}

    //! MH compensatory move on fitness parameters and hyper-parameters
    void CompMoveFitness()  {
        CompMoveFitness(1.0,10);
    }

    //! \brief MH compensatory move on fitness parameters and hyper-parameters
    //!
    //! for a given amino-acid.
    //! shift by a constant factor e all *active* gamma fitness parameters across all amino-acids and conditions
    double CompMoveFitness(double tuning, int nrep) {

        double nacc = 0;
        double ntot = 0;

        for (int rep = 0; rep < nrep; rep++) {
            for (int i = 0; i < Nsite; i++) {

                const vector<int>& mask = (*sitemaskarray)[i];

                double deltalogprob = 0;

                // calculate log prior for active fitness parameters across amino-acids and conditions before the move
                for (int k=0; k<Ncond; k++) {
                    for (int a=0; a<Naa; a++)   {
                        if ((mask[k]) && ((!k) || ((*toggle)(k-1,i)[a])))   {
                            double alpha = fitnessshape*fitnesscenter[a];
                            deltalogprob -= - Random::logGamma(alpha) + (alpha-1)*log((*fitness)(k,i)[a]) - (*fitness)(k,i)[a];
                        }
                    }
                }

                double m = tuning*(Random::Uniform() - 0.5);
                double e = exp(m);

                // multiply all active fitness parameters across amino-acids and conditions by e
                int n = 0;
                for (int k=0; k<Ncond; k++) {
                    for (int a=0; a<Naa; a++)   {
                        if ((mask[k]) && ((!k) || ((*toggle)(k-1,i)[a])))   {
                            (*fitness)(k,i)[a] *= e;
                            n++;
                        }
                    }
                }

                double loghastings = n * m;

                // calculate log prior for active fitness parameters across amino-acids and conditions after the move
                for (int k=0; k<Ncond; k++) {
                    for (int a=0; a<Naa; a++)   {
                        if ((mask[k]) && ((!k) || ((*toggle)(k-1,i)[a])))   {
                            double alpha = fitnessshape*fitnesscenter[a];
                            deltalogprob += - Random::logGamma(alpha) + (alpha-1)*log((*fitness)(k,i)[a]) - (*fitness)(k,i)[a];
                        }
                    }
                }

                deltalogprob += loghastings;

                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    // restore previous value for active fitness parameters
                    for (int k=0; k<Ncond; k++) {
                        for (int a=0; a<Naa; a++)   {
                            if ((mask[k]) && ((!k) || ((*toggle)(k-1,i)[a])))   {
                                (*fitness)(k,i)[a] /= e;
                            }
                        }
                    }
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    //! MH move schedule on baseline gamma fitness parameters (for condition k=0)
    void MoveBaselineFitness() {
        // if masks are not activated (all entries equal to 1), move a random subset of entries over the 20 amino-acids (2d parameter of call)
        if (maskmode == 3)  {
            MoveAllBaselineFitness(1.0, 3, 10);
            MoveAllBaselineFitness(1.0, 10, 10);
            MoveAllBaselineFitness(1.0, 20, 10);
            MoveAllBaselineFitness(0.3, 20, 10);
        }
        // if masks are activated, move all active entries
        else    {
            MoveBaselineFitness(1.0, 10);
            MoveBaselineFitness(0.3, 10);
        }
    }

    //! elementary MH move on baseline gamma fitness parameters (for condition k=0)
    double MoveAllBaselineFitness(double tuning, int n, int nrep) {

        double nacc = 0;
        double ntot = 0;
        vector<double> bk(Naa,0);


        for (int rep = 0; rep < nrep; rep++) {
            for (int i = 0; i < Nsite; i++) {

                vector<double>& x = (*fitness)(0,i);

                bk = x;

                double deltalogprob = -fitness->GetLogProb(0,i) - SiteSuffStatLogProb(i);
                double loghastings = Random::PosRealVectorProposeMove(x, Naa, tuning, n);
                deltalogprob += loghastings;

                UpdateSite(i);

                deltalogprob += fitness->GetLogProb(0,i) + SiteSuffStatLogProb(i);

                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    x = bk;
                    UpdateSite(i);
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    //! elementary MH move on baseline fitness parameters (for condition k=0): version used when masks are activated
    double MoveBaselineFitness(double tuning, int nrep) {

        double nacc = 0;
        double ntot = 0;
        vector<double> bk(Naa,0);


        for (int rep = 0; rep < nrep; rep++) {
            for (int i = 0; i < Nsite; i++) {

                vector<double>& fit = (*fitness)(0,i);
                const vector<int>& mask = (*sitemaskarray)[i];

                bk = fit;

                // calculate log prob before the move (for the prior factor, only those entries that are not masked are counted)
                double deltalogprob = -fitness->GetLogProb(0,i,mask) - SiteSuffStatLogProb(i);
                double loghastings = Random::PosRealVectorProposeMove(fit, Naa, tuning, mask);
                deltalogprob += loghastings;

                UpdateSite(i);

                deltalogprob += fitness->GetLogProb(0,i,mask) + SiteSuffStatLogProb(i);

                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    fit = bk;
                    UpdateSite(i);
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    //! MH move schedule on gamma fitness parameters (fitness shifts) for non-baseline conditions
    void MoveFitnessShifts(int nrep)    {
        for (int k=1; k<Ncond; k++) {
            MoveFitnessShifts(k,1,nrep);
            MoveFitnessShifts(k,0.3,nrep);
        }
    }

    //! elementary MH move on fitness shifts for non-baseline conditions
    double MoveFitnessShifts(int k, double tuning, int nrep) {

        double nacc = 0;
        double ntot = 0;
        vector<double> bk(Naa,0);

        for (int rep = 0; rep < nrep; rep++) {
            for (int i = 0; i < Nsite; i++) {

                vector<double>& x = (*fitness)(k,i);
                const vector<int>& t = (*toggle)(k-1,i);
                const vector<int>& m = sitemaskarray->GetVal(i);

                // compute condition-specific mask, which is the conjunction of baseline mask and condition-specific vector of toggles: s = m*t
                // this mask specifies which amino-acids are both active (across the tree) and undergoing a fitness shift in current condition
                vector<int> s(Naa,0);

                // nshift: number of amino-acids that are active in baseline and undergoing a fitness shift in current condition
                // nmask : number of amino-acids that are active in baseline
                int nshift = 0;
                int nmask = 0;
                for (int a=0; a<Naa; a++)   {
                    s[a] = t[a]*m[a];
                    nmask += m[a];
                    nshift += s[a];
                }

                // fitness shifts are impacting the final fitness only if 
                // (1) there are at least 2 active amino-acids
                // (2) at least one amino-acid is undergoing a fitness shift
                if ((nmask > 1) && nshift) {

                    bk = x;

                    // log prob before the move (for the prior factor, only those amino-acids that are concerned, such as specified by s, are counted)
                    double deltalogprob = -fitness->GetLogProb(k,i,s) - SiteSuffStatLogProb(i);

                    // propose move (only for the relevant amino-acids)
                    double loghastings = Random::PosRealVectorProposeMove(x, Naa, tuning, s);

                    deltalogprob += loghastings;

                    UpdateSite(i);

                    // log prob after the move
                    deltalogprob += fitness->GetLogProb(k,i,s) + SiteSuffStatLogProb(i);

                    int accepted = (log(Random::Uniform()) < deltalogprob);
                    if (accepted) {
                        nacc++;
                    } else {
                        x = bk;
                        UpdateSite(i);
                    }
                    ntot++;
                }
            }
        }
        return nacc / ntot;
    }

    //! MH moves on hyperparameters of distribution of fitness factors
    void MoveFitnessHyperParameters() {
        // collect suff stats across all active fitness parameters
        hyperfitnesssuffstat.Clear();
        hyperfitnesssuffstat.AddSuffStat(*fitness,*sitemaskarray,*toggle);

        if (fitnessshapemode < 2)   {
            ScalingMove(fitnessshape,1.0,100,&DiffSelDoublySparseModel::FitnessHyperLogProb,&DiffSelDoublySparseModel::NoUpdate,this);
            ScalingMove(fitnessshape,0.3,100,&DiffSelDoublySparseModel::FitnessHyperLogProb,&DiffSelDoublySparseModel::NoUpdate,this);
            ScalingMove(fitnessshape,0.1,100,&DiffSelDoublySparseModel::FitnessHyperLogProb,&DiffSelDoublySparseModel::NoUpdate,this);
        }
		fitness->SetShape(fitnessshape);

        if (fitnesscentermode < 2)  {
            ProfileMove(fitnesscenter,0.3,1,100,&DiffSelDoublySparseModel::FitnessHyperLogProb,&DiffSelDoublySparseModel::NoUpdate,this);
            ProfileMove(fitnesscenter,0.1,1,100,&DiffSelDoublySparseModel::FitnessHyperLogProb,&DiffSelDoublySparseModel::NoUpdate,this);
            ProfileMove(fitnesscenter,0.1,3,100,&DiffSelDoublySparseModel::FitnessHyperLogProb,&DiffSelDoublySparseModel::NoUpdate,this);
        }
    }

    //! Move schedule for Gibbs resampling of shifting probabilities
    void ResampleShiftProb()    {
        for (int k=1; k<Ncond; k++) {
            ResampleShiftProb(k);
        }
    }
    
    //! Gibbs resampling of shifting probability under condition k
    void ResampleShiftProb(int k)   {
        // pre-calculate parameters of the Beta distribution for non-zero case
        double alpha = shiftprobhypermean[k-1] / shiftprobhyperinvconc[k-1];
        double beta = (1-shiftprobhypermean[k-1]) / shiftprobhyperinvconc[k-1];

        // nshift: number of amino-acids that are active in baseline and undergoing a fitness shift in current condition
        // nmask : number of amino-acids that are active in baseline
        // both are summed across all sites: sufficient statistics for shiftprob
        int nshift = 0;
        int nmask = 0;
        for (int i=0; i<Nsite; i++) {
            const vector<int>& t = (*toggle)(k-1,i);
            const vector<int>& m = sitemaskarray->GetVal(i);
            int ns = 0;
            int nm = 0;
            for (int a=0; a<Naa; a++)   {
                ns += m[a]*t[a];
                nm += m[a];
            }

            // fitness shifts are counted (have an effect) only if there are at least 2 active amino-acids
            if (nm > 1) {
                nshift += ns;
                nmask += nm;
            }
        }

        if (nshift || (pi[k-1] == 1.0)) {
            shiftprob[k-1] = Random::BetaSample(alpha + nshift, beta + nmask - nshift);
        }
        else    {
            double logp0 = log(1-pi[k-1]);

            double logp1 = log(pi[k-1]);
            logp1 -= Random::logGamma(alpha) + Random::logGamma(beta) - Random::logGamma(alpha+beta);
            logp1 += Random::logGamma(alpha + nshift) + Random::logGamma(beta + nmask - nshift) - Random::logGamma(alpha+beta+nmask);

            double max = (logp0 > logp1) ? logp0 : logp1;
            double p0 = exp(logp0-max);
            double p1 = exp(logp1-max);
            double tot = p0+p1;
            p0/=tot;
            p1/=tot;

            if (Random::Uniform() < p0) {
                shiftprob[k-1] = 0;
            }
            else    {
                shiftprob[k-1] = Random::BetaSample(alpha + nshift, beta + nmask - nshift);
            }
        }
    }

    //! number of shifts in condition k
    double GetNShift(int k)  const {

        int nshift = 0;
        for (int i=0; i<Nsite; i++) {
            const vector<int>& t = (*toggle)(k-1,i);
            const vector<int>& m = sitemaskarray->GetVal(i);
            int ns = 0;
            int nm = 0;
            for (int a=0; a<Naa; a++)   {
                ns += m[a]*t[a];
                nm += m[a];
            }

            // fitness shifts are counted (have an effect) only if there are at least 2 active amino-acids
            if (nm > 1) {
                nshift += ns;
            }
        }
        return nshift;
    }

    //! number of amino acids allowed to undergo a shift
    double GetNTarget()  const {

        int nmask = 0;
        for (int i=0; i<Nsite; i++) {
            const vector<int>& m = sitemaskarray->GetVal(i);
            int nm = 0;
            for (int a=0; a<Naa; a++)   {
                nm += m[a];
            }

            // fitness shifts are counted (have an effect) only if there are at least 2 active amino-acids
            if (nm > 1) {
                nmask += nm;
            }
        }
        return nmask;
    }

    //! empirical fraction of allowed positions that undergo a shift
    double GetPropShift(int k)  const {

        // nshift: number of amino-acids that are active in baseline and undergoing a fitness shift in current condition
        // nmask : number of amino-acids that are active in baseline
        // both are summed across all sites: sufficient statistics for shiftprob
        int nshift = 0;
        int nmask = 0;
        for (int i=0; i<Nsite; i++) {
            const vector<int>& t = (*toggle)(k-1,i);
            const vector<int>& m = sitemaskarray->GetVal(i);
            int ns = 0;
            int nm = 0;
            for (int a=0; a<Naa; a++)   {
                ns += m[a]*t[a];
                nm += m[a];
            }

            // fitness shifts are counted (have an effect) only if there are at least 2 active amino-acids
            if (nm > 1) {
                nshift += ns;
                nmask += nm;
            }
        }

        return ((double) nshift) / nmask;
    }

    //! conditional posterior probability of the null hypothesis (null proportion of shifts) in condition k
    double GetProbNull(int k)  const {

        // pre-calculate parameters of the Beta distribution for non-zero case
        double alpha = shiftprobhypermean[k-1] / shiftprobhyperinvconc[k-1];
        double beta = (1-shiftprobhypermean[k-1]) / shiftprobhyperinvconc[k-1];

        // nshift: number of amino-acids that are active in baseline and undergoing a fitness shift in current condition
        // nmask : number of amino-acids that are active in baseline
        // both are summed across all sites: sufficient statistics for shiftprob
        int nshift = 0;
        int nmask = 0;
        for (int i=0; i<Nsite; i++) {
            const vector<int>& t = (*toggle)(k-1,i);
            const vector<int>& m = sitemaskarray->GetVal(i);
            int ns = 0;
            int nm = 0;
            for (int a=0; a<Naa; a++)   {
                ns += m[a]*t[a];
                nm += m[a];
            }

            // fitness shifts are counted (have an effect) only if there are at least 2 active amino-acids
            if (nm > 1) {
                nshift += ns;
                nmask += nm;
            }
        }

        if (nshift || (pi[k-1] == 1.0)) {
            return 0;
        }
        double logp0 = log(1-pi[k-1]);
        double logp1 = log(pi[k-1]);
        logp1 -= Random::logGamma(alpha) + Random::logGamma(beta) - Random::logGamma(alpha+beta);
        logp1 += Random::logGamma(alpha + nshift) + Random::logGamma(beta + nmask - nshift) - Random::logGamma(alpha+beta+nmask);

        double max = (logp0 > logp1) ? logp0 : logp1;
        double p0 = exp(logp0-max);
        double p1 = exp(logp1-max);
        double tot = p0+p1;
        p0/=tot;
        p1/=tot;
        return p0;
    }

    //! MH move schedule on mask hyperparameter (maskprob)
    void MoveMaskHyperParameters()  {
        SlidingMove(maskprob,1.0,10,0.05,0.975,&DiffSelDoublySparseModel::MaskLogProb,&DiffSelDoublySparseModel::UpdateMask,this);
        SlidingMove(maskprob,0.1,10,0.05,0.975,&DiffSelDoublySparseModel::MaskLogProb,&DiffSelDoublySparseModel::UpdateMask,this);
    }

    //! MH move schedule on background fitness (maskepsilon)
    void MoveMaskEpsilon()  {
        SlidingMove(maskepsilon,1.0,10,0,1.0,&DiffSelDoublySparseModel::MaskEpsilonLogProb,&DiffSelDoublySparseModel::UpdateAll,this);
        SlidingMove(maskepsilon,0.1,10,0,1.0,&DiffSelDoublySparseModel::MaskEpsilonLogProb,&DiffSelDoublySparseModel::UpdateAll,this);
        ScalingMove(maskepsilon,1.0,10,&DiffSelDoublySparseModel::MaskEpsilonLogProb,&DiffSelDoublySparseModel::UpdateAll,this);
        ScalingMove(maskepsilon,0.1,10,&DiffSelDoublySparseModel::MaskEpsilonLogProb,&DiffSelDoublySparseModel::UpdateAll,this);
    }

    //! MH move on masks
    double MoveMasks(int nrep)    {

		double nacc = 0;
		double ntot = 0;

        for (int i=0; i<Nsite; i++) {
            vector<int>& mask = (*sitemaskarray)[i];

            // compute number of active entries 
            int naa = 0;
            for (int a=0; a<Naa; a++)   {
                naa += mask[a];
            }

            for (int rep = 0; rep < nrep; rep++) {

                // randomly choose amino-acid
                int a = (int) (Naa * Random::Uniform());

                // don't propose move if this leads to 0 active entry in the end
                if ((!mask[a]) || (naa > 1))    {

                    // logprob before the move
                    double deltalogprob = -MaskLogPrior(i) - SiteSuffStatLogProb(i);

                    // do the move on the entry mask[k]
                    int oldnaa = naa;
                    naa -= mask[a];
                    mask[a] = 1-mask[a];
                    naa += mask[a];

                    // if move is from inactive to active
                    if (mask[a])    {
                        // resample baseline fitness
                        (*fitness)(0,i)[a] = Random::sGamma(fitnessshape * fitnesscenter[a]);
                        if (! (*fitness)(0,i)[a]) {
                            gammanullcount++;
                            (*fitness)(0,i)[a] = 1e-8;
                        }
                        // resample toggles and fitness shifts across all non-baseline conditions
                        for (int k=1; k<Ncond; k++) {
                            (*toggle)(k-1,i)[a] = (Random::Uniform() < pi[k-1]);
                            if ((*toggle)(k-1,i)[a])    {
                                (*fitness)(k,i)[a] = Random::sGamma(fitnessshape * fitnesscenter[a]);
                                if (! (*fitness)(k,i)[a]) {
                                    gammanullcount++;
                                    (*fitness)(k,i)[a] = 1e-8;
                                }
                            }  
                        }
                    }

                    // if move is from 1 to 2 active amino-acids,
                    // the other unmasked amino-acid should also redraw its toggles and fitness shifts (but not the baseline fitness)
                    if ((oldnaa == 1) && (naa == 2)) {
                        // grep the index of the other active amino-acid
                        int b = 0;
                        while ((b<Naa) && ((!mask[b]) || (b==a)))   {
                            b++;
                        }
                        if (b == Naa)   {
                            cerr << "error in MoveMasks, when choosing other amino-acid\n";
                            exit(1);
                        }
                        /*
                        int b = -1;
                        for (int c=0; c<Naa; c++)   {
                            if (mask[c] && (c!=a))  {
                                b = c;
                            }
                        }
                        if (b == -1)    {
                            cerr << "error in Move masks\n";
                            exit(1);
                        }
                        */

                        // resample toggles and fitness shifts across all non-baseline conditions
                        for (int k=1; k<Ncond; k++) {
                            (*toggle)(k-1,i)[b] = (Random::Uniform() < pi[k-1]);
                            if ((*toggle)(k-1,i)[b])    {
                                (*fitness)(k,i)[b] = Random::sGamma(fitnessshape * fitnesscenter[b]);
                                if (! (*fitness)(k,i)[b]) {
                                    gammanullcount++;
                                    (*fitness)(k,i)[b] = 1e-8;
                                }
                            }  
                        }
                    }

                    // logprob after the move
                    UpdateSite(i);
                    deltalogprob += MaskLogPrior(i) + SiteSuffStatLogProb(i);

                    // decision
                    int accepted = (log(Random::Uniform()) < deltalogprob);
                    if (accepted)	{
                        nacc ++;
                    }
                    else	{
                        // restore state (not necessary to restore all fitness parameters and toggles that have been modified during the move,
                        // since they all go back to inactive state anyway and will be resampled upon next activation)
                        naa -= mask[a];
                        mask[a] = 1-mask[a];
                        naa += mask[a];
                        UpdateSite(i);
                    }
                    ntot++;
                }
            }
        }
		return nacc/ntot;
	}

    //! MH move schedule on toggles
    void MoveShiftToggles(int nrep) {
        for (int k=1; k<Ncond; k++) {
            MoveShiftToggles(k,nrep);
        }
    }

    //! helper function: returns the marginal log prob of distribution of toggles for a condition,
    //! given number of toggles in active state and given hyperparameters
    double ToggleMarginalLogPrior(int nmask, int nshift, double pi, double alpha, double beta) const {
        double logp1 = log(pi) + Random::logGamma(alpha+beta) - Random::logGamma(alpha) - Random::logGamma(beta) + Random::logGamma(alpha+nshift) + Random::logGamma(beta + nmask - nshift) - Random::logGamma(alpha + beta + nmask);
        if (nshift || (pi == 1.0)) {
            return logp1;
        }
        double logp0 = log(1-pi);
        double max = (logp0 > logp1) ? logp0 : logp1;
        double tot = exp(logp0-max) + exp(logp1-max);
        double ret = log(tot) + max;
        return ret;
    }

    //! elementary MH move on toggles
    double MoveShiftToggles(int k, int nrep)  {

        // to achieve better MCMC mixing, shiftprob[k-1] is integrated out during this MH move on toggles (and Gibbs-resampled upon leaving this MH update)
        // nshift: number of amino-acids that are active in baseline and undergoing a fitness shift in current condition
        // nmask : number of amino-acids that are active in baseline
        // both are summed across all sites: sufficient statistics for shiftprob
        int nshift = 0;
        int nmask = 0;
        for (int i=0; i<Nsite; i++) {
            const vector<int>& t = (*toggle)(k-1,i);
            const vector<int>& m = sitemaskarray->GetVal(i);
            int ns = 0;
            int nm = 0;
            for (int a=0; a<Naa; a++)   {
                ns += m[a]*t[a];
                nm += m[a];
            }

            // fitness shifts are counted (have an effect) only if there are at least 2 active amino-acids
            if (nm > 1) {
                nmask += nm;
                nshift += ns;
            }
        }

        // pre-calculate parameters of the Beta distribution for non-zero case
        double alpha = shiftprobhypermean[k-1] / shiftprobhyperinvconc[k-1];
        double beta = (1-shiftprobhypermean[k-1]) / shiftprobhyperinvconc[k-1];
        double pp = pi[k-1];

        double ntot = 0;
        double nacc = 0;
        for (int rep = 0; rep < nrep; rep++) {
            for (int i = 0; i < Nsite; i++) {

                const vector<int>& t = (*toggle)(k-1,i);
                const vector<int>& m = sitemaskarray->GetVal(i);

                // compute nshift and nmask for this site only
                int ns = 0;
                int nm = 0;
                for (int a=0; a<Naa; a++)   {
                    ns += m[a]*t[a];
                    nm += m[a];
                }

                // do move only if there are at least 2 active amino-acids
                if (nm > 1) {

                    // randomly choose one active amino-acid (for which mask[a] == 1)
                    int b = (int) (nm * Random::Uniform()) + 1;
                    int a = 0;
                    while (b && (a<Naa))    {
                        if (m[a])   {
                            b--;
                        }
                        if (b) {
                            a++;
                        }
                    }
                    if (a == Naa)   {
                        cerr << "error in move shift toggles: overflow when choosing target amino-acid\n";
                        exit(1);
                    }
                    if (! m[a]) {
                        cerr << "error in move shift toggles: a is not within mask\n";
                        exit(1);
                    }

                    // 0 -> 1 case
                    if (!(*toggle)(k-1,i)[a])    {
                        double deltalogprob = -ToggleMarginalLogPrior(nmask,nshift,pp,alpha,beta) - SiteSuffStatLogProb(i);
                        (*toggle)(k-1,i)[a] = 1;
                        // redraw fitness parameter
                        (*fitness)(k,i)[a] = Random::sGamma(fitnessshape * fitnesscenter[a]);
                        if (! (*fitness)(k,i)[a]) {
                            gammanullcount++;
                            (*fitness)(k,i)[a] = 1e-8;
                        }
                        UpdateSite(i);
                        deltalogprob += ToggleMarginalLogPrior(nmask,nshift+1,pp,alpha,beta) + SiteSuffStatLogProb(i);

                        int accepted = (log(Random::Uniform()) < deltalogprob);
                        if (accepted) {
                            nacc++;
                            nshift++;
                        } else {
                            (*toggle)(k-1,i)[a] = 0;
                            UpdateSite(i);
                        }
                        ntot++;
                    }

                    // 1 -> 0 case
                    else    {
                        double deltalogprob = -ToggleMarginalLogPrior(nmask,nshift,pp,alpha,beta) - SiteSuffStatLogProb(i);
                        (*toggle)(k-1,i)[a] = 0;
                        UpdateSite(i);
                        deltalogprob += ToggleMarginalLogPrior(nmask,nshift-1,pp,alpha,beta) + SiteSuffStatLogProb(i);

                        int accepted = (log(Random::Uniform()) < deltalogprob);
                        if (accepted) {
                            nacc++;
                            nshift--;
                        } else {
                            (*toggle)(k-1,i)[a] = 1;
                            UpdateSite(i);
                        }
                        ntot++;
                    }
                }
            }
        }
        ResampleShiftProb(k);
        return nacc / ntot;
    }

    //-------------------
    // Accessors
    // ------------------

    //! const access to codon state space
	const CodonStateSpace* GetCodonStateSpace() const {
		return (CodonStateSpace*) codondata->GetStateSpace();
	}

    //! return number of aligned sites
    int GetNsite() const { return Nsite; }
    //! return number of conditions
    int GetNcond() const { return Ncond; }

    //-------------------
    // Traces and monitors
    // ------------------

    double GetMeanWidth() const {
        return sitemaskarray->GetMeanWidth();
    }

    void TraceHeader(ostream& os) const override {
        os << "#logprior\tlnL\tlength\t";
        os << "pi\t";
        os << "width\t";
        os << "epsilon\t";
        os << "shape\t";
        os << "center\t";
        for (int k = 1; k < Ncond; k++) {
            os << "prob" << k << '\t';
            os << "probnull" << k << '\t';
            os << "nshift" << k << '\t';
        }
        os << "statent\t";
        os << "rrent\t";
        os << "gammanulls\n";
    }

    void Trace(ostream& os) const override {
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood() << '\t';
        os << 3*branchlength->GetTotalLength() << '\t';
        os << maskprob << '\t';
        os << GetMeanWidth() << '\t';
        os << maskepsilon << '\t';
        os << fitnessshape << '\t';
        os << Random::GetEntropy(fitnesscenter) << '\t';
        for (int k=1; k<Ncond; k++) {
            os << shiftprob[k-1] << '\t';
            os << GetProbNull(k) << '\t';
            os << GetPropShift(k) << '\t';
        }
        os << Random::GetEntropy(nucstat) << '\t';
        os << Random::GetEntropy(nucrelrate) << '\t';
        os << gammanullcount << '\n';
    }

    //! trace the current value of toggles, across all sites and all amino-acids, under condition k (one single line in output stream)
    void TraceToggle(int k, ostream& os) const {
        for (int i=0; i<GetNsite(); i++)    {
            for (int a=0; a<Naa; a++)   {
                os << sitemaskarray->GetVal(i)[a] * toggle->GetVal(k-1,i)[a] << '\t';
            }
        }
        os << '\n';
    }

    //! trace the current value of fitness params, across all sites and all amino-acids, under condition k (one single line in output stream)
    void TraceFitness(int k, ostream& os) const {
        for (int i=0; i<GetNsite(); i++)    {
            for (int a=0; a<Naa; a++)   {
                os << fitness->GetVal(k,i)[a] << '\t';
            }
        }
        os << '\n';
    }

    void Monitor(ostream&) const override {}

    void FromStream(istream& is) override {
        if (blmode < 2) {
            is >> lambda;
            is >> *branchlength;
        }
        if (nucmode < 2)    {
            is >> nucrelrate;
            is >> nucstat;
        }
        if (fitnessshapemode < 2)   {
            is >> fitnessshape;
        }
        if (fitnesscentermode < 2)  {
            is >> fitnesscenter;
        }
        is >> *fitness;
        if (maskmode < 2)   {
            is >> maskprob;
        }
        if (maskmode < 3)   {
            is >> *sitemaskarray;
        }
        if (maskepsilonmode < 2)    {
            is >> maskepsilon;
        }
        is >> shiftprob;
        is >> *toggle;
    }

    void ToStream(ostream& os) const override {
        if (blmode < 2) {
            os << lambda << '\t';
            os << *branchlength << '\t';
        }
        if (nucmode < 2)    {
            os << nucrelrate << '\t';
            os << nucstat << '\t';
        }
        if (fitnessshapemode < 2)   {
            os << fitnessshape << '\t';
        }
        if (fitnesscentermode < 2)  {
            os << fitnesscenter << '\t';
        }
        os << *fitness << '\t';
        if (maskmode < 2)   {
            os << maskprob << '\t';
        }
        if (maskmode < 3)  {
            os << *sitemaskarray << '\t';
        }
        if (maskepsilonmode < 2)    {
            os << maskepsilon << '\t';
        }
        os << shiftprob << '\t';
        os << *toggle << '\t';
    }

    //! return size of model, when put into an MPI buffer (in multigene context -- only omegatree)
    unsigned int GetMPISize() const {
        int size = 0;
        if (blmode < 2) {
            size++;
            size += branchlength->GetMPISize();
        }
        if (nucmode < 2)    {
            size += nucrelrate.size();
            size += nucstat.size();
        }
        if (fitnessshapemode < 2)   {
            size ++;
        }
        if (fitnesscentermode < 2)  {
            size += fitnesscenter.size();
        }
        size += fitness->GetMPISize();
        if (maskmode < 2)   {
            size++;
        }
        if (maskmode < 3)   {
            size += sitemaskarray->GetMPISize();
        }
        if (maskepsilonmode < 2)    {
            size++;
        }
        size += shiftprob.size();
        size += toggle->GetMPISize();
        return size;
    }

    //! get array from MPI buffer
    void MPIGet(const MPIBuffer& is)    {
        if (blmode < 2) {
            is >> lambda;
            is >> *branchlength;
        }
        if (nucmode < 2)    {
            is >> nucrelrate;
            is >> nucstat;
        }
        if (fitnessshapemode < 2)   {
            is >> fitnessshape;
        }
        if (fitnesscentermode < 2)  {
            is >> fitnesscenter;
        }
        is >> *fitness;
        if (maskmode < 2)   {
            is >> maskprob;
        }
        if (maskmode < 3)   {
            is >> *sitemaskarray;
        }
        if (maskepsilonmode < 2)    {
            is >> maskepsilon;
        }
        is >> shiftprob;
        is >> *toggle;
    }

    //! write array into MPI buffer
    void MPIPut(MPIBuffer& os) const {
        if (blmode < 2) {
            os << lambda;
            os << *branchlength;
        }
        if (nucmode < 2)    {
            os << nucrelrate;
            os << nucstat;
        }
        if (fitnessshapemode < 2)   {
            os << fitnessshape;
        }
        if (fitnesscentermode < 2)  {
            os << fitnesscenter;
        }
        os << *fitness;
        if (maskmode < 2)   {
            os << maskprob;
        }
        if (maskmode < 3)   {
            os << *sitemaskarray;
        }
        if (maskepsilonmode < 2)    {
            os << maskepsilon;
        }
        os << shiftprob;
        os << *toggle;
    }
};
