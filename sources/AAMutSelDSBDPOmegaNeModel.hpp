#include "CodonSequenceAlignment.hpp"
#include "Tree.hpp"
#include "ProbModel.hpp"
#include "GTRSubMatrix.hpp"
#include "AAMutSelOmegaCodonSubMatrix.hpp"
#include "PhyloProcess.hpp"
#include "IIDGamma.hpp"
#include "GammaSuffStat.hpp"
#include "IIDDirichlet.hpp"
#include "CodonSuffStat.hpp"
#include "ProbModel.hpp"
#include "StickBreakingProcess.hpp"
#include "MultinomialAllocationVector.hpp"
#include "Chrono.hpp"
#include "Permutation.hpp"

/**
 * \brief The mutation-selection model with constant fitness landscape over the tree -- double Dirichlet process version, but with condition-dependent variation in Ne.
 *
 * The model is parameterized by
 * - a fixed unrooted phylogenetic tree tau, with branch lengths l = (l_j), for j running over branches
 * - a GTR nucleotide matrix Q = rho * pi, specifying the mutation process (assumed homogeneous across sites and lineages)
 * - an array of site-specific amino-acid fitness profiles F_ia, for site i and amino-acid a
 * - an omega multiplier, capturing deviations of the non-syn rate from the model (see Rodrigue and Lartillot, 2107); this parameter is fixed to 1 by default.
 * - a vector of effective population sizes Ne, one Ne per condition, capturing variations of the effective population sizes that correlate with the different conditions. One Ne (the first condition by convention) is fixed to 1.0. The number of conditions is fixed, as is their assignment to branches of the tree.

 *
 * Site-specific amino-acid fitness profiles are drawn from a Dirichlet process,
 * implemented using a truncated stick-breaking process, of concentration parameter kappa, and truncated at Ncat.
 *
 * The base distribution of this Dirichlet process, G_0, is itself a (truncated stick breaking) mixture of Dirichlet distributions,
 * parameterized by basekappa, and truncated at baseNcat. Each component of this mixture Dirichlet is parameterized by a center (a 20-frequency vector) and a concentration.
 *
 * Concerning the first-level stick-breaking process, by default, Ncat == 100 (can be changed using the  -ncat option).
 * As for baseNcat, it is equal to 1, in which case the base distribution G_0 is just a simple Dirichlet (as in Rodrigue et al, 2010). This simple model is probably the one that should be used by default for single-gene analyses. The more complex model with baseNcat > 1 is meant to be used in a multi-gene context (although, even in that case, mixing is still challenging, and not sure whether this model brings important improvement in the end). Thus, baseNcat = 1 appears to be the most reasonable model settings for now.
 *
 * Priors (in a single-gene context):
 * - branch lengths iid exponential, of rate lambda
 * - lambda exponential of rate 10
 * - rho and pi uniform Dirichlet
 * - omega: fixed to 1 or exponential of rate 1
 * - kappa: exponential of rate 0.1
 * - center of base distribution: uniform Dirichlet
 * - concentration of base distribution: exponential of mean 20
 * - Ne: fixed to 1 for the first condition, GammaIID for all other conditions.
 *
 * In a multi-gene context, shrinkage across genes can be applied to branch lengths, omega, nucleotide rate parameters (rho and pi), and to the parameters of the base distribution (center and concentration) -- see MultiGeneAAMutSelDSBDPModel.
 * Shrinkage across genes could also be applied to Ne, but will be more challenging to implement, because sufficient statistics need to be collected from the clients to the master for a move on Ne to be accepted. Current idea: server proposes M new Ne values, clients receive them all and compute M sets of sufficient statistics, and return all M of them. The master then combines them all, picks one (with the appropriate Hastings ratio), and returns the chosen one to the clients.
 *
 */

class AAMutSelDSBDPOmegaNeModel : public ProbModel {

	Tree* tree;
	FileSequenceAlignment* data;
	const TaxonSet* taxonset;
	CodonSequenceAlignment* codondata;

	int Nsite;
	int Ntaxa;
	int Nbranch;

	double lambda;
	BranchIIDGamma* branchlength;
	PoissonSuffStatBranchArray* lengthpathsuffstatarray;
	GammaSuffStat hyperlengthsuffstat;

	std::vector<double> nucstat;
	std::vector<double> nucrelrate;
	GTRSubMatrix* nucmatrix;

    // of mean omegahypermean and inverse shape parameter omegahyperinvshape
    double omegahypermean;
    double omegahyperinvshape;
	double omega;
	OmegaPathSuffStat omegapathsuffstat;

    // base distribution G0 is itself a stick-breaking mixture of Dirichlet distributions

    int baseNcat;
    double basekappa;
    StickBreakingProcess* baseweight;
    OccupancySuffStat* baseoccupancy;

    vector<double> basecenterhypercenter;
    double basecenterhyperinvconc;
    IIDDirichlet* basecenterarray;

    double baseconchypermean;
    double baseconchyperinvshape;
    IIDGamma* baseconcentrationarray;

    MultinomialAllocationVector* componentalloc;
    MixtureSelector<vector<double> >* componentcenterarray;
    MixtureSelector<double>* componentconcentrationarray;

    // aa fitness arrays across sites are a SBDP process of base G0 defined above
    int Ncat;
    double kappa;
    StickBreakingProcess* weight;
    OccupancySuffStat* occupancy;

    MultiDirichlet* componentaafitnessarray;
    DirichletSuffStatArray* basesuffstatarray;

	MultinomialAllocationVector* sitealloc;

    // an array of codon matrices (one for each distinct aa fitness profile)
	AAMutSelOmegaCodonSubMatrixArray* componentcodonmatrixarray;

	// number of diff Ne categories
	size_t Ncond;
	vector<double> Ne;

	// which branch is under which condition
	BranchAllocationSystem* branchalloc;



	// this one is used by PhyloProcess: has to be a Selector<SubMatrix>
	MixtureSelector<SubMatrix>* sitesubmatrixarray;

	PhyloProcess* phyloprocess;

	PathSuffStatArray* sitepathsuffstatarray;
	PathSuffStatArray* componentpathsuffstatarray;

    // 0: free wo shrinkage
    // 1: free with shrinkage
    // 2: shared across genes
    // 3: fixed

    // currently: shared across genes
    int blmode;
    // currently, free without shrinkage: shared across genes
    int nucmode;
    // currently, shared across genes.
    // free without shrinkage, only with baseNcat = 1
    // free with shrinkage: not really interesting
    int basemode;
    // currently: fixed or free with shrinkage
    int omegamode;

    Chrono aachrono;
    Chrono basechrono;
    Chrono totchrono;

    double acca1,acca2,acca3,acca4;
    double tota1,tota2,tota3,tota4;
    double accb1,accb2,accb3,accb4;
    double totb1,totb2,totb3,totb4;

	public:

    //-------------------
    // Construction and allocation
    // ------------------

    //! \brief constructor
    //!
    //! parameters:
    //! - datafile: name of file containing codon sequence alignment
    //! - treefile: name of file containing tree topology (and branch conditions, such as specified by branch names)
    //! - Ncat: truncation of the first-level stick-breaking process (by default: 100)
    //! - baseNcat: truncation of the second-level stick-breaking process (by default: 1)
    //! - Ncond: number of conditions (by default: 2)
	AAMutSelDSBDPOmegaNeModel(string datafile, string treefile, int inNcat, int inbaseNcat, int inNcond)   {

        blmode = 0;
        nucmode = 0;
        basemode = 0;
        omegamode = 3;

		data = new FileSequenceAlignment(datafile);
		codondata = new CodonSequenceAlignment(data, true);

		Nsite = codondata->GetNsite();    // # columns
		Ntaxa = codondata->GetNtaxa();

        Ncat = inNcat;
        if (Ncat == -1) {
            Ncat = Nsite;
            if (Ncat > 100)    {
                Ncat = 100;
            }
        }

        baseNcat = inbaseNcat;

		std::cerr << "-- Number of sites: " << Nsite << std::endl;
        cerr << "ncat : " << Ncat << '\n';
        cerr << "basencat : " << baseNcat << '\n';

		taxonset = codondata->GetTaxonSet();

		// get tree from file (newick format)
		tree = new Tree(treefile);

		// check whether tree and data fits together
		tree->RegisterWith(taxonset);

		tree->SetIndices();
		Nbranch = tree->GetNbranch();

		Ncond = inNcond;
		// specifies which condition for which branch
		branchalloc = new BranchAllocationSystem(*tree,Ncond);

    acca1=acca2=acca3=acca4=0;
    tota1=tota2=tota3=tota4=0;
    accb1=accb2=accb3=accb4=0;
    totb1=totb2=totb3=totb4=0;

		// Allocate();
	}

    //! \brief set estimation method for branch lengths
    //!
    //! - mode == 2: shared and estimated across genes
    //! - mode == 1: gene specific, with hyperparameters estimated across genes (with shrinkage)
    //! - mode == 0: gene-specific, with fixed hyperparameters (without shrinkage)
    //!
    //! for single-gene analyses, only mode 0 can be used.
    void SetBLMode(int mode)    {
        blmode = mode;
    }

    //! \brief set estimation method for nuc rates
    //!
    //! - mode == 2: shared and estimated across genes
    //! - mode == 1: gene specific, with hyperparameters estimated across genes (with shrinkage)
    //! - mode == 0: gene-specific, with fixed hyperparameters (without shrinkage)
    //!
    //! for single-gene analyses, only mode 0 can be used.
    void SetNucMode(int mode)   {
        nucmode = mode;
    }

    //! \brief set estimation method for nuc rates
    //!
    //! - mode == 3: fixed to 1
    //! - mode == 2: shared and estimated across genes: currently not implemented
    //! - mode == 1: gene specific, with hyperparameters estimated across genes (with shrinkage)
    //! - mode == 0: gene-specific, with fixed hyperparameters (without shrinkage)
    //!
    //! for single-gene analyses, either mode 3 and mode 0 can be used -- default mode is 3.
    void SetOmegaMode(int mode) {
        omegamode = mode;
    }

    //! \brief set estimation method for center and concentration parameters of base distribution
    //!
    //! - mode == 3: shared across genes and fixed to externally given empirical values: currently not implemented
    //! - mode == 2: shared and estimated across genes
    //! - mode == 1: gene specific, with hyperparameters estimated across genes (with shrinkage): currently not implemented
    //! - mode == 0: gene-specific, with fixed hyperparameters (without shrinkage)
    //!
    //! for single-gene analyses, either mode 3 and mode 0 can be used -- default mode is 3.
    void SetBaseMode(int mode)  {
        basemode = mode;
    }

    //! allocate the model (data structures)
	void Allocate()	{

        // branch lengths
		lambda = 10;
		branchlength = new BranchIIDGamma(*tree,1.0,lambda);
		lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);

        // nucleotide mutation matrix
		nucrelrate.assign(Nrr,0);
        Random::DirichletSample(nucrelrate,vector<double>(Nrr,1.0/Nrr),((double) Nrr));
		nucstat.assign(Nnuc,0);
        Random::DirichletSample(nucstat,vector<double>(Nnuc,1.0/Nnuc),((double) Nnuc));
		nucmatrix = new GTRSubMatrix(Nnuc,nucrelrate,nucstat,true);

        // base distribution (can be skipped)
        basekappa = 1.0;
        baseweight = new StickBreakingProcess(baseNcat,basekappa);
        baseoccupancy = new OccupancySuffStat(baseNcat);

        basecenterhypercenter.assign(Naa,1.0/Naa);
        basecenterhyperinvconc = 1.0/Naa;

        basecenterarray = new IIDDirichlet(baseNcat,basecenterhypercenter,1.0/basecenterhyperinvconc);
        basecenterarray->SetUniform();

        baseconchypermean = Naa;
        baseconchyperinvshape = 1.0;
        double alpha = 1.0 / baseconchyperinvshape;
        double beta = alpha / baseconchypermean;

        baseconcentrationarray = new IIDGamma(baseNcat,alpha,beta);
        for (int k=0; k<baseNcat; k++)  {
            (*baseconcentrationarray)[k] = 20.0;
        }
        // suff stats for component aa fitness arrays
        basesuffstatarray = new DirichletSuffStatArray(baseNcat,Naa);
        componentalloc = new MultinomialAllocationVector(Ncat,baseweight->GetArray());
        componentcenterarray = new MixtureSelector<vector<double> >(basecenterarray,componentalloc);
        componentconcentrationarray = new MixtureSelector<double>(baseconcentrationarray,componentalloc);


				//
        // Vector of Ne
        //

				Ne.push_back(1);
				// We set the shape and scales to 1.0
				for(size_t i = 1; i< Ncond; i++) {
					Ne.push_back(Random::GammaSample(1.0,1.0));
				}

        //
        // (truncated) Dirichlet mixture of aa fitness profiles
        //

        // Ncat fitness profiles iid from the base distribution
        componentaafitnessarray = new MultiDirichlet(componentcenterarray,componentconcentrationarray);


        // mixture weights (truncated stick breaking process)
        kappa = 1.0;
        weight = new StickBreakingProcess(Ncat,kappa);

        // site allocations to the mixture (multinomial allocation)
        sitealloc = new MultinomialAllocationVector(Nsite,weight->GetArray());

        // occupancy suff stats of site allocations (for resampling weights)
        occupancy = new OccupancySuffStat(Ncat);

        // global omega (fixed to 1 by default)
        omegahypermean = 1.0;
        omegahyperinvshape = 1.0;
		omega = 1.0;

        // Ncat mut sel codon matrices (based on the Ncat fitness profiles of the mixture)
        componentcodonmatrixarray = new AAMutSelOmegaCodonSubMatrixArray(GetCodonStateSpace(), nucmatrix, componentaafitnessarray, omega);

        // selector, specifying which codon matrix should be used for each site
        sitesubmatrixarray = new MixtureSelector<SubMatrix>(componentcodonmatrixarray,sitealloc);

		phyloprocess = new PhyloProcess(tree,codondata,branchlength,0,sitesubmatrixarray);
		phyloprocess->Unfold();

		sitepathsuffstatarray = new PathSuffStatArray(Nsite);
        componentpathsuffstatarray = new PathSuffStatArray(Ncat);
	}

    //-------------------
    // Accessors
    // ------------------

    //! const access to codon state space
	CodonStateSpace* GetCodonStateSpace() const {
		return (CodonStateSpace*) codondata->GetStateSpace();
	}

    //! return number of aligned sites
    int GetNsite() const    {
        return Nsite;
    }

    //! return current omega value
    double GetOmega() const {
        return omega;
    }

    //! \brief const access to array of length-pathsuffstats across branches
    const PoissonSuffStatBranchArray* GetLengthPathSuffStatArray() const {
        return lengthpathsuffstatarray;
    }

    //! \brief const access to array of suff stats for base mixture
    const DirichletSuffStatArray* GetBaseSuffStatArray() const   {
        return basesuffstatarray;
    }

    //! \brief const access to array of occupancy suff stats for base mixture
    const OccupancySuffStat* GetBaseOccupancies() const {
        return baseoccupancy;
    }

    //-------------------
    // Setting and updating
    // ------------------

    //! set branch lengths to a new value (multi-gene analyses)
    void SetBranchLengths(const BranchSelector<double>& inbranchlength)    {
        branchlength->Copy(inbranchlength);
    }

    //! set omega to new value (multi-gene analyses)
    void SetOmega(double inomega)   {
        omega = inomega;
        UpdateCodonMatrices();
    }

    //! set omega hyperparams to new value (multi-gene analyses)
    void SetOmegaHyperParameters(double inomegahypermean, double inomegahyperinvshape)   {
        omegahypermean = inomegahypermean;
        omegahyperinvshape = inomegahyperinvshape;
    }

    //! set nuc rates to new value (multi-gene analyses)
    void SetNucRates(const std::vector<double>& innucrelrate, const std::vector<double>& innucstat) {
        nucrelrate = innucrelrate;
        nucstat = innucstat;
        UpdateMatrices();
    }

    //! set base mixture concentration and center parameters to new value (multi-gene analyses)
    void SetBaseMixture(const Selector<vector<double> >& inbasecenterarray, const Selector<double>& inbaseconcentrationarray, const Selector<double>& inbaseweight, const Selector<int>& inpermut) {

        basecenterarray->Copy(inbasecenterarray);
        baseconcentrationarray->Copy(inbaseconcentrationarray);
        baseweight->Copy(inbaseweight);
        componentalloc->Permute(inpermut);
    }

    //! \brief tell the nucleotide matrix that its parameters have changed and that it should be updated
	void UpdateNucMatrix()	{
		nucmatrix->CopyStationary(nucstat);
		nucmatrix->CorruptMatrix();
	}

    //! \brief tell the codon matrices that their parameters have changed and that they should be updated
	void UpdateCodonMatrices()	{
        componentcodonmatrixarray->SetOmega(omega);
		componentcodonmatrixarray->UpdateCodonMatrices();
	}

    //! \brief tell codon matrix k that its parameters have changed and that it should be updated
    void UpdateCodonMatrix(int k)    {
        (*componentcodonmatrixarray)[k].CorruptMatrix();
    }

    //! \brief tell the nucleotide and the codon matrices that their parameters have changed and that it should be updated
    void UpdateMatrices()   {
        UpdateNucMatrix();
        UpdateCodonMatrices();
    }

    //! \brief dummy function that does not do anything.
    //!
    //! Used for the templates of ScalingMove, SlidingMove and ProfileMove (defined in ProbModel),
    //! all of which require a void (*f)(void) function pointer to be called after changing the value of the focal parameter.
    void NoUpdate() {}

    void Update() override {
        branchlength->SetScale(lambda);
        baseweight->SetKappa(basekappa);
        weight->SetKappa(kappa);
        UpdateBaseOccupancies();
        UpdateOccupancies();
        UpdateMatrices();
	    ResampleSub(1.0);
    }

    //-------------------
    // Priors and likelihood
    //-------------------

    //! \brief return total log prior (up to some constant)
    double GetLogPrior() const {
        double total = 0;
        if (blmode < 2) {
            total += BranchLengthsHyperLogPrior();
            total += BranchLengthsLogPrior();
        }
        if (nucmode < 2)    {
            total += NucRatesLogPrior();
        }
        if (basemode < 2)   {
            if (baseNcat > 1)   {
                total += BaseStickBreakingHyperLogPrior();
                total += BaseStickBreakingLogPrior();
            }
            total += BaseLogPrior();
        }
        total += StickBreakingHyperLogPrior();
        total += StickBreakingLogPrior();
        total += AALogPrior();
        if (omegamode < 2)  {
            total += OmegaLogPrior();
        }
        return total;
    }

    //! return log likelihood
	double GetLogLikelihood() const {
		return phyloprocess->GetLogLikelihood();
	}

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const   {
        return GetLogPrior() + GetLogLikelihood();
    }

    //! \brief log prior over hyperparameter of prior over branch lengths (here, lambda ~ exponential of rate 10)
	double BranchLengthsHyperLogPrior() const {
		return -lambda / 10;
	}

    //! log prior over branch lengths (iid exponential of rate lambda)
	double BranchLengthsLogPrior() const {
		return branchlength->GetLogProb();
	}

    //! log prior over omega (gamma of mean omegahypermean and inverse shape omegahyperinvshape)
	double OmegaLogPrior() const {
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
		return alpha * log(beta) - Random::logGamma(alpha) + (alpha-1) * log(omega) - beta*omega;
	}

    //! log prior over nuc rates rho and pi (uniform)
    double NucRatesLogPrior() const {
        return 0;
    }

    //! log prior over concentration parameters basekappa of mixture of base distribution
    double BaseStickBreakingHyperLogPrior() const   {
        return -basekappa/10;
    }

    //! log prior over weights of stick breaking process of base distribution
    double BaseStickBreakingLogPrior() const    {
        return baseweight->GetLogProb(basekappa);
    }

    //! log prior over concentration parameters kappa of stick-breaking mixture of amino-acid profiles
    double StickBreakingHyperLogPrior() const   {
        return -kappa/10;
    }

    //! log prior over weights of stick breaking process of amino-acid profiles
    double StickBreakingLogPrior() const    {
        return weight->GetLogProb(kappa);
    }

    //! log prior over base center and concentration parameters
    double BaseLogPrior() const {
        double total = 0;
        total += basecenterarray->GetLogProb();
        total += baseconcentrationarray->GetLogProb();
        if (std::isinf(total))   {
            cerr << "in BaseLogPrior: inf\n";
            exit(1);
        }
        return total;
    }

    //! log prior over base center and concentration parameters of component k of base distribution
    double BaseLogPrior(int k) const {
        double total = 0;
        total += basecenterarray->GetLogProb(k);
        total += baseconcentrationarray->GetLogProb(k);
        return total;
    }

    //! log prior of amino-acid fitness profiles
    double AALogPrior() const {
        return componentaafitnessarray->GetLogProb();
    }

    //! log prior of amino-acid fitness profile k
    double AALogPrior(int k) const {
        return componentaafitnessarray->GetLogProb(k);
    }

    //-------------------
    // Suff Stat and suffstatlogprobs
    //-------------------

    //! return log prob of the current substitution mapping, as a function of the current codon substitution process
	double PathSuffStatLogProb() const {
        return componentpathsuffstatarray->GetLogProb(*componentcodonmatrixarray);
	}

    //! return log prob of the substitution mappings over sites allocated to component k of the mixture
    double PathSuffStatLogProb(int k) const {
        return componentpathsuffstatarray->GetVal(k).GetLogProb(componentcodonmatrixarray->GetVal(k));
    }

    //! return log prob of current branch lengths, as a function of branch lengths hyperparameter lambda
	double BranchLengthsHyperSuffStatLogProb() const {
		return hyperlengthsuffstat.GetLogProb(1.0,lambda);
	}

    //! return log prob of first-level mixture components (i.e. all amino-acid profiles drawn from component k of the base distribution), as a function of the center and concentration parameters of this component
    double BaseSuffStatLogProb(int k) const   {
        return basesuffstatarray->GetVal(k).GetLogProb(basecenterarray->GetVal(k),baseconcentrationarray->GetVal(k));
    }

    //-------------------
    //  Log probs for MH moves
    //-------------------

    //! log prob factor to be recomputed when moving branch lengths hyperparameters (here, lambda)
    double BranchLengthsHyperLogProb() const {
        return BranchLengthsHyperLogPrior() + BranchLengthsHyperSuffStatLogProb();
    }

    //! log prob factor to be recomputed when moving nucleotide mutation rate parameters (nucrelrate and nucstat)
    double NucRatesLogProb() const {
        return NucRatesLogPrior() + PathSuffStatLogProb();
    }

    //! log prob factor to be recomputed when moving aa hyper params (center and concentration) for component k of base distribution
    double BaseLogProb(int k) const   {
        return BaseLogPrior(k) + BaseSuffStatLogProb(k);
    }

    //! log prob factor to be recomputed when moving basekappa, the concentration parameter of the second-level strick breaking process (base distribution)
    double BaseStickBreakingHyperLogProb() const {
        return BaseStickBreakingHyperLogPrior() + BaseStickBreakingLogPrior();
    }

    //! log prob factor to be recomputed when moving kappa, the concentration parameter of the first-level strick breaking process
    double StickBreakingHyperLogProb() const {
        return StickBreakingHyperLogPrior() + StickBreakingLogPrior();
    }

    //-------------------
    //  Collecting Suff Stats
    //-------------------

    //! collect sufficient statistics if substitution mappings across sites
	void CollectSitePathSuffStat()	{
		sitepathsuffstatarray->Clear();
        sitepathsuffstatarray->AddSuffStat(*phyloprocess);
	}

    //! gather site-specific sufficient statistics component-wise
	void CollectComponentPathSuffStat()	{
        componentpathsuffstatarray->Clear();
        componentpathsuffstatarray->Add(*sitepathsuffstatarray,*sitealloc);
    }

    //! collect sufficient statistics for moving branch lengths (directly from the substitution mappings)
    void CollectLengthSuffStat()    {
		lengthpathsuffstatarray->Clear();
		lengthpathsuffstatarray->AddLengthPathSuffStat(*phyloprocess);
    }

    //-------------------
    //  Moves
    //-------------------

    //! complete MCMC move schedule
	double Move()	{
        ResampleSub(1.0);
        MoveParameters(30);
        return 1.0;
    }

    //! Gibbs resample substitution mappings conditional on current parameter configuration
    void ResampleSub(double frac)   {
        UpdateMatrices();
		phyloprocess->Move(frac);
    }

    //! complete series of MCMC moves on all parameters (repeated nrep times)
    void MoveParameters(int nrep)   {
		for (int rep=0; rep<nrep; rep++)	{

            totchrono.Start();
            if (blmode < 2) {
                ResampleBranchLengths();
                MoveBranchLengthsHyperParameter();
            }

			CollectSitePathSuffStat();
            CollectComponentPathSuffStat();

            if (nucmode < 2)    {
                MoveNucRates();
            }

            if (omegamode < 2)  {
                MoveOmega();
            }

            aachrono.Start();
            MoveAAMixture(3);
            aachrono.Stop();

            basechrono.Start();
            if (basemode < 2)   {
                MoveBase(3);
            }
            basechrono.Stop();

            totchrono.Stop();
		}
	}

    //! MH move on base mixture
    void MoveBase(int nrep) {
        if (baseNcat > 1)   {
            ResampleBaseAlloc();
        }
        MoveBaseMixture(nrep);
    }

    //! Gibbs resample branch lengths (based on sufficient statistics and current value of lambda)
	void ResampleBranchLengths()	{
        CollectLengthSuffStat();
		branchlength->GibbsResample(*lengthpathsuffstatarray);
	}

    //! MH move on branch lengths hyperparameters (here, scaling move on lambda, based on suffstats for branch lengths)
	void MoveBranchLengthsHyperParameter()	{
		hyperlengthsuffstat.Clear();
		hyperlengthsuffstat.AddSuffStat(*branchlength);
        ScalingMove(lambda,1.0,10,&AAMutSelDSBDPOmegaNeModel::BranchLengthsHyperLogProb,&AAMutSelDSBDPOmegaNeModel::NoUpdate,this);
        ScalingMove(lambda,0.3,10,&AAMutSelDSBDPOmegaNeModel::BranchLengthsHyperLogProb,&AAMutSelDSBDPOmegaNeModel::NoUpdate,this);
		branchlength->SetScale(lambda);
	}

    //! MH move on omega
	void MoveOmega()	{

		omegapathsuffstat.Clear();
		omegapathsuffstat.AddSuffStat(*componentcodonmatrixarray,*componentpathsuffstatarray);
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
		omega = Random::GammaSample(alpha + omegapathsuffstat.GetCount(), beta + omegapathsuffstat.GetBeta());
		UpdateCodonMatrices();
	}

    //! MH move on nucleotide rate parameters
	void MoveNucRates()	{

        ProfileMove(nucrelrate,0.1,1,3,&AAMutSelDSBDPOmegaNeModel::NucRatesLogProb,&AAMutSelDSBDPOmegaNeModel::UpdateMatrices,this);
        ProfileMove(nucrelrate,0.03,3,3,&AAMutSelDSBDPOmegaNeModel::NucRatesLogProb,&AAMutSelDSBDPOmegaNeModel::UpdateMatrices,this);
        ProfileMove(nucrelrate,0.01,3,3,&AAMutSelDSBDPOmegaNeModel::NucRatesLogProb,&AAMutSelDSBDPOmegaNeModel::UpdateMatrices,this);

        ProfileMove(nucstat,0.1,1,3,&AAMutSelDSBDPOmegaNeModel::NucRatesLogProb,&AAMutSelDSBDPOmegaNeModel::UpdateMatrices,this);
        ProfileMove(nucstat,0.01,1,3,&AAMutSelDSBDPOmegaNeModel::NucRatesLogProb,&AAMutSelDSBDPOmegaNeModel::UpdateMatrices,this);
	}

    //! MCMC module for the mixture amino-acid fitness profiles
    void MoveAAMixture(int nrep)    {
        for (int rep=0; rep<nrep; rep++)  {
            MoveAAProfiles();
            ResampleEmptyComponents();
            ResampleAlloc();
            LabelSwitchingMove();
            ResampleWeights();
            MoveKappa();
            CollectComponentPathSuffStat();
            UpdateCodonMatrices();
        }
    }

    //! resample empty components of the mixture from prior
    void ResampleEmptyComponents()  {
        componentaafitnessarray->PriorResample(*occupancy);
        componentcodonmatrixarray->UpdateCodonMatrices(*occupancy);
    }

    //! MH move on amino-acid fitness profiles (occupied components only)
    void MoveAAProfiles()   {
        CompMoveAAProfiles(3);
        MulMoveAAProfiles(3);
    }

    //! MH move on amino-acid fitness profiles: additive compensated move on pairs of entries of the vector
    double CompMoveAAProfiles(int nrep) {
        accb1 += MoveAA(1.0,1,nrep);
        // accb2 += MoveAA(1.0,3,nrep);
        // accb3 += MoveAA(0.3,3,nrep);
        accb4 += MoveAA(0.1,3,nrep);
        totb1++;
        totb2++;
        totb3++;
        totb4++;
        return 1.0;
    }

    //! MH move on amino-acid fitness profiles: multiplicative move (using the Gamma representation of the Dirichlet)
    double MulMoveAAProfiles(int nrep) {
        acca1 += MoveAAGamma(3.0,nrep);
        acca2 += MoveAAGamma(1.0,nrep);
        // acca3 += MoveAAGamma(0.3,nrep);
        // acca4 += MoveAAGamma(0.1,nrep);
        tota1++;
        tota2++;
        tota3++;
        tota4++;
        return 1.0;
    }

    //! MH move on amino-acid fitness profiles: additive compensated move on pairs of entries of the vector
	double MoveAA(double tuning, int n, int nrep)	{
		double nacc = 0;
		double ntot = 0;
		double bk[Naa];
        for (int i=0; i<Ncat; i++) {
            if (occupancy->GetVal(i))   {
                vector<double>& aa = (*componentaafitnessarray)[i];
                for (int rep=0; rep<nrep; rep++)	{
                    for (int l=0; l<Naa; l++)	{
                        bk[l] = aa[l];
                    }
                    double deltalogprob = -AALogPrior(i) - PathSuffStatLogProb(i);
                    double loghastings = Random::ProfileProposeMove(aa,Naa,tuning,n);
                    deltalogprob += loghastings;
                    UpdateCodonMatrix(i);
                    deltalogprob += AALogPrior(i) + PathSuffStatLogProb(i);
                    int accepted = (log(Random::Uniform()) < deltalogprob);
                    if (accepted)	{
                        nacc ++;
                    }
                    else	{
                        for (int l=0; l<Naa; l++)	{
                            aa[l] = bk[l];
                        }
                        UpdateCodonMatrix(i);
                    }
                    ntot++;
                }
            }
        }
		return nacc/ntot;
	}

    //! helper function: log density of 20 gammas
    double GammaAALogPrior(const vector<double>& x, const vector<double>& aacenter, double aaconc) {
        double total = 0;
        for (int l=0; l<Naa; l++)   {
            total += (aaconc*aacenter[l] -1)*log(x[l]) - x[l] - Random::logGamma(aaconc*aacenter[l]);
        }
        return total;
    }

    //! MH move on amino-acid fitness profiles: multiplicative move (using the Gamma representation of the Dirichlet)
	double MoveAAGamma(double tuning, int nrep)	{

		double nacc = 0;
		double ntot = 0;
        for (int i=0; i<Ncat; i++) {
            if (occupancy->GetVal(i))   {

                double aaconc = componentconcentrationarray->GetVal(i);
                const vector<double>& aacenter = componentcenterarray->GetVal(i);

                vector<double>& aa = (*componentaafitnessarray)[i];
                vector<double> x(Naa,0);
                double z = Random::sGamma(aaconc);
                for (int l=0; l<Naa; l++)   {
                    x[l] = z*aa[l];
                }

                double bkz = z;
                vector<double> bkx = x;
                vector<double> bkaa = aa;

                for (int rep=0; rep<nrep; rep++)	{

                    double deltalogprob = -GammaAALogPrior(x,aacenter,aaconc) - PathSuffStatLogProb(i);

                    double loghastings = 0;
                    z = 0;
                    for (int l=0; l<Naa; l++)   {
                        double m = tuning * (Random::Uniform() - 0.5);
                        double e = exp(m);
                        x[l] *= e;
                        z += x[l];
                        loghastings += m;
                    }
                    for (int l=0; l<Naa; l++)   {
                        aa[l] = x[l]/z;
                        if (aa[l] < 1e-50)  {
                            aa[l] = 1e-50;
                        }
                    }

                    deltalogprob += loghastings;

                    UpdateCodonMatrix(i);

                    deltalogprob += GammaAALogPrior(x,aacenter,aaconc) + PathSuffStatLogProb(i);

                    int accepted = (log(Random::Uniform()) < deltalogprob);
                    if (accepted)	{
                        nacc ++;
                        bkaa = aa;
                        bkx = x;
                        bkz = z;
                    }
                    else	{
                        aa = bkaa;
                        x = bkx;
                        z = bkz;
                        UpdateCodonMatrix(i);
                    }
                    ntot++;
                }
            }
        }
        return nacc/ntot;
	}

    //! Gibbs resample mixture allocations
    void ResampleAlloc()    {
        vector<double> postprob(Ncat,0);
        for (int i=0; i<Nsite; i++) {
            GetAllocPostProb(i,postprob);
            sitealloc->GibbsResample(i,postprob);
        }
        UpdateOccupancies();
    }

    //! update mixture occupancy suff stats (for resampling mixture weights)
    void UpdateOccupancies()    {
        occupancy->Clear();
        occupancy->AddSuffStat(*sitealloc);
    }

    //! get allocation posterior probabilities for a given site
    void GetAllocPostProb(int site, vector<double>& postprob)    {

        double max = 0;
        const vector<double>& w = weight->GetArray();
        const PathSuffStat& suffstat = sitepathsuffstatarray->GetVal(site);
        for (int i=0; i<Ncat; i++) {
            double tmp = suffstat.GetLogProb(componentcodonmatrixarray->GetVal(i));
            postprob[i] = tmp;
            if ((!i) || (max < tmp))    {
                max = tmp;
            }
        }

        double total = 0;
        for (int i=0; i<Ncat; i++) {
            postprob[i] = w[i] * exp(postprob[i] - max);
            total += postprob[i];
        }

        for (int i=0; i<Ncat; i++) {
            postprob[i] /= total;
        }
    }

    //! MCMC sequence for label switching moves
    void LabelSwitchingMove()   {
        Permutation permut(Ncat);
        weight->LabelSwitchingMove(5,*occupancy,permut);
        sitealloc->Permute(permut);
        componentaafitnessarray->Permute(permut);
    }

    //! Gibbs resample mixture weights (based on occupancy suff stats)
    void ResampleWeights()  {
        weight->GibbsResample(*occupancy);
    }

    //! MH move on kappa, concentration parameter of the mixture
    void MoveKappa()    {
        ScalingMove(kappa,1.0,10,&AAMutSelDSBDPOmegaNeModel::StickBreakingHyperLogProb,&AAMutSelDSBDPOmegaNeModel::NoUpdate,this);
        ScalingMove(kappa,0.3,10,&AAMutSelDSBDPOmegaNeModel::StickBreakingHyperLogProb,&AAMutSelDSBDPOmegaNeModel::NoUpdate,this);
        weight->SetKappa(kappa);
    }

    //! MCMC module for the base mixture
    void MoveBaseMixture(int nrep)    {
        for (int rep=0; rep<nrep; rep++)  {
            MoveBaseComponents(10);
            ResampleBaseEmptyComponents();
            if (baseNcat > 1)   {
                BaseLabelSwitchingMove();
                ResampleBaseWeights();
                MoveBaseKappa();
            }
        }
    }

    //! MCMC module for moving the center and concentration parameters of the components of the the base mixture
    void MoveBaseComponents(int nrep) {

        CollectBaseSuffStat();
        for (int rep=0; rep<nrep; rep++)    {
            MoveBaseCenters(1.0,1);
            MoveBaseCenters(1.0,3);
            MoveBaseCenters(0.3,3);
            MoveBaseConcentrations(1.0);
            MoveBaseConcentrations(0.3);
        }
    }

    //! collect suff stats for moving center and concentration parameters of the base mixture
    void CollectBaseSuffStat()   {
        basesuffstatarray->Clear();
        componentaafitnessarray->AddSuffStat(*basesuffstatarray,*componentalloc);
    }

    //! MCMC module for moving the center parameters of the components of the the base mixture
    double MoveBaseCenters(double tuning, int n) {
		double nacc = 0;
		double ntot = 0;
        vector<double> bk(Naa,0);
        for (int k=0; k<baseNcat; k++)  {
            if (baseoccupancy->GetVal(k))  {
                vector<double>& aa = (*basecenterarray)[k];
                bk = aa;
                double deltalogprob = -BaseLogProb(k);
                double loghastings = Random::ProfileProposeMove(aa,Naa,tuning,n);
                deltalogprob += loghastings;
                deltalogprob += BaseLogProb(k);
                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted)	{
                    nacc ++;
                }
                else	{
                    aa = bk;
                }
                ntot++;
            }
        }
		return nacc/ntot;
	}

    //! MCMC module for moving the concentration parameters of the components of the the base mixture
    double MoveBaseConcentrations(double tuning)  {
		double nacc = 0;
		double ntot = 0;
        for (int k=0; k<baseNcat; k++)  {
            if (baseoccupancy->GetVal(k))  {
                double& c = (*baseconcentrationarray)[k];
                double bk = c;
                double deltalogprob = -BaseLogProb(k);
                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                c *= e;
                deltalogprob += m;
                deltalogprob += BaseLogProb(k);
                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted)	{
                    nacc ++;
                }
                else	{
                    c = bk;
                }
                ntot++;
            }
        }
		return nacc/ntot;
    }

    //! resample empty components of the base mixture from the prior
    void ResampleBaseEmptyComponents()  {
        basecenterarray->PriorResample(*baseoccupancy);
        baseconcentrationarray->PriorResample(*baseoccupancy);
    }

    //! Gibbs resample base mixture allocations
    void ResampleBaseAlloc()    {
        vector<double> postprob(baseNcat,0);
        for (int i=0; i<Ncat; i++) {
            GetBaseAllocPostProb(i,postprob);
            componentalloc->GibbsResample(i,postprob);
        }
        UpdateBaseOccupancies();
    }

    //! update base occupancy suff stats (for moving base weights)
    void UpdateBaseOccupancies()    {
        baseoccupancy->Clear();
        baseoccupancy->AddSuffStat(*componentalloc);
    }

    //! get allocation posterior probability of a component of the first-level mixture to the components of the second-level mixture
    void GetBaseAllocPostProb(int cat, vector<double>& postprob)    {

        double max = 0;
        const vector<double>& w = baseweight->GetArray();
        for (int i=0; i<baseNcat; i++) {
            double tmp = Random::logDirichletDensity(componentaafitnessarray->GetVal(cat),basecenterarray->GetVal(i),baseconcentrationarray->GetVal(i));
            postprob[i] = tmp;
            if ((!i) || (max < tmp))    {
                max = tmp;
            }
        }

        double total = 0;
        for (int i=0; i<baseNcat; i++) {
            postprob[i] = w[i] * exp(postprob[i] - max);
            total += postprob[i];
        }

        for (int i=0; i<baseNcat; i++) {
            postprob[i] /= total;
        }
    }

    //! MCMC sequence for label switching moves of the base mixture
    void BaseLabelSwitchingMove()   {
        Permutation permut(baseNcat);
        baseweight->LabelSwitchingMove(5,*baseoccupancy,permut);
        componentalloc->Permute(permut);
        basecenterarray->Permute(permut);
        baseconcentrationarray->Permute(permut);
        basesuffstatarray->Permute(permut);
    }

    //! Gibbs resample base mixture weights (based on occupancy suff stats)
    void ResampleBaseWeights()  {
        baseweight->GibbsResample(*baseoccupancy);
    }

    //! MH move on basekappa, concentration parameter of the base mixture
    void MoveBaseKappa()    {
        ScalingMove(basekappa,1.0,10,&AAMutSelDSBDPOmegaNeModel::BaseStickBreakingHyperLogProb,&AAMutSelDSBDPOmegaNeModel::NoUpdate,this);
        ScalingMove(basekappa,0.3,10,&AAMutSelDSBDPOmegaNeModel::BaseStickBreakingHyperLogProb,&AAMutSelDSBDPOmegaNeModel::NoUpdate,this);
        baseweight->SetKappa(basekappa);
    }

    //-------------------
    // Traces and Monitors
    // ------------------

    //! return number of occupied components in first-level mixture (mixture of amino-acid fitness profiles)
    int GetNcluster() const {

        int n = 0;
        for (int i=0; i<Ncat; i++)  {
            if (occupancy->GetVal(i))    {
                n++;
            }
        }
        return n;
    }

    //! return number of occupied components in base distribution
    int GetBaseNcluster() const {

        int n = 0;
        for (int i=0; i<baseNcat; i++)  {
            if (baseoccupancy->GetVal(i))    {
                n++;
            }
        }
        return n;
    }

    //! return mean entropy of amino-acd fitness profiles
    double GetMeanAAEntropy() const {
        return componentaafitnessarray->GetMeanEntropy();
    }

    //! return mean of concentration parameters of base distribution
    double GetMeanComponentAAConcentration() const {

        double tot = 0;
        for (int i=0; i<baseNcat; i++)  {
            tot += baseoccupancy->GetVal(i) * baseconcentrationarray->GetVal(i);
        }
        return tot / Ncat;
    }

    //! return mean entropy of centers of base distribution
    double GetMeanComponentAAEntropy() const {

        double tot = 0;
        for (int i=0; i<baseNcat; i++)  {
            tot += baseoccupancy->GetVal(i) * Random::GetEntropy(basecenterarray->GetVal(i));
        }
        return tot / Ncat;
    }

    //! return entropy of vector of nucleotide exchange rates
    double GetNucRREntropy() const  {
        return Random::GetEntropy(nucrelrate);
    }

    //! return entropy of vector of equilibrium nucleotide composition
    double GetNucStatEntropy() const    {
        return Random::GetEntropy(nucrelrate);
    }

	void TraceHeader(ostream& os) const override {
		os << "#logprior\tlnL\tlength\t";
		os << "omega\t";
        os << "ncluster\t";
        os << "kappa\t";
        if (baseNcat > 1)   {
            os << "basencluster\t";
            os << "basekappa\t";
        }
        os << "aaent\t";
        os << "meanaaconc\t";
        os << "aacenterent\t";
		os << "statent\t";
		os << "rrent\n";
	}

	void Trace(ostream& os) const override {
		os << GetLogPrior() << '\t';
		os << GetLogLikelihood() << '\t';
        // 3x: per coding site (and not per nucleotide site)
        os << 3*branchlength->GetTotalLength() << '\t';
		os << omega << '\t';
        os << GetNcluster() << '\t';
        os << kappa << '\t';
        if (baseNcat > 1)   {
            os << GetBaseNcluster() << '\t';
            os << basekappa << '\t';
        }
        os << GetMeanAAEntropy() << '\t';
		os << GetMeanComponentAAConcentration() << '\t';
        os << GetMeanComponentAAEntropy() << '\t';
		os << Random::GetEntropy(nucstat) << '\t';
		os << Random::GetEntropy(nucrelrate) << '\n';
	}

	void Monitor(ostream& os) const override {
        os << totchrono.GetTime() << '\t' << aachrono.GetTime() << '\t' << basechrono.GetTime() << '\n';
        os << "prop time in aa moves  : " << aachrono.GetTime() / totchrono.GetTime() << '\n';
        os << "prop time in base moves: " << basechrono.GetTime() / totchrono.GetTime() << '\n';
    }

	void FromStream(istream& is) override {

        if (blmode < 2) {
            is >> lambda;
            is >> *branchlength;
        }
        if (nucmode < 2)    {
            is >> nucrelrate;
            is >> nucstat;
        }
        if (basemode < 2)   {
            is >> basekappa;
            baseweight->FromStreamSB(is);
            // is >> *baseweight;
            is >> *componentalloc;
            is >> *basecenterarray;
            is >> *baseconcentrationarray;
        }
        is >> kappa;
        weight->FromStreamSB(is);
        // is >> *weight;
        is >> *componentaafitnessarray;
        is >> *sitealloc;
        if (omegamode < 2)  {
            is >> omega;
        }
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
        if (basemode < 2)   {
            os << basekappa << '\t';
            baseweight->ToStreamSB(os);
            // os << *baseweight << '\t';
            os << *componentalloc << '\t';
            os << *basecenterarray << '\t';
            os << *baseconcentrationarray << '\t';
        }
        os << kappa << '\t';
        weight->ToStreamSB(os);
        // os << *weight << '\t';
        os << *componentaafitnessarray << '\t';
        os << *sitealloc << '\t';
        if (omega < 2)  {
            os << omega << '\t';
        }
    }

    //! return size of model, when put into an MPI buffer (in multigene context -- only omegatree)
    unsigned int GetMPISize() const {
        int size = 0;
        if (blmode < 2) {
            size ++;
            size += branchlength->GetMPISize();
        }
        if (nucmode < 2)    {
            size += nucrelrate.size();
            size += nucstat.size();
        }
        if (basemode < 2)   {
            size ++;
            size += baseweight->GetMPISizeSB();
            size += componentalloc->GetMPISize();
            size += basecenterarray->GetMPISize();
            size += baseconcentrationarray->GetMPISize();
        }
        size ++;
        size += weight->GetMPISizeSB();
        size += componentaafitnessarray->GetMPISize();
        size += sitealloc->GetMPISize();
        if (omega < 2)  {
            size++;
        }
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
        if (basemode < 2)   {
            is >> basekappa;
            baseweight->MPIGetSB(is);
            is >> *componentalloc;
            is >> *basecenterarray;
            is >> *baseconcentrationarray;
        }
        is >> kappa;
        weight->MPIGetSB(is);
        is >> *componentaafitnessarray;
        is >> *sitealloc;
        if (omegamode < 2)  {
            is >> omega;
        }
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
        if (basemode < 2)   {
            os << basekappa;
            baseweight->MPIPutSB(os);
            os << *componentalloc;
            os << *basecenterarray;
            os << *baseconcentrationarray;
        }
        os << kappa;
        weight->MPIPutSB(os);
        os << *componentaafitnessarray;
        os << *sitealloc;
        if (omega < 2)  {
            os << omega;
        }
    }
};
