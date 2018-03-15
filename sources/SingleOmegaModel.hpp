#include "CodonSequenceAlignment.hpp"
#include "Tree.hpp"
#include "ProbModel.hpp"
#include "GTRSubMatrix.hpp"
#include "CodonSubMatrix.hpp"
#include "PhyloProcess.hpp"
#include "IIDGamma.hpp"
#include "CodonSuffStat.hpp"
#include "ProbModel.hpp"
#include "GammaSuffStat.hpp"

/**
 * \brief A standard site- and branch-homogeneous Muse and Gaut omega-codon model
 *
 * The model has the following structure:
 * - branch lengths iid Exponential of rate lambda
 * - nucleotide relative exchangeabilities and stationaries are uniform Dirichlet
 * - there is one single omega=dN/dS for all sites and across all branches
 * - prior over omega ~ Gamma(omegahypermean,omegahyperinvshape).
 * 
 * The 2 hyperparameters omegahypermean and hyperinvshape are fixed when this model is used in isolation.
 * On the other hand, this model can be used in a multigene context (see MultiGeneSingeOmegaModel),
 * in which case omegahypermean and hyperinvshape are estimated across genes.
 */


class SingleOmegaModel : public ProbModel {

    // tree and data
	Tree* tree;
	FileSequenceAlignment* data;
	const TaxonSet* taxonset;
	CodonSequenceAlignment* codondata;

	int Nsite;
	int Ntaxa;
	int Nbranch;

    // branch lengths iid expo (gamma of shape 1 and scale lambda)
    // where lambda is a hyperparameter
	double lambda;
	BranchIIDGamma* branchlength;

    // nucleotide exchange rates and equilibrium frequencies (stationary probabilities)
	std::vector<double> nucrelrate;
	std::vector<double> nucstat;

    // omega has a Gamma prior
    // of mean omegahypermean and inverse shape parameter omegahyperinvshape
    double omegahypermean;
    double omegahyperinvshape;
	double omega;

    // a nucleotide matrix (parameterized by nucrelrate and nucstat)
	GTRSubMatrix* nucmatrix;
    // a codon matrix (parameterized by nucmatrix and omega)
	MGOmegaCodonSubMatrix* codonmatrix;
	
	PhyloProcess* phyloprocess;

	// suffstats

    // suff stats for substitution paths
    // summed over all branches and over all sites
	PathSuffStat pathsuffstat;

    // path suff stat can be summarized in terms of 4x4 suff stats, as a function of nucleotide rates
    NucPathSuffStat nucpathsuffstat;

    // or, alternatively, collected as a simple Poisson suff stat, as a function of omega
	OmegaPathSuffStat omegapathsuffstat;

    // Poisson suffstats for substitution histories, as a function of branch lengths
	PoissonSuffStatBranchArray* lengthpathsuffstatarray;

    // suff stats branch lengths, as a function of their hyper parameter lambda
    // (bl are iid gamma, of scale parameter lambda)
	GammaSuffStat hyperlengthsuffstat;

	public:

    //-------------------
    // Construction and allocation
    // ------------------

    //! \brief constructor, parameterized by names of data and tree files
    //!
    //! Note: in itself, the constructor does not allocate the model;
    //! It only reads the data and tree file and register them together.
	SingleOmegaModel(string datafile, string treefile)  {

		data = new FileSequenceAlignment(datafile);
		codondata = new CodonSequenceAlignment(data, true);

		Nsite = codondata->GetNsite();    // # columns
		Ntaxa = codondata->GetNtaxa();

		taxonset = codondata->GetTaxonSet();

		// get tree from file (newick format)
		tree = new Tree(treefile);

		// check whether tree and data fits together
		tree->RegisterWith(taxonset);

		tree->SetIndices();
		Nbranch = tree->GetNbranch();
	}

    //! model allocation
	void Allocate()	{

		lambda = 10;
		branchlength = new BranchIIDGamma(*tree,1.0,lambda);
		lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);

		nucrelrate.assign(Nrr,0);
        Random::DirichletSample(nucrelrate,vector<double>(Nrr,1.0/Nrr),((double) Nrr));

		nucstat.assign(Nnuc,0);
        Random::DirichletSample(nucstat,vector<double>(Nnuc,1.0/Nnuc),((double) Nnuc));

		nucmatrix = new GTRSubMatrix(Nnuc,nucrelrate,nucstat,true);

        omegahypermean = 1.0;
        omegahyperinvshape = 1.0;
		omega = 1.0;
		codonmatrix = new MGOmegaCodonSubMatrix(GetCodonStateSpace(), nucmatrix, omega);

		phyloprocess = new PhyloProcess(tree,codondata,branchlength,0,codonmatrix);
        phyloprocess->Unfold();
	}

    //-------------------
    // Accessors
    // ------------------

    //! const access to codon state space
	CodonStateSpace* GetCodonStateSpace() const {
		return (CodonStateSpace*) codondata->GetStateSpace();
	}

    //! return current value of omega
    double GetOmega() const {
        return omega;
    }

    //-------------------
    // Setting and updating
    // ------------------

    //! \brief set omega to a new value 
    //! 
    //! Used in a multigene context.
    //! Notifies corruption to the codon matrix.
    void SetOmega(double inomega)   {
        omega = inomega;
        TouchCodonMatrix();
    }

    //! \brief set the hyperparameters of the gamma prior over omega
    //! 
    //! Used in a multigene context.
    void SetOmegaHyperParameters(double inomegahypermean, double inomegahyperinvshape)   {
        omegahypermean = inomegahypermean;
        omegahyperinvshape = inomegahyperinvshape;
    }

    //! \brief set branch lengths to a new value
    //! 
    //! Used in a multigene context.
    void SetBranchLengths(const BranchSelector<double>& inbranchlength)    {
        branchlength->Copy(inbranchlength);
    }

    //! \brief set nucleotide relative exchangeabilities (nucrelrate) and equilibrium frequencies (nucstat) to a new value
    //!
    //! Notifies corruption to the nucleotide and the codon matrices
    void SetNucRates(const std::vector<double>& innucrelrate, const std::vector<double>& innucstat) {
        nucrelrate = innucrelrate;
        nucstat = innucstat;
        TouchMatrices();
    }

    //! \brief tell the nucleotide matrix that its parameters have changed and that it should be updated
    //!
    //! The matrix is not directly updated at that step. Instead, corruption is notified,
    //! such that the matrix knows that it will have to recalculate whichever component is requested later on upon demand.
	void TouchNucMatrix()	{
		nucmatrix->CopyStationary(nucstat);
		nucmatrix->CorruptMatrix();
	}

    //! \brief tell the codon matrix that its parameters have changed and that it should be updated
    //!
    //! The matrix is not directly updated at that step. Instead, corruption is notified,
    //! such that the matrix knows that it will have to recalculate whichever component is requested later on upon demand.
	void TouchCodonMatrix()	{
		codonmatrix->SetOmega(omega);
		codonmatrix->CorruptMatrix();
	}
		
    //! \brief tell the nucleotide and the codon matrices that their parameters have changed and that they should be updated
    //!
    //! Just successive calls to TouchNucMatrix() and then TouchCodonMatrix();
    void TouchMatrices()   {
        TouchNucMatrix();
        TouchCodonMatrix();
    }

    //! \brief dummy function that does not do anything.
    //! 
    //! Used for the templates of ScalingMove, SlidingMove and ProfileMove (defined in ProbModel),
    //! all of which require a void (*f)(void) function pointer to be called after changing the value of the focal parameter.
    void NoUpdate() {}

    void Update() override {
        branchlength->SetScale(lambda);
	    TouchMatrices();
	    ResampleSub(1.0);
    }

    //-------------------
    // Posterior Predictive
    // ------------------

    void PostPred(string name)  {
        branchlength->SetScale(lambda);
	    TouchMatrices();
        phyloprocess->PostPredSample(name);
    }

    //-------------------
    // Priors and likelihood
    //-------------------

    //! \brief return total log prior
    //!
    //! Note: up to some multiplicative constant
    double GetLogPrior() const {
        double total = 0;
        total += BranchLengthsHyperLogPrior();
        total += BranchLengthsLogPrior();
        total += NucRatesLogPrior();
        total += OmegaLogPrior();
        return total;
    }

    //! return current value of likelihood (pruning-style, i.e. integrated over all substitution histories)
	double GetLogLikelihood() const {
		return phyloprocess->GetLogLikelihood();
	}

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const   {
        return GetLogPrior() + GetLogLikelihood();
    }

    //! \brief log prior over hyperparameters of prior over branch lengths (here, lambda ~ exponential of rate 10)
	double BranchLengthsHyperLogPrior()	const {
        // exponential of mean 10
		return -lambda / 10;
	}

    //! log prior over branch lengths (iid exponential of rate lambda)
	double BranchLengthsLogPrior()	const {
		return branchlength->GetLogProb();
	}

    //! log prior over nucleotide relative exchangeabilities (nucrelrate) and eq. freqs. (nucstat) -- uniform Dirichlet in both cases
    double NucRatesLogPrior() const {
        return 0;
    }

    //! log prior over omega (gamma of mean omegahypermean and shape 1/omegahyperinvshape)
	double OmegaLogPrior()	const {
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
        return Random::logGammaDensity(omega,alpha,beta);
	}

    //-------------------
    // Suff Stat and suffstatlogprobs
    //-------------------

    //! \brief const access to array of length-pathsuffstats across branches
    //!
    //! Useful for resampling branch lengths conditional on the current substitution mapping
    const PoissonSuffStatBranchArray* GetLengthPathSuffStatArray() const {
        return lengthpathsuffstatarray;
    }

    //! \brief const acess to nuc-pathsuffstat
    //! 
    //! Useful for resampling nucleotide relative exchangeabilities (nucrelrate) and equilibrium frequencies (nucstat)
    //! conditional on the current substitution mapping.
    const NucPathSuffStat& GetNucPathSuffStat() const {
        return nucpathsuffstat;
    }

    //! \brief return log prob of the current substitution mapping, as a function of the current codon substitution process
    //!
    //! Calculated using pathsuffstat (which summarizes all information about the substitution mapping)
    //! and the codonmatrix.
    //! Both pathsuffstat and codonmatrix are assumed to be updated.
	double PathSuffStatLogProb() const {
		return pathsuffstat.GetLogProb(*codonmatrix);
	}

    //! \brief return log prob of current substitution mapping, as a function of branch lengths
    //!
    //! Calculated using the lengthpathsuffstat
    //! (which summarizes all information about how the prob of the substitution mapping depends on branch lengths).
    //! lengthpathsuffstat is assumed to be updated.
	double BranchLengthsHyperSuffStatLogProb()	const {
		return hyperlengthsuffstat.GetLogProb(1.0,lambda);
	}

    //! \brief return log prob of current substitution mapping, as a function of nucleotide parameters (nucrelrate and nucstat)
    //!
    //! Calculated using nucpathsuffstat 
    //! (which summarizes all information about how the probability of the substitution mapping depends on nucleotide mutation rates)
    //! and the nucmatrix.
    //! Both nucpathsuffstat and nucmatrix are assumed to be updated.
    double NucRatesSuffStatLogProb() const {
        return nucpathsuffstat.GetLogProb(*nucmatrix,*GetCodonStateSpace());
    }

    //-------------------
    //  Log probs for MH moves
    //-------------------

    //! \brief log prob factor to be recomputed when moving branch lengths hyperparameters (here, lambda)
    //!
    //! simply: BranchLengthsHyperLogPrior() + BranchLengthsHyperSuffStatLogProb()
    double BranchLengthsHyperLogProb() const {
        return BranchLengthsHyperLogPrior() + BranchLengthsHyperSuffStatLogProb();
    }

    //! \brief log prob factor to be recomputed when moving nucleotide mutation rate parameters (nucrelrate and nucstat)
    //!
    //! simply: NucRatesLogPrior() + NucRatesSuffStatLogProb();
    double NucRatesLogProb() const {
        return NucRatesLogPrior() + NucRatesSuffStatLogProb();
    }

    //-------------------
    //  Moves 
    //-------------------

    //! \brief complete MCMC move schedule
	double Move()	{
        ResampleSub(1.0);
        MoveParameters(30);
        return 1.0;
	}

    //! Gibbs resample substitution mappings conditional on current parameter configuration
    void ResampleSub(double frac)   {
        TouchMatrices();
		phyloprocess->Move(frac);
    }

    //! complete series of MCMC moves on all parameters (repeated nrep times)
    void MoveParameters(int nrep)   {
		for (int rep=0; rep<nrep; rep++)	{

			ResampleBranchLengths();
			MoveBranchLengthsHyperParameter();

			CollectPathSuffStat();

			MoveOmega();
			MoveNucRates();
		}
    }

    //! collect sufficient statistics for moving branch lengths (directly from the substitution mappings)
    void CollectLengthSuffStat()    {
		lengthpathsuffstatarray->Clear();
        lengthpathsuffstatarray->AddLengthPathSuffStat(*phyloprocess);
    }

    //! Gibbs resample branch lengths (based on sufficient statistics and current value of lambda)
	void ResampleBranchLengths()	{
        CollectLengthSuffStat();
		branchlength->GibbsResample(*lengthpathsuffstatarray);
	}

    //! MH move on branch lengths hyperparameters (here, scaling move on lambda, based on suffstats for branch lengths)
	void MoveBranchLengthsHyperParameter()	{

		hyperlengthsuffstat.Clear();
		hyperlengthsuffstat.AddSuffStat(*branchlength);;
        ScalingMove(lambda,1.0,10,&SingleOmegaModel::BranchLengthsHyperLogProb,&SingleOmegaModel::NoUpdate,this);
        ScalingMove(lambda,0.3,10,&SingleOmegaModel::BranchLengthsHyperLogProb,&SingleOmegaModel::NoUpdate,this);
		branchlength->SetScale(lambda);
	}

    //! collect generic sufficient statistics from substitution mappings
	void CollectPathSuffStat()	{
		pathsuffstat.Clear();
        pathsuffstat.AddSuffStat(*phyloprocess);
	}

    //! Gibbs resample omega (based on sufficient statistics of current substitution mapping)
	void MoveOmega()	{

		omegapathsuffstat.Clear();
		omegapathsuffstat.AddSuffStat(*codonmatrix,pathsuffstat);
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
		omega = Random::GammaSample(alpha + omegapathsuffstat.GetCount(), beta + omegapathsuffstat.GetBeta());
		TouchCodonMatrix();
	}

    //! collect sufficient statistics for moving nucleotide rates (based on generic sufficient statistics stored in pathsuffstat)
    void CollectNucPathSuffStat()    {
        TouchMatrices();
        nucpathsuffstat.Clear();
        nucpathsuffstat.AddSuffStat(*codonmatrix,pathsuffstat);
    }

    //! MH moves on nucleotide rate parameters (nucrelrate and nucstat: using ProfileMove)
	void MoveNucRates()	{

        CollectNucPathSuffStat();

        ProfileMove(nucrelrate,0.1,1,3,&SingleOmegaModel::NucRatesLogProb,&SingleOmegaModel::TouchNucMatrix,this);
        ProfileMove(nucrelrate,0.03,3,3,&SingleOmegaModel::NucRatesLogProb,&SingleOmegaModel::TouchNucMatrix,this);
        ProfileMove(nucrelrate,0.01,3,3,&SingleOmegaModel::NucRatesLogProb,&SingleOmegaModel::TouchNucMatrix,this);

        ProfileMove(nucstat,0.1,1,3,&SingleOmegaModel::NucRatesLogProb,&SingleOmegaModel::TouchNucMatrix,this);
        ProfileMove(nucstat,0.01,1,3,&SingleOmegaModel::NucRatesLogProb,&SingleOmegaModel::TouchNucMatrix,this);

        TouchMatrices();
	}

    //-------------------
    // Traces and Monitors
    // ------------------

	void TraceHeader(ostream& os) const override {
		os << "#logprior\tlnL\tlength\tlambda\t";
		os << "omega\t";
		os << "pi_a\tpi_c\tpi_g\tpi_t\t";
        os << "r_ac\tr_ag\tr_at\tr_cg\tr_ct\tr_gt\n";
	}

	void Trace(ostream& os) const override {	
		os << GetLogPrior() << '\t';
		os << GetLogLikelihood() << '\t';
        os << branchlength->GetTotalLength() << '\t';
		os << lambda << '\t';
		os << omega << '\t';
		os << nucstat << nucrelrate << '\n';
	}

	void Monitor(ostream& os) const {}

	void ToStream(ostream& os) const {
        os << omega << '\t';
        os << nucstat << '\t';
        os << nucrelrate << '\t';
        os << lambda << '\t';
        os << *branchlength << '\n';
    }

	void FromStream(istream& is) {
        is >> omega;
        is >> nucstat;
        is >> nucrelrate;
        is >> lambda;
        is >> *branchlength;
    }
    
};


