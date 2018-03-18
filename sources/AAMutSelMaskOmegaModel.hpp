
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
 * \brief The mutation-selection model with constant fitness landscape over the tree -- masked IID Dirichlet version
 *
 */


class IIDProfileMask : public SimpleArray<vector<int> >     {

    public:

    IIDProfileMask(int size, int indim, double pi) : SimpleArray(size,vector<int>(indim,1)), dim(indim), pi(0.1) {}

    int GetDim() const  {
        return dim;
    }

    void SetPi(double inpi) {
        pi = inpi;
    }

    double GetLogProb() const   {
        double total = 0;
        for (int i=0; i<GetSize(); i++) {
            total += GetLogProb(i);
        }
        return total;
    }

    double GetLogProb(int i) const  {
        int naa = 0;
        const vector<int>& x = GetVal(i);
        for (int k=0; k<GetDim(); k++)  {
            naa += x[k];
        }
        if (! naa)  {
            cerr << "error in IIDProfileMask: all entries are null\n";
            exit(1);
        }
        // probability is conditional on at least one entry being 1
        return naa*log(pi) + (GetDim()-naa)*log(1.0-pi) - log(1.0 - exp(GetDim()*log(1.0-pi)));
    }

    double GetMeanWidth() const {
        double mean = 0;
        for (int i=0; i<GetSize(); i++) {
            const vector<int>& x = GetVal(i);
            for (int k=0; k<GetDim(); k++)  {
                mean += x[k];
            }
        }
        mean /= GetSize();
        return mean;
    }

    private:
    int dim;
    double pi;
};

class MultiMaskWeight : public SimpleArray<vector<double> > {

    public: 
    
    MultiMaskWeight(const Selector<vector<int> >& inmaskarray, double inepsilon, const vector<double>& incenter, double inconcentration) : SimpleArray<vector<double> >(inmaskarray.GetSize(),vector<double>(inmaskarray.GetVal(0).size(),0)), dim(inmaskarray.GetVal(0).size()), maskarray(inmaskarray), epsilon(inepsilon), center(incenter), concentration(inconcentration) {
        Update();
    }

    int GetDim() const  {
        return dim;
    }

    void SetEpsilon(double inepsilon)   {
        epsilon = inepsilon;
    }

    void SetConcentration(double inconcentration)   {
        concentration = inconcentration;
    }

    void Update()   {
        for (int i=0; i<GetSize(); i++) {
            Update(i);
        }
    }

    void Update(int i)  {
        const vector<int>& mask = maskarray.GetVal(i);
        vector<double>& x = (*this)[i];
        for (int k=0; k<GetDim(); k++)  {
            if (mask[k] == 1)   {
                x[k] = concentration*center[k]*(1-epsilon);
            }
            else    {
                // x[k] = epsilon;
                x[k] = concentration*center[k]*epsilon;
            }
            if (x[k] < 0)   {
                cerr << "error: negative weight\n";
                cerr << mask[k] << '\t' << concentration << '\t' << center[k] << '\t' << epsilon << '\n';
                exit(1);
            }
        }
    }

    private:

    int dim;
    const Selector<vector<int> >& maskarray;
    double epsilon;
    const vector<double>& center;
    double concentration;
};

class AAMutSelMaskOmegaModel : public ProbModel {

	Tree* tree;
	FileSequenceAlignment* data;
	const TaxonSet* taxonset;
	CodonSequenceAlignment* codondata;

	int Nsite;
	int Ntaxa;
	int Nbranch;

	double lambda;
	BranchIIDGamma* blhypermean;
    double blhyperinvshape;
    GammaWhiteNoise* branchlength;
	PoissonSuffStatBranchArray* lengthpathsuffstatarray;
	GammaSuffStat hyperlengthsuffstat;

    // nucleotide rates hyperparameters
    vector<double> nucrelratehypercenter;
    double nucrelratehyperinvconc;
    vector<double> nucstathypercenter;
    double nucstathyperinvconc;

	std::vector<double> nucstat;
	std::vector<double> nucrelrate;
	GTRSubMatrix* nucmatrix;

    // of mean omegahypermean and inverse shape parameter omegahyperinvshape
    double omegahypermean;
    double omegahyperinvshape;
	double omega;
	OmegaPathSuffStat omegapathsuffstat;
	
    double pi;
    // 0.01
    double epsilon0;
    // uniform (0,1)
    double epsilon;
    vector<double> aacenter;
    double aaconcentration;

    IIDProfileMask* sitemaskarray;
    MultiMaskWeight* siteweightarray;
    MultiDirichlet* siteaafitnessarray;

    // an array of site-specific codon matrices
	AAMutSelOmegaCodonSubMatrixArray* sitecodonmatrixarray;

	PhyloProcess* phyloprocess;

	PathSuffStatArray* sitepathsuffstatarray;

    // 0: free wo shrinkage
    // 1: free with shrinkage
    // 2: shared across genes
    // 3: fixed

    // currently: shared across genes
    int blmode;
    // currently, free without shrinkage: shared across genes
    int nucmode;
    // currently: fixed or free with shrinkage
    int omegamode;

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
	AAMutSelMaskOmegaModel(string datafile, string treefile, int inomegamode)   {

        blmode = 0;
        nucmode = 0;
        omegamode = inomegamode;

		data = new FileSequenceAlignment(datafile);
		codondata = new CodonSequenceAlignment(data, true);

		Nsite = codondata->GetNsite();    // # columns
		Ntaxa = codondata->GetNtaxa();

		std::cerr << "-- Number of sites: " << Nsite << std::endl;

		taxonset = codondata->GetTaxonSet();

		// get tree from file (newick format)
		tree = new Tree(treefile);

		// check whether tree and data fits together
		tree->RegisterWith(taxonset);

		tree->SetIndices();
		Nbranch = tree->GetNbranch();
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

    //! allocate the model (data structures)
	void Allocate()	{

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

        // base distribution for masks: parameterized by pi
        pi = 0.1;
        // total prior weight of minor amino-acids
        epsilon0 = 0.1;
        epsilon = 0.05;
        // base center and concentration for the 20 amino-acids
        aacenter.assign(Naa,1.0/Naa);
        aaconcentration = Naa;

        // vector of Ncat masks iid from G0
        sitemaskarray = new IIDProfileMask(Nsite,Naa,pi);
        // 20-dimensional weight vectors (for Dirichlet distributions)
        siteweightarray = new MultiMaskWeight(*sitemaskarray,epsilon,aacenter,aaconcentration);
        // site-specific fitness profiles Dirichlet from the corresponding weights
        siteaafitnessarray = new MultiDirichlet(siteweightarray);

        // global omega (fixed to 1 by default)
        omegahypermean = 1.0;
        omegahyperinvshape = 1.0;
		omega = 1.0;

        // mut sel codon matrices (based on the fitness profiles of the mixture)
        sitecodonmatrixarray = new AAMutSelOmegaCodonSubMatrixArray(GetCodonStateSpace(), nucmatrix, siteaafitnessarray, omega);

		phyloprocess = new PhyloProcess(tree,codondata,branchlength,0,sitecodonmatrixarray);
		phyloprocess->Unfold();

		sitepathsuffstatarray = new PathSuffStatArray(Nsite);
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

    //-------------------
    // Setting and updating
    // ------------------

    //! set branch lengths to a new value (multi-gene analyses)
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
        UpdateMatrices();
    }

    //! copy nucleotide rates into vectors given as arguments (multi-gene analyses)
    void GetNucRates(std::vector<double>& innucrelrate, std::vector<double>& innucstat) const {
        innucrelrate = nucrelrate;
        innucstat = nucstat;
    }

    //! \brief tell the nucleotide matrix that its parameters have changed and that it should be updated
	void UpdateNucMatrix()	{
		nucmatrix->CopyStationary(nucstat);
		nucmatrix->CorruptMatrix();
	}

    //! \brief tell the codon matrices that their parameters have changed and that they should be updated
	void UpdateCodonMatrices()	{
        sitecodonmatrixarray->SetOmega(omega);
		sitecodonmatrixarray->UpdateCodonMatrices();
	}

    //! \brief tell codon matrix i that its parameters have changed and that it should be updated
    void UpdateCodonMatrix(int i)    {
        (*sitecodonmatrixarray)[i].CorruptMatrix();
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
        if (blmode == 0)    {
            blhypermean->SetAllBranches(1.0/lambda);
        }
        UpdateMask();
        UpdateAA();
        UpdateMatrices();
	    ResampleSub(1.0);
    }

    void UpdateMask()   {
        sitemaskarray->SetPi(pi);
    }

    void UpdateAA() {
        siteweightarray->SetEpsilon(epsilon);
        siteweightarray->SetConcentration(aaconcentration);
        siteweightarray->Update();
    }

    void UpdateAA(int i)    {
        siteweightarray->Update(i);
    }

    //-------------------
    // Priors and likelihood
    //-------------------

    //! \brief return total log prior (up to some constant)
    double GetLogPrior() const {
        double total = 0;
        if (blmode < 2) {
            total += BranchLengthsLogPrior();
        }
        if (nucmode < 2)    {
            total += NucRatesLogPrior();
        }
        total += MaskHyperLogPrior();
        total += MaskLogPrior();
        total += AAHyperLogPrior();
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
	double LambdaHyperLogPrior() const {
		return -lambda / 10;
	}

    //! log prior over branch lengths (iid exponential of rate lambda)
	double BranchLengthsLogPrior() const {
		double ret = branchlength->GetLogProb();
        if (blmode == 0)    {
            ret += LambdaHyperLogPrior();
        }
        return ret;
	}

    //! log prior over omega (gamma of mean omegahypermean and inverse shape omegahyperinvshape)
	double OmegaLogPrior() const {
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
		return alpha * log(beta) - Random::logGamma(alpha) + (alpha-1) * log(omega) - beta*omega;
	}

    //! log prior over nuc rates rho and pi (uniform)
    double NucRatesLogPrior() const {
        double total = 0;
        total += Random::logDirichletDensity(nucrelrate,nucrelratehypercenter,1.0/nucrelratehyperinvconc);
        total += Random::logDirichletDensity(nucstat,nucstathypercenter,1.0/nucstathyperinvconc);
        return total;
    }

    double MaskHyperLogPrior() const    {
        return 0;
    }

    double MaskLogPrior() const {
        return sitemaskarray->GetLogProb();
    }

    double MaskLogPrior(int i) const {
        return sitemaskarray->GetLogProb(i);
    }

    double AAHyperLogPrior() const  {
        double tot = 0;
        // uniform over epsilon
        // uniform over aacenter
        // exponential of mean 1 over concentration
        tot -= aaconcentration / Naa;
        return tot;
    }

    //! log prior of amino-acid fitness profiles
    double AALogPrior() const {
        return siteaafitnessarray->GetLogProb();
    }

    //! log prior of amino-acid fitness profile i
    double AALogPrior(int i) const {
        return siteaafitnessarray->GetLogProb(i);
    }

    //-------------------
    // Suff Stat and suffstatlogprobs
    //-------------------

    //! return log prob of the current substitution mapping, as a function of the current codon substitution process
	double PathSuffStatLogProb() const {
        return sitepathsuffstatarray->GetLogProb(*sitecodonmatrixarray);
	}

    //! return log prob of the substitution mappings for site i
    double PathSuffStatLogProb(int i) const {
        return sitepathsuffstatarray->GetVal(i).GetLogProb(sitecodonmatrixarray->GetVal(i));
    }

    //! return log prob of current branch lengths, as a function of branch lengths hyperparameter lambda
	double LambdaHyperSuffStatLogProb() const {
		return hyperlengthsuffstat.GetLogProb(1.0,lambda);
	}

    // when moving aaconcentration or aacenter or epsilon
    double AALogProb() const    {
        return AAHyperLogPrior() + AALogPrior();
    }

    double MaskLogProb() const  {
        return MaskHyperLogPrior() + MaskLogPrior();
    }

    //-------------------
    //  Log probs for MH moves
    //-------------------

    //! log prob factor to be recomputed when moving branch lengths hyperparameter lambda
    double LambdaHyperLogProb() const {
        return LambdaHyperLogPrior() + LambdaHyperSuffStatLogProb();
    }

    //! log prob factor to be recomputed when moving nucleotide mutation rate parameters (nucrelrate and nucstat)
    double NucRatesLogProb() const {
        return NucRatesLogPrior() + PathSuffStatLogProb();
    }

    //-------------------
    //  Collecting Suff Stats
    //-------------------

    //! collect sufficient statistics if substitution mappings across sites
	void CollectSitePathSuffStat()	{
		sitepathsuffstatarray->Clear();
        sitepathsuffstatarray->AddSuffStat(*phyloprocess);
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

            if (blmode < 2) {
                MoveBranchLengths();
            }

			CollectSitePathSuffStat();

            if (nucmode < 2)    {
                MoveNucRates();
            }

            if (omegamode < 2)  {
                MoveOmega();
            }

            MoveAAMixture(3);
		}
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
        ScalingMove(lambda,1.0,10,&AAMutSelMaskOmegaModel::LambdaHyperLogProb,&AAMutSelMaskOmegaModel::NoUpdate,this);
        ScalingMove(lambda,0.3,10,&AAMutSelMaskOmegaModel::LambdaHyperLogProb,&AAMutSelMaskOmegaModel::NoUpdate,this);
        blhypermean->SetAllBranches(1.0/lambda);
	}

    //! MH move on omega
	void MoveOmega()	{

		omegapathsuffstat.Clear();
		omegapathsuffstat.AddSuffStat(*sitecodonmatrixarray,*sitepathsuffstatarray);
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
		omega = Random::GammaSample(alpha + omegapathsuffstat.GetCount(), beta + omegapathsuffstat.GetBeta());
		UpdateCodonMatrices();
	}

    //! MH move on nucleotide rate parameters 
	void MoveNucRates()	{

        ProfileMove(nucrelrate,0.1,1,3,&AAMutSelMaskOmegaModel::NucRatesLogProb,&AAMutSelMaskOmegaModel::UpdateMatrices,this);
        ProfileMove(nucrelrate,0.03,3,3,&AAMutSelMaskOmegaModel::NucRatesLogProb,&AAMutSelMaskOmegaModel::UpdateMatrices,this);
        ProfileMove(nucrelrate,0.01,3,3,&AAMutSelMaskOmegaModel::NucRatesLogProb,&AAMutSelMaskOmegaModel::UpdateMatrices,this);

        ProfileMove(nucstat,0.1,1,3,&AAMutSelMaskOmegaModel::NucRatesLogProb,&AAMutSelMaskOmegaModel::UpdateMatrices,this);
        ProfileMove(nucstat,0.01,1,3,&AAMutSelMaskOmegaModel::NucRatesLogProb,&AAMutSelMaskOmegaModel::UpdateMatrices,this);
	}

    // pi, epsilon, aacenter, aaconcentration, masks, fitnessprofiles
    //! MCMC module for the mixture amino-acid fitness profiles
    void MoveAAMixture(int nrep)    {
        for (int rep=0; rep<nrep; rep++)  {
            MoveAAProfiles(3);
            MoveMasks();
            MoveAAHyperParameters(10);
            MoveMaskHyperParameters(10);
        }
        UpdateCodonMatrices();
    }

    void MoveAAHyperParameters(int nrep)    {
        for (int rep=0; rep<nrep; rep++)  {
            ProfileMove(aacenter,0.1,1,10,&AAMutSelMaskOmegaModel::AALogProb,&AAMutSelMaskOmegaModel::UpdateAA,this);
            ProfileMove(aacenter,0.03,3,10,&AAMutSelMaskOmegaModel::AALogProb,&AAMutSelMaskOmegaModel::UpdateAA,this);
            ProfileMove(aacenter,0.01,3,10,&AAMutSelMaskOmegaModel::AALogProb,&AAMutSelMaskOmegaModel::UpdateAA,this);
            ScalingMove(aaconcentration,1.0,10,&AAMutSelMaskOmegaModel::AALogProb,&AAMutSelMaskOmegaModel::UpdateAA,this);
            ScalingMove(aaconcentration,0.3,10,&AAMutSelMaskOmegaModel::AALogProb,&AAMutSelMaskOmegaModel::UpdateAA,this);
            SlidingMove(epsilon,1.0,10,0,epsilon0,&AAMutSelMaskOmegaModel::AALogProb,&AAMutSelMaskOmegaModel::UpdateAA,this);
            SlidingMove(epsilon,0.1,10,0,epsilon0,&AAMutSelMaskOmegaModel::AALogProb,&AAMutSelMaskOmegaModel::UpdateAA,this);
        }
    }

    void MoveMaskHyperParameters(int nrep)  {
        for (int rep=0; rep<nrep; rep++)  {
            SlidingMove(pi,1.0,10,0.05,0.975,&AAMutSelMaskOmegaModel::MaskLogProb,&AAMutSelMaskOmegaModel::UpdateMask,this);
            SlidingMove(pi,0.1,10,0.05,0.975,&AAMutSelMaskOmegaModel::MaskLogProb,&AAMutSelMaskOmegaModel::UpdateMask,this);
        }
    }

    double MoveMasks()    {
		double nacc = 0;
		double ntot = 0;
        for (int i=0; i<Nsite; i++) {
            vector<int>& mask = (*sitemaskarray)[i];
            int naa = 0;
            for (int k=0; k<Naa; k++)   {
                naa += mask[k];
            }
            for (int k=0; k<Naa; k++)   {
                if ((!mask[k]) || (naa > 1))    {
                    double deltalogprob = -AALogPrior(i) - MaskLogPrior(i);
                    naa -= mask[k];
                    mask[k] = 1-mask[k];
                    naa += mask[k];
                    UpdateAA(i);
                    deltalogprob += AALogPrior(i) + MaskLogPrior(i);
                    int accepted = (log(Random::Uniform()) < deltalogprob);
                    if (accepted)	{
                        nacc ++;
                    }
                    else	{
                        naa -= mask[k];
                        mask[k] = 1-mask[k];
                        naa += mask[k];
                        UpdateAA(i);
                    }
                    ntot++;
                }
            }
        }
		return nacc/ntot;
	}

    //! MH move on amino-acid fitness profiles
    void MoveAAProfiles(int nrep)   {
        CompMoveAAProfiles(nrep);
        MulMoveAAProfiles(nrep);
    }

    //! MH move on amino-acid fitness profiles: additive compensated move on pairs of entries of the vector
    double CompMoveAAProfiles(int nrep) {
        MoveAA(1.0,1,nrep);
        MoveAA(0.1,3,nrep);
        return 1.0;
    }

    //! MH move on amino-acid fitness profiles: multiplicative move (using the Gamma representation of the Dirichlet)
    double MulMoveAAProfiles(int nrep) {
        MoveAAGamma(3.0,nrep);
        MoveAAGamma(1.0,nrep);
        return 1.0;
    }

    //! MH move on amino-acid fitness profiles: additive compensated move on pairs of entries of the vector
	double MoveAA(double tuning, int n, int nrep)	{
		double nacc = 0;
		double ntot = 0;
		double bk[Naa];
        for (int i=0; i<Nsite; i++) {
                vector<double>& aa = (*siteaafitnessarray)[i];
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
		return nacc/ntot;
	}

    //! helper function: log density of 20 gammas
    double GammaAALogPrior(const vector<double>& x, const vector<double>& weight)   {
        double total = 0;
        for (int l=0; l<Naa; l++)   {
            total += (weight[l] -1)*log(x[l]) - x[l] - Random::logGamma(weight[l]);
        }
        return total;
    }

    //! MH move on amino-acid fitness profiles: multiplicative move (using the Gamma representation of the Dirichlet)
	double MoveAAGamma(double tuning, int nrep)	{

		double nacc = 0;
		double ntot = 0;
        for (int i=0; i<Nsite; i++) {

                vector<double>& aa = (*siteaafitnessarray)[i];
                const vector<double>& weight = siteweightarray->GetVal(i);
                double conc = 0;
                for (int k=0; k<Naa; k++)   {
                    conc += weight[k];
                }
                vector<double> x(Naa,0);
                double z = Random::sGamma(conc);
                for (int l=0; l<Naa; l++)   {
                    x[l] = z*aa[l];
                }

                double bkz = z;
                vector<double> bkx = x;
                vector<double> bkaa = aa;

                for (int rep=0; rep<nrep; rep++)	{

                    double deltalogprob = -GammaAALogPrior(x,weight) - PathSuffStatLogProb(i);

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

                    deltalogprob += GammaAALogPrior(x,weight) + PathSuffStatLogProb(i);

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
        return nacc/ntot;
	}

    //-------------------
    // Traces and Monitors
    // ------------------

	void TraceHeader(ostream& os) const override {
		os << "#logprior\tlnL\tlength\t";
		os << "omega\t";
        os << "aaent\t";
        os << "aaconc\t";
        os << "pi\t";
        os << "epsilon\t";
        os << "width\t";
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
        os << siteaafitnessarray->GetMeanEntropy() << '\t';
		os << aaconcentration << '\t';
        os << pi << '\t';
        os << epsilon << '\t';
        os << sitemaskarray->GetMeanWidth() << '\t';
        os << Random::GetEntropy(aacenter) << '\t';
		os << Random::GetEntropy(nucstat) << '\t';
		os << Random::GetEntropy(nucrelrate) << '\n';
	}

	void Monitor(ostream& os) const override {
        os << GetNsite() << '\t' << Naa << '\n';
        for (int i=0; i<GetNsite(); i++)    {
            os << i;
            for (int a=0; a<Naa; a++)   {
                os << '\t' << siteaafitnessarray->GetVal(i)[a];
            }
            os << '\n';
        }
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
        if (omegamode < 2)  {
            is >> omega;
        }
        is >> pi;
        is >> epsilon;
        is >> aacenter;
        is >> aaconcentration;
        is >> *sitemaskarray;
        is >> *siteaafitnessarray;
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
        if (omegamode < 2)  {
            os << omega << '\t';
        }
        os << pi << '\t';
        os << epsilon << '\t';
        os << aacenter << '\t';
        os << aaconcentration << '\t';
        os << *sitemaskarray;
        os << *siteaafitnessarray;
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
        if (omegamode < 2)  {
            size++;
        }
        size +=3 + Naa;
        size += Naa*Nsite*2;
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
        if (omegamode < 2)  {
            is >> omega;
        }
        is >> pi;
        is >> epsilon;
        is >> aacenter;
        is >> aaconcentration;
        is >> *sitemaskarray;
        is >> *siteaafitnessarray;
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
        if (omegamode < 2)  {
            os << omega;
        }
        os << pi;
        os << epsilon;
        os << aacenter;
        os << aaconcentration;
        os << *sitemaskarray;
        os << *siteaafitnessarray;
    }
};


