
#include "CodonSequenceAlignment.hpp"
#include "Tree.hpp"
#include "ProbModel.hpp"
#include "GTRSubMatrix.hpp"
#include "CodonSubMatrix.hpp"
#include "AAMutSelOmegaCodonSubMatrix.hpp"
#include "PhyloProcess.hpp"
#include "IIDGamma.hpp"
#include "IIDDirichlet.hpp"
#include "CodonSuffStat.hpp"
#include "CodonSubMatrixArray.hpp"
#include "MechM2aMix.hpp"
#include "MultinomialAllocationVector.hpp"
#include "ProbModel.hpp"

class AAMutSelM2aModel : public ProbModel {
    
	public:

    //-------------------
    // Constructors
    // ------------------

	AAMutSelM2aModel(string datafile, string treefile, double inpi);
	void Allocate();
    void Unfold();

    //-------------------
    // Accessors
    // ------------------

	int GetNsite() const {return codondata->GetNsite();}

    CodonStateSpace* GetCodonStateSpace() const {
		return (CodonStateSpace*) codondata->GetStateSpace();
    }

    double GetDPosOm() const {
        return dposom;
    }

    double GetPosW() const {
        return posw;
    }

    bool FixedBranchLengths() const    {
        return blmode == 2;
    }

    bool FixedNucRates() const  {
        return nucmode == 2;
    }

    //-------------------
    // Setting and updating
    // ------------------

    // Setting model features and (hyper) parameters

    // called upon constructing the model
    // mode == 2: global
    // mode == 1: gene specific, with hyperparameters estimated across genes
    // mode == 0: gene-specific, with fixed hyperparameters
    //void SetAcrossGenesModes(int inblmode, int innucmode)   {
    void SetAcrossGenesModes(int inblmode)   {
        blmode = inblmode;
        //nucmode = innucmode;
    }

    void SetBranchLengths(const BranchSelector<double>& inbranchlength);
    void GetBranchLengths(BranchArray<double>& inbranchlength) const;

    void SetBranchLengthsHyperParameters(const BranchSelector<double>& inblmean, double inblinvshape);

    void SetNucRates(const std::vector<double>& innucrelrate, const std::vector<double>& innucstat);
    void GetNucRates(std::vector<double>& innucrelrate, std::vector<double>& innucstat) const;

    void SetNucRatesHyperParameters(const std::vector<double>& innucrelratehypercenter, double innucrelratehyperinvconc, const std::vector<double>& innucstathypercenter, double innucstathyperinvconc);

    //void SetMixtureParameters(double inpurom, double indposom, double inpurw, double inposw);
    //void GetMixtureParameters(double& inpurom, double& indposom, double& inpurw, double& inposw) const;
    void SetMixtureParameters(double indposom, double inposw);
    void GetMixtureParameters(double& indposom, double& inposw) const;

    //void SetMixtureHyperParameters(double inpuromhypermean, double inpuromhyperinvconc, double indposomhypermean, double indposomhyperinvshape, double inpi, double inpurwhypermean, double inpurwhyperinvconc, double inposwhypermean, double inposwhyperinvconc);
    void SetMixtureHyperParameters(double indposomhypermean, double indposomhyperinvshape, double inpi, double inposwhypermean, double inposwhyperinvconc);

    // void Update() override {}
    void NoUpdate() {}

    //-------------------
    // Matrices
    //-------------------

    void Update() override;
	void UpdateNucMatrix();
	void UpdateCodonMatrices();
    void UpdateCodonMatrix(int i);
	void UpdateMatrices();

    //-------------------
    // Traces and Monitors
    // ------------------

	double GetMeanOmega() const;

	void TraceHeader(std::ostream& os) const override;
	void Trace(ostream& os) const override;
	void TracePostProb(ostream& os) const ;
    void GetSitesPostProb(double* array) const;

	// void Monitor(ostream& os) const override {}
	void FromStream(istream& is) override;
	void ToStream(ostream& os) const override;

    //-------------------
    // Likelihood
    //-------------------

    double GetLogProb() const override {
        //return GetLogPrior() + GetIntegratedLogLikelihood();
        return GetLogPrior() + GetLogLikelihood();
    }

	double GetLogLikelihood() const;
    //double GetIntegratedLogLikelihood() const;

    //-------------------
    // Suff Stat and suffstatlogprobs
    //-------------------

    const PoissonSuffStatBranchArray* GetLengthPathSuffStatArray() const;
    const NucPathSuffStat& GetNucPathSuffStat() const;

	double PathSuffStatLogProb() const;
	double PathSuffStatLogProb(int site) const;
	double LambdaHyperSuffStatLogProb() const;
    //double NucRatesSuffStatLogProb() const;
	double OmegaSuffStatLogProb() const;
    double AAHyperSuffStatLogProb() const;
	
    void CollectPathSuffStat();

    //-------------------
    // Priors
    //-------------------

	double GetLogPrior() const;

	double LambdaHyperLogPrior() const;
	double BranchLengthsLogPrior() const;
    double NucRatesLogPrior() const;

    double OmegaLogPrior() const;

    // Beta prior for purifmean
    //double PurOmegaLogProb() const;
    // Gamma prior for dposom
	double PosOmegaLogProb() const;
    // Beta prior for purw
	double PurWeightLogProb() const;
    // mixture of point mass at 0 (with prob pi) and Beta distribution (with prob 1 - pi) for posw
	double PosWeightLogProb() const;
    // Bernoulli for whether posw == 0 or > 0
    double PosSwitchLogProb() const;

    double AAHyperLogPrior() const;
    double AALogPrior() const;
    double AALogPrior(int i) const;

    //-------------------
    //  Log probs for MH moves
    //-------------------

    double LambdaHyperLogProb() const {
        return LambdaHyperLogPrior() + LambdaHyperSuffStatLogProb();
    }

    double NucRatesLogProb() const {
        //return NucRatesLogPrior() + NucRatesSuffStatLogProb();
        return NucRatesLogPrior() + PathSuffStatLogProb();
    }

    double OmegaLogProb() const {
        return OmegaLogPrior() + OmegaSuffStatLogProb();
    }

    double AAHyperLogProb() const   {
        return AAHyperLogPrior() + AAHyperSuffStatLogProb();
    }

    //-------------------
    //  Moves 
    //-------------------

	double Move();
    void ResampleSub(double frac);
    void MoveParameters(int nrep);

    //
    // Branch Lengths and hyperparam lambda
    //

    void MoveBranchLengths();
	void ResampleBranchLengths();
    void CollectLengthSuffStat();
	void MoveLambda();


    //
    // Omega mixture 
    //

	//void CollectComponentPathSuffStat();
	void CollectSiteOmegaPathSuffStat();
	void MoveOmega();
	void CollectOmegaSuffStat();
	void ResampleAlloc();
    double DrawBetaPosWeight();
	double SwitchPosWeight(int nrep);

    //
    // nucleotide parameters
    //

    void CollectNucPathSuffStat();
	void MoveNucRates();

    //
    // AA fitness
    //

    void MoveAAHyperParameters();
    double MoveAA();
    double MoveAA(double tuning, int n, int rep);

    //-------------------
    // Data structures
    // ------------------

    private:

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
	
    //
    // parameters of the distribution of omega across sites
    //

    // omega1 = 1: weight * (1-posw)
    // omega2 > 1: weight posw

	double dposom;
	double posw;

	MechM2aMix* componentomegaarray;
    MixtureSelector<double>* siteomegaarray;
	MultinomialAllocationVector* sitealloc;
	mutable vector<vector<double> > sitepostprobarray;

    // 
    // hyperparameters of the priors over the mixture parameters
    //

    // prior probability for the gene to be under positive selection (i.e. prior prob that posw > 0)
    double pi;

    // Beta prior for purom (with hypermean and hyper inverse concentration)
    //double puromhypermean;
	//double puromhyperinvconc;

    // Gamma prior for dposom = omega_pos - 1 (with hyper mean and inverse shape parameter)
    double dposomhypermean;
    double dposomhyperinvshape;

    // Beta prior for purw
    //double purwhypermean;
    //double purwhyperinvconc;

    // Beta prior for posw (assuming posw>0)
    double poswhypermean;
    double poswhyperinvconc;

    // nucleotide rates hyperparameters
    vector<double> nucrelratehypercenter;
    double nucrelratehyperinvconc;
    vector<double> nucstathypercenter;
    double nucstathyperinvconc;

    // nucleotide rate parameters
	vector<double> nucrelrate;
	vector<double> nucstat;
	GTRSubMatrix* nucmatrix;

	// array of matrices across alignments
    vector<double> aacenter;
    double aainvconc;
    IIDDirichlet* aafitnessarray;
	// used for collecting omega suffstats: need to have access to the *codon* matrix for each site
	AAMutSelOmegaCodonSubMatrixArray* sitecodonmatrixarray;
    DirichletSuffStat aahypersuffstat;


	// used by PhyloProcess: has to be a Selector<SubMatrix>
	MixtureSelector<SubMatrix>* sitesubmatrixarray;

	PhyloProcess* phyloprocess;

	// suffstats

	PoissonSuffStatBranchArray* lengthpathsuffstatarray;
	GammaSuffStat lambdasuffstat;
	OmegaSuffStatArray* siteomegasuffstatarray;
	PathSuffStatArray* sitepathsuffstatarray;
	//PathSuffStatArray* componentpathsuffstatarray;

    NucPathSuffStat nucpathsuffstat;

    int blmode;
    int nucmode;


};

