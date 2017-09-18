
#include "CodonSequenceAlignment.hpp"
#include "Tree.hpp"
#include "ProbModel.hpp"
#include "GTRSubMatrix.hpp"
#include "CodonSubMatrix.hpp"
#include "PhyloProcess.hpp"
#include "IIDGamma.hpp"
#include "CodonSuffStat.hpp"
#include "CodonSubMatrixArray.hpp"
#include "M2aMix.hpp"
#include "MultinomialAllocationVector.hpp"
#include "Model.hpp"

const int Nrr = Nnuc * (Nnuc-1) / 2;
const int Nstate = 61;

class CodonM2aModel : public Model {

	public:

    //-------------------
    // Constructors
    // ------------------

	CodonM2aModel(string datafile, string treefile, double inpi);
	void Allocate();
    void Unfold();

    //-------------------
    // Accessors
    // ------------------

	int GetNsite() const {return codondata->GetNsite();}

    CodonStateSpace* GetCodonStateSpace() const {
		return (CodonStateSpace*) codondata->GetStateSpace();
    }

    double GetPurOm() const {
        return purom;
    }

    double GetDPosOm() const {
        return dposom;
    }

    double GetPurW() const {
        return purw;
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
    void SetAcrossGenesModes(int inblmode, int innucmode)   {
        blmode = inblmode;
        nucmode = innucmode;
    }

    void SetBranchLengths(const ConstBranchArray<double>& inbranchlength);
    void GetBranchLengths(BranchArray<double>& inbranchlength) const;

    void SetBranchLengthsHyperParameters(const ConstBranchArray<double>& inblmean, double inblinvshape);

    void SetNucRates(const std::vector<double>& innucrelrate, const std::vector<double>& innucstat);
    void GetNucRates(std::vector<double>& innucrelrate, std::vector<double>& innucstat) const;

    void SetNucRatesHyperParameters(const std::vector<double>& innucrelratehypercenter, double innucrelratehyperinvconc, const std::vector<double>& innucstathypercenter, double innucstathyperinvconc);

    void SetMixtureParameters(double inpurom, double indposom, double inpurw, double inposw);
    void GetMixtureParameters(double& inpurom, double& indposom, double& inpurw, double& inposw) const;

    void SetMixtureHyperParameters(double inpuromhypermean, double inpuromhyperinvconc, double indposomhypermean, double indposomhyperinvshape, double inpi, double inpurwhypermean, double inpurwhyperinvconc, double inposwhypermean, double inposwhyperinvconc);

    void NoUpdate() {}

    //-------------------
    // Matrices
    //-------------------

	void UpdateNucMatrix();
	void UpdateCodonMatrices();
	void UpdateMatrices();

    //-------------------
    // Traces and Monitors
    // ------------------

	double GetMeanOmega() const;

	void TraceHeader(std::ostream& os);
	void Trace(ostream& os);
	void TracePostProb(ostream& os);
    void GetSitesPostProb(double* array) const;

	void Monitor(ostream& os) {}
	void FromStream(istream& is) {}
	void ToStream(ostream& os) {}

    //-------------------
    // Likelihood
    //-------------------

	double GetLogLikelihood();
    double GetIntegratedLogLikelihood();

    //-------------------
    // Suff Stat and suffstatlogprobs
    //-------------------

    const PoissonSuffStatBranchArray* GetLengthSuffStatArray();
    const NucPathSuffStat& GetNucPathSuffStat();

	double PathSuffStatLogProb();
	double LambdaHyperSuffStatLogProb();
    double NucRatesSuffStatLogProb();
	double OmegaSuffStatLogProb();

    //-------------------
    // Priors
    //-------------------

	double GetLogPrior() const;

	double LambdaHyperLogPrior() const;
	double BranchLengthsLogPrior() const;
    double NucRatesLogPrior() const;

    double OmegaLogPrior() const;

    // Beta prior for purifmean
    double PurOmegaLogProb() const;
    // Gamma prior for dposom
	double PosOmegaLogProb() const;
    // Beta prior for purw
	double PurWeightLogProb() const;
    // mixture of point mass at 0 (with prob pi) and Beta distribution (with prob 1 - pi) for posw
	double PosWeightLogProb() const;
    // Bernoulli for whether posw == 0 or > 0
    double PosSwitchLogProb() const;

    //-------------------
    //  Log probs for MH moves
    //-------------------

    double LambdaHyperLogProb() {
        return LambdaHyperLogPrior() + LambdaHyperSuffStatLogProb();
    }

    double NucRatesLogProb()    {
        return NucRatesLogPrior() + NucRatesSuffStatLogProb();
    }

    double OmegaLogProb()   {
        return OmegaLogPrior() + OmegaSuffStatLogProb();
    }

    //-------------------
    //  Moves 
    //-------------------

	void Move();
    void ResampleSub(double frac);
    void MoveParameters(int nrep);

    //
    // Branch Lengths and hyperparam lambda
    //

    void MoveBranchLengths();
	void ResampleBranchLengths();
    void CollectLengthSuffStat();
	void MoveLambda();

	void CollectPathSuffStat();

    //
    // Omega mixture 
    //

	void CollectComponentPathSuffStat();
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

    // omega0 < 1: weight purw * (1 - posw)
    // omega1 = 1: weight (1-purw) * (1-posw)
    // omega2 > 1: weight posw

    // omega0 = purom
    double purom;

    // omega2 = 1 + dposom
	double dposom;

	double posw;
    double purw;

	M2aMix* componentomegaarray;
	MultinomialAllocationVector* sitealloc;
	vector<vector<double> > sitepostprobarray;

    // 
    // hyperparameters of the priors over the mixture parameters
    //

    // prior probability for the gene to be under positive selection (i.e. prior prob that posw > 0)
    double pi;

    // Beta prior for purom (with hypermean and hyper inverse concentration)
    double puromhypermean;
	double puromhyperinvconc;

    // Gamma prior for dposom = omega_pos - 1 (with hyper mean and inverse shape parameter)
    double dposomhypermean;
    double dposomhyperinvshape;

    // Beta prior for purw
    double purwhypermean;
    double purwhyperinvconc;

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

	// array of matrices across components of the mixture
	MGOmegaCodonSubMatrixArray* componentcodonmatrixarray;

	// arrays of matrices across sites (such as determined by the site allocations to the mixture components)
	// two versions differing only by their exact type

	// used for collecting omega suffstats: need to have access to the *codon* matrix for each site
	ConstMixtureArray<MGOmegaCodonSubMatrix>* sitecodonmatrixarray;

	// used by PhyloProcess: has to be a ConstArray<SubMatrix>
	ConstMixtureArray<SubMatrix>* sitesubmatrixarray;

	PhyloProcess* phyloprocess;

	// suffstats

	PoissonSuffStatBranchArray* lengthsuffstatarray;
	GammaSuffStat lambdasuffstat;
	OmegaSuffStatArray* siteomegasuffstatarray;
	PathSuffStatArray* sitepathsuffstatarray;
	PathSuffStatArray* componentpathsuffstatarray;

    NucPathSuffStat nucpathsuffstat;

    int blmode;
    int nucmode;


};

