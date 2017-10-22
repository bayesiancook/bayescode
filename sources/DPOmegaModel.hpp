
#include "CodonSequenceAlignment.hpp"
#include "Tree.hpp"
#include "ProbModel.hpp"
#include "GTRSubMatrix.hpp"
#include "CodonSubMatrixArray.hpp"
#include "PhyloProcess.hpp"
#include "IIDGamma.hpp"
#include "CodonSuffStat.hpp"
#include "ProbModel.hpp"
#include "DiscGamma.hpp"
#include "MultinomialAllocationVector.hpp"
#include "StickBreakingProcess.hpp"

const int Nrr = Nnuc * (Nnuc-1) / 2;
const int Nstate = 61;

class DPOmegaModel : public ProbModel {

    // tree and data
	Tree* tree;
	FileSequenceAlignment* data;
	const TaxonSet* taxonset;
	CodonSequenceAlignment* codondata;

	int Nsite;
	int Ntaxa;
	int Nbranch;
    int Ncat; // number of categories of the stick breaking

    // branch lengths iid expo (gamma of shape 1 and scale lambda)
    // where lambda is a hyperparameter
	double lambda;
	BranchIIDGamma* branchlength;

    // nucleotide exchange rates and equilibrium frequencies (stationary probabilities)
	std::vector<double> nucrelrate;
	std::vector<double> nucstat;

    // a nucleotide matrix (parameterized by nucrelrate and nucstat)
	GTRSubMatrix* nucmatrix;

    // omega across sites: Dirichlet process
    // made in three steps

    // (1) component weights (truncated stick breaking)
    double kappa;
    StickBreakingProcess* weight;

    // (2) component values: independent omega values from Gamma distribution
    // of mean omegamean and inverse shape parameter omegainvshape
    double omegamean;
    double omegainvshape;
    IIDGamma* componentomegaarray;

    // (3) site allocations: multinomial given the weights
    // multinomial allocation of sites to components of omega distribution
	MultinomialAllocationVector* sitealloc;

    // an array of codon matrices (one for each discrete value of omega)
	MGOmegaCodonSubMatrixArray* componentcodonmatrixarray;
	
    // 2 arrays of matrix pointers across sites
    // obtained from component matrix array and site allocations
    // the 2 arrays are identical in value, they just differ by their type

	// this one is used by PhyloProcess: has to be a ConstArray<SubMatrix>
	ConstMixtureArray<SubMatrix>* sitesubmatrixarray;

	// this one is used for collecting omega suffstats: need to have access to the *codon* matrix for each site
	ConstMixtureArray<MGOmegaCodonSubMatrix>* sitecodonmatrixarray;

	PhyloProcess* phyloprocess;

	// suffstats

    // generic suff stats for substitution paths 
    // per site
	PathSuffStatArray* sitepathsuffstatarray;
    // per component of the mixture
	PathSuffStatArray* componentpathsuffstatarray;

    // which can be collected across all sites and branches
    // and summarized in terms of 4x4 suff stats, as a function of nucleotide rates
    NucPathSuffStat nucpathsuffstat;

    // or, alternatively, collected as a simple Poisson suff stat, as a function of omega
    // per site
	OmegaSuffStatArray* siteomegasuffstatarray;
    // per component
	OmegaSuffStatArray* componentomegasuffstatarray;

    // for moving omegamean and invshape
    GammaSuffStat omegahypersuffstat;

    // Poisson suffstats for substitution histories, as a function of branch lengths
	PoissonSuffStatBranchArray* lengthsuffstatarray;

    // suff stats for branch lengths, as a function of their hyper parameter lambda
    // (bl are iid gamma, of scale parameter lambda)
	GammaSuffStat lambdasuffstat;

	public:

    //-------------------
    // Construction and allocation
    // ------------------

	DPOmegaModel(string datafile, string treefile, int inNcat)  {

		data = new FileSequenceAlignment(datafile);
		codondata = new CodonSequenceAlignment(data, true);
        Ncat = inNcat;

		Nsite = codondata->GetNsite();    // # columns
		Ntaxa = codondata->GetNtaxa();

        if (Ncat == -1) {
            Ncat = Nsite;
            if (Ncat > 100)    {
                Ncat = 100;
            }
        }

		std::cerr << "-- Number of sites: " << Nsite << std::endl;

		taxonset = codondata->GetTaxonSet();

		// get tree from file (newick format)
		tree = new Tree(treefile);

		// check whether tree and data fits together
		tree->RegisterWith(taxonset);

		tree->SetIndices();
		Nbranch = tree->GetNbranch();
	}

	void Allocate()	{

		lambda = 10;
		branchlength = new BranchIIDGamma(*tree,1.0,lambda);

		nucrelrate.assign(Nrr,0);
        Random::DirichletSample(nucrelrate,vector<double>(Nrr,1.0/Nrr),((double) Nrr));

		nucstat.assign(Nnuc,0);
        Random::DirichletSample(nucstat,vector<double>(Nnuc,1.0/Nnuc),((double) Nnuc));

		nucmatrix = new GTRSubMatrix(Nnuc,nucrelrate,nucstat,true);

        kappa = 1.0;
        weight = new StickBreakingProcess(Ncat,kappa);

        sitealloc = new MultinomialAllocationVector(Nsite,weight->GetArray());

        omegamean = 1.0;
        omegainvshape = 1.0;
        componentomegaarray = new IIDGamma(Ncat,omegamean,omegainvshape);

        componentcodonmatrixarray = new MGOmegaCodonSubMatrixArray((CodonStateSpace*) codondata->GetStateSpace(),nucmatrix,componentomegaarray);

        sitesubmatrixarray = new ConstMixtureArray<SubMatrix>(componentcodonmatrixarray,sitealloc);
        sitecodonmatrixarray = new ConstMixtureArray<MGOmegaCodonSubMatrix>(componentcodonmatrixarray,sitealloc);

        phyloprocess = new PhyloProcess(tree,codondata,branchlength,0,sitesubmatrixarray);

        // suff stats
		lengthsuffstatarray = new PoissonSuffStatBranchArray(*tree);
        sitepathsuffstatarray = new PathSuffStatArray(Nsite);
        componentpathsuffstatarray = new PathSuffStatArray(Ncat);
        siteomegasuffstatarray = new OmegaSuffStatArray(Nsite);
        componentomegasuffstatarray = new OmegaSuffStatArray(Ncat);
	}

    void Unfold()   {

		cerr << "-- unfold\n";
		phyloprocess->Unfold();
		cerr << phyloprocess->GetLogProb() << '\n';
		std::cerr << "-- mapping substitutions\n";
		phyloprocess->ResampleSub();
    }

    void Update()   {
        UpdateMatrices();
    }

    //-------------------
    // Accessors
    // ------------------

	CodonStateSpace* GetCodonStateSpace() const {
		return (CodonStateSpace*) codondata->GetStateSpace();
	}

    //-------------------
    // Setting and updating
    // ------------------

    void SetBranchLengths(const ConstBranchArray<double>& inbranchlength)    {
        branchlength->Copy(inbranchlength);
    }

    void SetNucRates(const std::vector<double>& innucrelrate, const std::vector<double>& innucstat) {
        nucrelrate = innucrelrate;
        nucstat = innucstat;
        UpdateMatrices();
    }

	void UpdateNucMatrix()	{
		nucmatrix->CopyStationary(nucstat);
		nucmatrix->CorruptMatrix();
	}

    void UpdateCodonMatrices()  {
        componentcodonmatrixarray->UpdateCodonMatrices(sitealloc->GetOccupancies());
    }

    void UpdateMatrices()   {
        UpdateNucMatrix();
        UpdateCodonMatrices();
    }

    void NoUpdate() {}

    //-------------------
    // Priors and likelihood
    //-------------------

    double GetLogPrior() const {
        double total = 0;
        total += BranchLengthsHyperLogPrior();
        total += BranchLengthsLogPrior();
        total += NucRatesLogPrior();
        total += StickBreakingHyperLogPrior();
        total += StickBreakingLogPrior();
        total += OmegaHyperLogPrior();
        total += OmegaLogPrior();
        return total;
    }

    // conditional on site allocations
	double GetLogLikelihood() const {
		return phyloprocess->GetLogProb();
	}

    double GetLogProb() const   {
        return GetLogPrior() + GetLogLikelihood();
    }

	double BranchLengthsHyperLogPrior()	const {
        // exponential of mean 10
		return -lambda / 10;
	}

	double BranchLengthsLogPrior()	const {
		return branchlength->GetLogProb();
	}

    // uniform prior
    double NucRatesLogPrior() const {
        return 0;
    }

    double StickBreakingHyperLogPrior() const   {
        return -kappa/10;
    }

    double StickBreakingLogPrior() const    {
        return weight->GetLogProb(kappa);
        // return weight->GetMarginalLogProb(sitealloc->GetOccupancies());
    }

    double OmegaHyperLogPrior() const   {
        double total = 0;
        total -= omegamean;
        total -= omegainvshape;
        return total;
    }

    double OmegaLogPrior() const    {
        return componentomegaarray->GetLogProb();
    }

    //-------------------
    // Suff Stat and suffstatlogprobs
    //-------------------

    const PoissonSuffStatBranchArray* GetLengthSuffStatArray() const {
        return lengthsuffstatarray;
    }

    const NucPathSuffStat& GetNucPathSuffStat() const {
        return nucpathsuffstat;
    }

	double PathSuffStatLogProb() const {
        return componentpathsuffstatarray->GetLogProb(*componentcodonmatrixarray);
	}

	double BranchLengthsHyperSuffStatLogProb()	const {
		return lambdasuffstat.GetLogProb(1.0,lambda);
	}

    double NucRatesSuffStatLogProb() const {
        return nucpathsuffstat.GetLogProb(*nucmatrix,*GetCodonStateSpace());
    }

    double OmegaHyperSuffStatLogProb() const    {
        double alpha = 1.0 / omegainvshape;
        double beta = alpha / omegamean;
        return omegahypersuffstat.GetLogProb(alpha,beta);
    }

    /*
    double OmegaSuffStatLogProb() const {
        return componentomegaarray->GetLogProb(*componentomegasuffstatarray);
    }
    */

    //-------------------
    //  Log probs for MH moves
    //-------------------

    // for moving branch lengths hyperparameter lambda
    double BranchLengthsHyperLogProb() const {
        return BranchLengthsHyperLogPrior() + BranchLengthsHyperSuffStatLogProb();
    }

    // for moving nuc rates
    double NucRatesLogProb() const {
        return NucRatesLogPrior() + NucRatesSuffStatLogProb();
    }

    // for moving omegamean and invshape
    double OmegaHyperLogProb() const {
        return OmegaHyperLogPrior() + OmegaHyperSuffStatLogProb();
    }

    double StickBreakingHyperLogProb() const {
        return StickBreakingHyperLogPrior() + StickBreakingLogPrior();
    }

    //-------------------
    //  Moves 
    //-------------------

	double Move()	{
        ResampleSub(1.0);
        MoveParameters(30);
        return 1.0;
	}

    void ResampleSub(double frac)   {
        UpdateMatrices();
		phyloprocess->Move(frac);
    }

    void MoveParameters(int nrep)   {
		for (int rep=0; rep<nrep; rep++)	{

			ResampleBranchLengths();
			MoveBranchLengthsHyperParameter();

			CollectPathSuffStat();

			MoveOmega();
            MoveOmegaHyperParameters();
            MoveKappa();
			MoveNucRates();
		}
    }

    void CollectLengthSuffStat()    {
		lengthsuffstatarray->Clear();
		phyloprocess->AddLengthSuffStat(*lengthsuffstatarray);
    }

	void ResampleBranchLengths()	{
        CollectLengthSuffStat();
		branchlength->GibbsResample(*lengthsuffstatarray);
	}

	void MoveBranchLengthsHyperParameter()	{

		lambdasuffstat.Clear();
		branchlength->AddSuffStat(lambdasuffstat);
        ScalingMove(lambda,1.0,10,&DPOmegaModel::BranchLengthsHyperLogProb,&DPOmegaModel::NoUpdate,this);
        ScalingMove(lambda,0.3,10,&DPOmegaModel::BranchLengthsHyperLogProb,&DPOmegaModel::NoUpdate,this);
		branchlength->SetScale(lambda);
	}

    // per site
	void CollectPathSuffStat()	{
		sitepathsuffstatarray->Clear();
		phyloprocess->AddPathSuffStat(*sitepathsuffstatarray);
	}

    // per component of the mixture
	void CollectComponentPathSuffStat()	{
        componentpathsuffstatarray->Clear();
        sitepathsuffstatarray->AddToComponents(*componentpathsuffstatarray,*sitealloc);
    }

    void MoveOmega() {
        CollectSiteOmegaSuffStat();
        for (int rep=0; rep<10; rep++)  {
            ResampleOmega();
            ResampleAlloc();
            ResampleWeights();
        }
    }

    void CollectSiteOmegaSuffStat() {
        siteomegasuffstatarray->Clear();
        siteomegasuffstatarray->AddSuffStat(*sitecodonmatrixarray,*sitepathsuffstatarray);
    }

    void CollectComponentOmegaSuffStat()    {
        componentomegasuffstatarray->Clear();
        siteomegasuffstatarray->AddToComponents(*componentomegasuffstatarray,*sitealloc);
    }

    void ResampleOmega()    {
        CollectComponentOmegaSuffStat();
        componentomegaarray->GibbsResample(*componentomegasuffstatarray);
    }

    void ResampleAlloc()    {
        vector<double> postprob(Ncat,0);
        for (int i=0; i<Nsite; i++) {
            componentomegaarray->GetAllocPostProb(siteomegasuffstatarray->GetVal(i),weight->GetArray(),postprob);
            sitealloc->GibbsResample(i,postprob);
        }
        sitealloc->UpdateOccupancies();
    }

    void ResampleWeights()  {
        weight->GibbsResample(sitealloc->GetOccupancies());
    }

    void MoveOmegaHyperParameters() {

        omegahypersuffstat.Clear();
        componentomegaarray->AddSuffStat(omegahypersuffstat,sitealloc->GetOccupancies());
        ScalingMove(omegamean,1.0,10,&DPOmegaModel::OmegaHyperLogProb,&DPOmegaModel::NoUpdate,this);
        ScalingMove(omegamean,0.3,10,&DPOmegaModel::OmegaHyperLogProb,&DPOmegaModel::NoUpdate,this);
        ScalingMove(omegainvshape,1.0,10,&DPOmegaModel::OmegaHyperLogProb,&DPOmegaModel::NoUpdate,this);
        ScalingMove(omegainvshape,0.3,10,&DPOmegaModel::OmegaHyperLogProb,&DPOmegaModel::NoUpdate,this);
        double alpha = 1.0 / omegainvshape;
        double beta = alpha / omegamean;
        componentomegaarray->SetShape(alpha);
        componentomegaarray->SetScale(beta);
        componentomegaarray->PriorResample(sitealloc->GetOccupancies());
    }

    void MoveKappa()    {
        ScalingMove(kappa,1.0,10,&DPOmegaModel::StickBreakingHyperLogProb,&DPOmegaModel::NoUpdate,this);
        ScalingMove(kappa,0.3,10,&DPOmegaModel::StickBreakingHyperLogProb,&DPOmegaModel::NoUpdate,this);
    }

    void CollectNucPathSuffStat()    {
        UpdateMatrices();
        nucpathsuffstat.Clear();
        nucpathsuffstat.AddSuffStat(*componentcodonmatrixarray,*componentpathsuffstatarray);
    }

	void MoveNucRates()	{

        CollectComponentPathSuffStat();
        CollectNucPathSuffStat();

        ProfileMove(nucrelrate,0.1,1,3,&DPOmegaModel::NucRatesLogProb,&DPOmegaModel::UpdateNucMatrix,this);
        ProfileMove(nucrelrate,0.03,3,3,&DPOmegaModel::NucRatesLogProb,&DPOmegaModel::UpdateNucMatrix,this);
        ProfileMove(nucrelrate,0.01,3,3,&DPOmegaModel::NucRatesLogProb,&DPOmegaModel::UpdateNucMatrix,this);

        ProfileMove(nucstat,0.1,1,3,&DPOmegaModel::NucRatesLogProb,&DPOmegaModel::UpdateNucMatrix,this);
        ProfileMove(nucstat,0.01,1,3,&DPOmegaModel::NucRatesLogProb,&DPOmegaModel::UpdateNucMatrix,this);

        UpdateMatrices();
	}

    //-------------------
    // Traces and Monitors
    // ------------------

    int GetNcluster() const {

        int n = 0;
        const vector<int>& occupancy = sitealloc->GetOccupancies();
        for (int i=0; i<Ncat; i++)  {
            if (occupancy[i])    {
                n++;
            }
        }
        return n;
    }

	void TraceHeader(std::ostream& os) const {
		os << "#logprior\tlnL\tlength\tlambda\t";
        os << "ncluster\t";
		os << "omegamean\tinvshape\t";
		os << "statent\t";
		os << "rrent\n";
	}

	void Trace(ostream& os) const {	
		os << GetLogPrior() << '\t';
		os << GetLogLikelihood() << '\t';
        os << branchlength->GetTotalLength() << '\t';
		os << lambda << '\t';
        os << GetNcluster() << '\t';
		os << omegamean << '\t';
        os << omegainvshape << '\t';
		os << Random::GetEntropy(nucstat) << '\t';
		os << Random::GetEntropy(nucrelrate) << '\n';
        cerr << *componentomegaarray << '\n';
	}

	void Monitor(ostream& os) const {}

	void ToStream(ostream& os) const {
        os << lambda << '\n';
        os << *branchlength << '\n';
        os << kappa << '\n';
        os << omegamean << '\t' << omegainvshape << '\n';
        os << nucrelrate << '\n';
        os << nucstat << '\n';
    }

	void FromStream(istream& is) {
        is >> lambda;
        is >> *branchlength;
        is >> kappa;
        is >> omegamean >> omegainvshape;
        is >> nucrelrate;
        is >> nucstat;
    }

};


