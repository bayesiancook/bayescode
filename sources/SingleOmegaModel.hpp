
// this is a codon model, with the following structure
//
// - branch lengths iid Exponential of rate lambda
// - nucleotide relative exchangeabilities and stationaries are uniform Dirichlet
// - there is one single omega=dN/dS for all sites and across all branches
//
// prior over omega:
// omega ~ Gamma(omegahypermean,omegahyperinvshape)
// these 2 hyperparameters are fixed when this model is used in isolation
// on the other hand, this model can be used in a multigene context (see MultiGeneSingeOmegaModel)
// in which case omegahypermean and hyperinvshape are estimated across genes

#include "CodonSequenceAlignment.hpp"
#include "Tree.hpp"
#include "ProbModel.hpp"
#include "GTRSubMatrix.hpp"
#include "CodonSubMatrix.hpp"
#include "PhyloProcess.hpp"
#include "IIDGamma.hpp"
#include "CodonSuffStat.hpp"
#include "ProbModel.hpp"

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
	OmegaSuffStat omegasuffstat;

    // Poisson suffstats for substitution histories, as a function of branch lengths
	PoissonSuffStatBranchArray* lengthsuffstatarray;

    // suff stats branch lengths, as a function of their hyper parameter lambda
    // (bl are iid gamma, of scale parameter lambda)
	GammaSuffStat lambdasuffstat;

	public:

    //-------------------
    // Construction and allocation
    // ------------------

	SingleOmegaModel(string datafile, string treefile)  {

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

	void Allocate()	{

		lambda = 10;
		branchlength = new BranchIIDGamma(*tree,1.0,lambda);
		lengthsuffstatarray = new PoissonSuffStatBranchArray(*tree);

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
	}

    void Unfold()   {

		cerr << "-- unfold\n";
		phyloprocess->Unfold();
		cerr << phyloprocess->GetLogProb() << '\n';
		std::cerr << "-- mapping substitutions\n";
		phyloprocess->ResampleSub();
    }

    //-------------------
    // Accessors
    // ------------------

	CodonStateSpace* GetCodonStateSpace() const {
		return (CodonStateSpace*) codondata->GetStateSpace();
	}

    double GetOmega() const {
        return omega;
    }

    //-------------------
    // Setting and updating
    // ------------------

    void SetOmega(double inomega)   {
        omega = inomega;
        UpdateCodonMatrix();
    }

    void SetOmegaHyperParameters(double inomegahypermean, double inomegahyperinvshape)   {
        omegahypermean = inomegahypermean;
        omegahyperinvshape = inomegahyperinvshape;
    }

    void SetBranchLengths(const BranchSelector<double>& inbranchlength)    {
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

	void UpdateCodonMatrix()	{
		codonmatrix->SetOmega(omega);
		codonmatrix->CorruptMatrix();
	}
		
    void UpdateMatrices()   {
        UpdateNucMatrix();
        UpdateCodonMatrix();
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
        total += OmegaLogPrior();
        return total;
    }

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

	double OmegaLogPrior()	const {
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
        return Random::logGammaDensity(omega,alpha,beta);
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
		return pathsuffstat.GetLogProb(*codonmatrix);
	}

	double BranchLengthsHyperSuffStatLogProb()	const {
		return lambdasuffstat.GetLogProb(1.0,lambda);
	}

    double NucRatesSuffStatLogProb() const {
        return nucpathsuffstat.GetLogProb(*nucmatrix,*GetCodonStateSpace());
    }

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
        ScalingMove(lambda,1.0,10,&SingleOmegaModel::BranchLengthsHyperLogProb,&SingleOmegaModel::NoUpdate,this);
        ScalingMove(lambda,0.3,10,&SingleOmegaModel::BranchLengthsHyperLogProb,&SingleOmegaModel::NoUpdate,this);
		branchlength->SetScale(lambda);
	}

	void CollectPathSuffStat()	{

		pathsuffstat.Clear();
		phyloprocess->AddPathSuffStat(pathsuffstat);
	}

	void MoveOmega()	{

		omegasuffstat.Clear();
		omegasuffstat.AddSuffStat(*codonmatrix,pathsuffstat);
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
		omega = Random::GammaSample(alpha + omegasuffstat.GetCount(), beta + omegasuffstat.GetBeta());
		UpdateCodonMatrix();
	}

    void CollectNucPathSuffStat()    {
        UpdateMatrices();
        nucpathsuffstat.Clear();
        nucpathsuffstat.AddSuffStat(*codonmatrix,pathsuffstat);
    }

	void MoveNucRates()	{

        CollectNucPathSuffStat();

        ProfileMove(nucrelrate,0.1,1,3,&SingleOmegaModel::NucRatesLogProb,&SingleOmegaModel::UpdateNucMatrix,this);
        ProfileMove(nucrelrate,0.03,3,3,&SingleOmegaModel::NucRatesLogProb,&SingleOmegaModel::UpdateNucMatrix,this);
        ProfileMove(nucrelrate,0.01,3,3,&SingleOmegaModel::NucRatesLogProb,&SingleOmegaModel::UpdateNucMatrix,this);

        ProfileMove(nucstat,0.1,1,3,&SingleOmegaModel::NucRatesLogProb,&SingleOmegaModel::UpdateNucMatrix,this);
        ProfileMove(nucstat,0.01,1,3,&SingleOmegaModel::NucRatesLogProb,&SingleOmegaModel::UpdateNucMatrix,this);

        UpdateMatrices();
	}

    //-------------------
    // Traces and Monitors
    // ------------------

	void TraceHeader(std::ostream& os) const {
		os << "#logprior\tlnL\tlength\tlambda\t";
		os << "omega\t";
		os << "statent\t";
		os << "rrent\n";
	}

	void Trace(ostream& os) const {	
		os << GetLogPrior() << '\t';
		os << GetLogLikelihood() << '\t';
        os << branchlength->GetTotalLength() << '\t';
		os << lambda << '\t';
		os << omega << '\t';
		os << Random::GetEntropy(nucstat) << '\t';
		os << Random::GetEntropy(nucrelrate) << '\n';
	}

	void Monitor(ostream& os) const {}

	void ToStream(ostream& os) const {
        os << lambda << '\n';
        os << *branchlength << '\n';
        os << omega << '\n';
        os << nucrelrate << '\n';
        os << nucstat << '\n';
    }

	void FromStream(istream& is) {
        is >> lambda;
        is >> *branchlength;
        is >> omega;
        is >> nucrelrate;
        is >> nucstat;
    }

};


