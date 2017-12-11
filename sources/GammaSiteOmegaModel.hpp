
#include "CodonSequenceAlignment.hpp"
#include "Tree.hpp"
#include "ProbModel.hpp"
#include "GTRSubMatrix.hpp"
#include "CodonSubMatrixArray.hpp"
#include "PhyloProcess.hpp"
#include "IIDGamma.hpp"
#include "GammaSuffStat.hpp"
#include "CodonSuffStat.hpp"
#include "ProbModel.hpp"

class GammaSiteOmegaModel : public ProbModel {

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

    // omega across sites: Gamma distribution
    // of mean omegamean and inverse shape parameter omegainvshape
    double omegamean;
    double omegainvshape;

    IIDGamma* omegaarray;

    // a nucleotide matrix (parameterized by nucrelrate and nucstat)
	GTRSubMatrix* nucmatrix;
    // a codon matrix (parameterized by nucmatrix and omega)
	MGOmegaCodonSubMatrixArray* codonmatrixarray;
	
	PhyloProcess* phyloprocess;

	// suffstats

    // generic suff stats for substitution paths
	PathSuffStatArray* pathsuffstatarray;

    // which can be collected across all sites and branches
    // and summarized in terms of 4x4 suff stats, as a function of nucleotide rates
    NucPathSuffStat nucpathsuffstat;

    // or, alternatively, collected as a simple Poisson suff stat, as a function of omega
	OmegaPathSuffStatArray* omegapathsuffstatarray;

    // suffstat for resampling mean and invshape of distribution of omega's across sites
    GammaSuffStat omegahypersuffstat;

    // Poisson suffstats for substitution histories, as a function of branch lengths
	PoissonSuffStatBranchArray* lengthpathsuffstatarray;

    // suff stats for branch lengths, as a function of their hyper parameter lambda
    // (bl are iid gamma, of scale parameter lambda)
	GammaSuffStat hyperlengthsuffstat;

	public:

    //-------------------
    // Construction and allocation
    // ------------------

	GammaSiteOmegaModel(string datafile, string treefile)  {

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
		lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);

		nucrelrate.assign(Nrr,0);
        Random::DirichletSample(nucrelrate,vector<double>(Nrr,1.0/Nrr),((double) Nrr));

		nucstat.assign(Nnuc,0);
        Random::DirichletSample(nucstat,vector<double>(Nnuc,1.0/Nnuc),((double) Nnuc));

		nucmatrix = new GTRSubMatrix(Nnuc,nucrelrate,nucstat,true);

        omegamean = 1.0;
        omegainvshape = 1.0;
        double alpha = 1.0 / omegainvshape;
        double beta = alpha / omegamean;
        omegaarray = new IIDGamma(Nsite,alpha,beta);
		omegapathsuffstatarray = new OmegaPathSuffStatArray(Nsite);

		codonmatrixarray = new MGOmegaCodonSubMatrixArray(GetCodonStateSpace(), nucmatrix, omegaarray);

		phyloprocess = new PhyloProcess(tree,codondata,branchlength,0,codonmatrixarray);
		pathsuffstatarray = new PathSuffStatArray(Nsite);
	}

    void Unfold()   {

		cerr << "-- unfold\n";
		phyloprocess->Unfold();
		cerr << phyloprocess->GetLogLikelihood() << '\n';
		std::cerr << "-- mapping substitutions\n";
		phyloprocess->ResampleSub();
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

    void SetOmegaHyperParameters(double inomegamean, double inomegainvshape)   {
        omegamean = inomegamean;
        omegainvshape = inomegainvshape;
        double alpha = 1.0 / omegainvshape;
        double beta = alpha / omegamean;
        omegaarray->SetShape(alpha);
        omegaarray->SetScale(beta);
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

    void UpdateCodonMatrices()  {
        codonmatrixarray->UpdateCodonMatrices();
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
        total += OmegaHyperLogPrior();
        total += OmegaLogPrior();
        return total;
    }

	double GetLogLikelihood() const {
		return phyloprocess->GetLogLikelihood();
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

    double OmegaHyperLogPrior() const   {
        double total = 0;
        total -= omegamean;
        total -= omegainvshape;
        return total;
    }

	double OmegaLogPrior()	const {
        return omegaarray->GetLogProb();
	}

    //-------------------
    // Suff Stat and suffstatlogprobs
    //-------------------

    const PoissonSuffStatBranchArray* GetLengthPathSuffStatArray() const {
        return lengthpathsuffstatarray;
    }

    const NucPathSuffStat& GetNucPathSuffStat() const {
        return nucpathsuffstat;
    }

	double PathSuffStatLogProb() const {
		return pathsuffstatarray->GetLogProb(*codonmatrixarray);
	}

	double BranchLengthsHyperSuffStatLogProb()	const {
		return hyperlengthsuffstat.GetLogProb(1.0,lambda);
	}

    double NucRatesSuffStatLogProb() const {
        return nucpathsuffstat.GetLogProb(*nucmatrix,*GetCodonStateSpace());
    }

    double OmegaHyperSuffStatLogProb() const    {
        double alpha = 1.0 / omegainvshape;
        double beta = alpha / omegamean;
        return omegahypersuffstat.GetLogProb(alpha,beta);
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

    // for moving omegamean and omegainvshape
    double OmegaHyperLogProb() const    {
        return OmegaHyperLogPrior() + OmegaHyperSuffStatLogProb();
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
			MoveNucRates();
		}
    }

    void CollectLengthSuffStat()    {
		lengthpathsuffstatarray->Clear();
        lengthpathsuffstatarray->AddLengthPathSuffStat(*phyloprocess);
    }

	void ResampleBranchLengths()	{
        CollectLengthSuffStat();
		branchlength->GibbsResample(*lengthpathsuffstatarray);
	}

	void MoveBranchLengthsHyperParameter()	{

		hyperlengthsuffstat.Clear();
		hyperlengthsuffstat.AddSuffStat(*branchlength);
        ScalingMove(lambda,1.0,10,&GammaSiteOmegaModel::BranchLengthsHyperLogProb,&GammaSiteOmegaModel::NoUpdate,this);
        ScalingMove(lambda,0.3,10,&GammaSiteOmegaModel::BranchLengthsHyperLogProb,&GammaSiteOmegaModel::NoUpdate,this);
		branchlength->SetScale(lambda);
	}

	void CollectPathSuffStat()	{

		pathsuffstatarray->Clear();
        pathsuffstatarray->AddSuffStat(*phyloprocess);
	}

	void MoveOmega()	{

		omegapathsuffstatarray->Clear();
		omegapathsuffstatarray->AddSuffStat(*codonmatrixarray,*pathsuffstatarray);
		omegaarray->GibbsResample(*omegapathsuffstatarray);
		UpdateCodonMatrices();
	}

    void MoveOmegaHyperParameters() {

        omegahypersuffstat.Clear();
        omegahypersuffstat.AddSuffStat(*omegaarray);
        ScalingMove(omegamean,1.0,10,&GammaSiteOmegaModel::OmegaHyperLogProb,&GammaSiteOmegaModel::NoUpdate,this);
        ScalingMove(omegamean,0.3,10,&GammaSiteOmegaModel::OmegaHyperLogProb,&GammaSiteOmegaModel::NoUpdate,this);
        ScalingMove(omegainvshape,1.0,10,&GammaSiteOmegaModel::OmegaHyperLogProb,&GammaSiteOmegaModel::NoUpdate,this);
        ScalingMove(omegainvshape,0.3,10,&GammaSiteOmegaModel::OmegaHyperLogProb,&GammaSiteOmegaModel::NoUpdate,this);
        double alpha = 1.0 / omegainvshape;
        double beta = alpha / omegamean;
        omegaarray->SetShape(alpha);
        omegaarray->SetScale(beta);
    }

    void CollectNucPathSuffStat()    {
        UpdateMatrices();
        nucpathsuffstat.Clear();
        nucpathsuffstat.AddSuffStat(*codonmatrixarray,*pathsuffstatarray);
    }

	void MoveNucRates()	{

        CollectNucPathSuffStat();

        ProfileMove(nucrelrate,0.1,1,3,&GammaSiteOmegaModel::NucRatesLogProb,&GammaSiteOmegaModel::UpdateNucMatrix,this);
        ProfileMove(nucrelrate,0.03,3,3,&GammaSiteOmegaModel::NucRatesLogProb,&GammaSiteOmegaModel::UpdateNucMatrix,this);
        ProfileMove(nucrelrate,0.01,3,3,&GammaSiteOmegaModel::NucRatesLogProb,&GammaSiteOmegaModel::UpdateNucMatrix,this);

        ProfileMove(nucstat,0.1,1,3,&GammaSiteOmegaModel::NucRatesLogProb,&GammaSiteOmegaModel::UpdateNucMatrix,this);
        ProfileMove(nucstat,0.01,1,3,&GammaSiteOmegaModel::NucRatesLogProb,&GammaSiteOmegaModel::UpdateNucMatrix,this);

        UpdateMatrices();
	}

    //-------------------
    // Traces and Monitors
    // ------------------

    double GetEmpiricalPosFrac() const {
        double tot = 0;
        for (int i=0; i<Nsite; i++) {
            if ((*omegaarray)[i] > 1.0) {
                tot++;
            }
        }
        return tot / Nsite;
    }

	void TraceHeader(std::ostream& os) const {
		os << "#logprior\tlnL\tlength\t";
		os << "omegamean\tinvshape\t";
        os << "posfrac\t";
		os << "statent\t";
		os << "rrent\n";
	}

	void Trace(ostream& os) const {	
		os << GetLogPrior() << '\t';
		os << GetLogLikelihood() << '\t';
        os << branchlength->GetTotalLength() << '\t';
		os << omegamean << '\t';
        os << omegainvshape << '\t';
        os << GetEmpiricalPosFrac() << '\t';
		os << Random::GetEntropy(nucstat) << '\t';
		os << Random::GetEntropy(nucrelrate) << '\n';
	}

    void TraceSiteOmega(ostream& os) const  {
        for (int i=0; i<Nsite; i++) {
            os << (*omegaarray)[i] << '\t';
        }
        os << '\n';
        os.flush();
    }

	void Monitor(ostream& os) const {}

	void ToStream(ostream& os) const {
        os << lambda << '\n';
        os << *branchlength << '\n';
        os << omegamean << '\t' << omegainvshape << '\n';
        os << *omegaarray << '\n';
        os << nucrelrate << '\n';
        os << nucstat << '\n';
    }

	void FromStream(istream& is) {
        is >> lambda;
        is >> *branchlength;
        is >> omegamean >> omegainvshape;
        is >> *omegaarray;
        is >> nucrelrate;
        is >> nucstat;
    }

};


