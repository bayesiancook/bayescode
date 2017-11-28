
#include "CodonSequenceAlignment.hpp"
#include "Tree.hpp"
#include "ProbModel.hpp"
#include "GTRSubMatrix.hpp"
#include "AAMutSelOmegaCodonSubMatrix.hpp"
#include "PhyloProcess.hpp"
#include "IIDGamma.hpp"
#include "IIDDirichlet.hpp"
#include "CodonSuffStat.hpp"
#include "ProbModel.hpp"
#include "Chrono.hpp"

class AAMutSelOmegaModel : public ProbModel {

	Tree* tree;
	FileSequenceAlignment* data;
	const TaxonSet* taxonset;
	CodonSequenceAlignment* codondata;

	int Nsite;
	int Ntaxa;
	int Nbranch;

	double lambda;
	BranchIIDGamma* branchlength;
	PoissonSuffStatBranchArray* lengthsuffstatarray;
	GammaSuffStat lambdasuffstat;

	std::vector<double> nucstat;
	std::vector<double> nucrelrate;
	GTRSubMatrix* nucmatrix;

    // of mean omegahypermean and inverse shape parameter omegahyperinvshape
    double omegahypermean;
    double omegahyperinvshape;
	double omega;
	OmegaSuffStat omegasuffstat;
	
    vector<double> aacenter;
    double aaconc;
    // double aainvconc;
    IIDDirichlet* aafitnessarray;
	AAMutSelOmegaCodonSubMatrixArray* codonmatrixarray;
    DirichletSuffStat aahypersuffstat;

	PhyloProcess* phyloprocess;

	PathSuffStatArray* pathsuffstatarray;

    Chrono aachrono;
    Chrono totchrono;
    double acc1,acc2,acc3,acc4,acc5;
    double tot1,tot2,tot3,tot4,tot5;

	public:

    //-------------------
    // Construction and allocation
    // ------------------

	AAMutSelOmegaModel(string datafile, string treefile) : aahypersuffstat(20) {

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

		// Allocate();
        totchrono.Reset();
        aachrono.Reset();
        acc1=acc2=acc3=acc4=acc5=0;
        tot1=tot2=tot3=tot4=tot5=0;
	}

    void Unfold()   {

		cerr << "-- unfold\n";
		phyloprocess->Unfold();
		cerr << phyloprocess->GetLogProb() << '\n';
		std::cerr << "-- mapping substitutions\n";
		phyloprocess->ResampleSub();
		// Trace(cerr);
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

        aacenter.assign(Naa,1.0/Naa);
        aaconc = ((double) Naa);
        // aainvconc = 1.0/Naa;
        aafitnessarray = new IIDDirichlet(Nsite,aacenter,aaconc);
        // aafitnessarray = new IIDDirichlet(Nsite,aacenter,1.0/aainvconc);

        omegahypermean = 1.0;
        omegahyperinvshape = 1.0;
		omega = 1.0;

        codonmatrixarray = new AAMutSelOmegaCodonSubMatrixArray(GetCodonStateSpace(), nucmatrix, aafitnessarray, omega);
		phyloprocess = new PhyloProcess(tree,codondata,branchlength,0,codonmatrixarray);
		pathsuffstatarray = new PathSuffStatArray(Nsite);
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

    const PoissonSuffStatBranchArray* GetLengthSuffStatArray() const {
        return lengthsuffstatarray;
    }

    //-------------------
    // Setting and updating
    // ------------------

    void SetBranchLengths(const ConstBranchArray<double>& inbranchlength)    {
        branchlength->Copy(inbranchlength);
    }

    void SetOmega(double inomega)   {
        omega = inomega;
        UpdateCodonMatrices();
    }

    void SetOmegaHyperParameters(double inomegahypermean, double inomegahyperinvshape)   {
        omegahypermean = inomegahypermean;
        omegahyperinvshape = inomegahyperinvshape;
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

	void UpdateCodonMatrices()	{
        codonmatrixarray->SetOmega(omega);
		codonmatrixarray->UpdateCodonMatrices();
	}

    void UpdateCodonMatrix(int site)    {
        (*codonmatrixarray)[site].CorruptMatrix();
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
        total += AAHyperLogPrior();
        total += AALogPrior();
        total += OmegaLogPrior();
        return total;
    }

	double GetLogLikelihood() const {
		return phyloprocess->GetLogProb();
	}

    double GetLogProb() const   {
        return GetLogPrior() + GetLogLikelihood();
    }

	double BranchLengthsHyperLogPrior() const {
		return -lambda / 10;
	}

	double BranchLengthsLogPrior() const {
		return branchlength->GetLogProb();
	}

	// exponential of mean 1
	double OmegaLogPrior() const {
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
		return alpha * log(beta) - Random::logGamma(alpha) + (alpha-1) * log(omega) - beta*omega;
	}

    double AAHyperLogPrior() const {
        return Random::logGammaDensity(aaconc,20.0,1.0);
        // return -aainvconc;
    }

    double AALogPrior() const {
        return aafitnessarray->GetLogProb();
    }

    double AALogPrior(int i) const {
        return aafitnessarray->GetLogProb(i);
    }

    double NucRatesLogPrior() const {
        return 0;
    }

    //-------------------
    // Suff Stat and suffstatlogprobs
    //-------------------

	double PathSuffStatLogProb() const {
		return pathsuffstatarray->GetLogProb(*codonmatrixarray);
	}

    double PathSuffStatLogProb(int site) const {
        return pathsuffstatarray->GetVal(site).GetLogProb(codonmatrixarray->GetVal(site));
    }

	double BranchLengthsHyperSuffStatLogProb() const {
		return lambdasuffstat.GetLogProb(1.0,lambda);
	}

    double AAHyperSuffStatLogProb() const   {
        return aahypersuffstat.GetLogProb(aacenter,aaconc);
        // return aahypersuffstat.GetLogProb(aacenter,1.0/aainvconc);
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
        return NucRatesLogPrior() + PathSuffStatLogProb();
    }

    // for moving aa hyper params (aacenter and aainvconc)
    double AAHyperLogProb() const   {
        return AAHyperLogPrior() + AAHyperSuffStatLogProb();
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

            totchrono.Start();
			ResampleBranchLengths();
			MoveBranchLengthsHyperParameter();

			CollectPathSuffStat();

            aachrono.Start();
            MoveAA();
            aachrono.Stop();
            MoveAAHyperParameters();
			MoveOmega();
			MoveNucRates();
            totchrono.Stop();
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
        ScalingMove(lambda,1.0,10,&AAMutSelOmegaModel::BranchLengthsHyperLogProb,&AAMutSelOmegaModel::NoUpdate,this);
        ScalingMove(lambda,0.3,10,&AAMutSelOmegaModel::BranchLengthsHyperLogProb,&AAMutSelOmegaModel::NoUpdate,this);
		branchlength->SetScale(lambda);
	}

	void CollectPathSuffStat()	{

		pathsuffstatarray->Clear();
		phyloprocess->AddPathSuffStat(*pathsuffstatarray);
	}

	void MoveOmega()	{

		omegasuffstat.Clear();
		omegasuffstat.AddSuffStat(*codonmatrixarray,*pathsuffstatarray);
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
		omega = Random::GammaSample(alpha + omegasuffstat.GetCount(), beta + omegasuffstat.GetBeta());
		UpdateCodonMatrices();
	}

	void MoveNucRates()	{

        UpdateMatrices();

        ProfileMove(nucrelrate,0.1,1,3,&AAMutSelOmegaModel::NucRatesLogProb,&AAMutSelOmegaModel::UpdateMatrices,this);
        ProfileMove(nucrelrate,0.03,3,3,&AAMutSelOmegaModel::NucRatesLogProb,&AAMutSelOmegaModel::UpdateMatrices,this);
        ProfileMove(nucrelrate,0.01,3,3,&AAMutSelOmegaModel::NucRatesLogProb,&AAMutSelOmegaModel::UpdateMatrices,this);

        ProfileMove(nucstat,0.1,1,3,&AAMutSelOmegaModel::NucRatesLogProb,&AAMutSelOmegaModel::UpdateMatrices,this);
        ProfileMove(nucstat,0.01,1,3,&AAMutSelOmegaModel::NucRatesLogProb,&AAMutSelOmegaModel::UpdateMatrices,this);

        UpdateMatrices();
	}

    void MoveAAHyperParameters()    {
        aahypersuffstat.Clear();
        aafitnessarray->AddSuffStat(aahypersuffstat);
        for (int rep=0; rep<10; rep++)  {
            ProfileMove(aacenter,0.1,1,10,&AAMutSelOmegaModel::AAHyperLogProb,&AAMutSelOmegaModel::NoUpdate,this);
            ProfileMove(aacenter,0.03,3,10,&AAMutSelOmegaModel::AAHyperLogProb,&AAMutSelOmegaModel::NoUpdate,this);
            ProfileMove(aacenter,0.01,3,10,&AAMutSelOmegaModel::AAHyperLogProb,&AAMutSelOmegaModel::NoUpdate,this);
            ScalingMove(aaconc,1.0,10,&AAMutSelOmegaModel::AAHyperLogProb,&AAMutSelOmegaModel::NoUpdate,this);
            ScalingMove(aaconc,0.3,10,&AAMutSelOmegaModel::AAHyperLogProb,&AAMutSelOmegaModel::NoUpdate,this);
            /*
            ScalingMove(aainvconc,1.0,10,&AAMutSelOmegaModel::AAHyperLogProb,&AAMutSelOmegaModel::NoUpdate,this);
            ScalingMove(aainvconc,0.3,10,&AAMutSelOmegaModel::AAHyperLogProb,&AAMutSelOmegaModel::NoUpdate,this);
            */
        }
        aafitnessarray->SetCenter(aacenter);
        aafitnessarray->SetConcentration(aaconc);
        // aafitnessarray->SetConcentration(1.0/aainvconc);
    }

    double MoveAA() {
        acc1 += MoveAA(1.0,1,10);
        acc2 += MoveAA(1.0,3,10);
        acc3 += MoveAA(0.3,3,10);
        acc4 += MoveAA(0.1,3,10);
        acc5 += MoveAA(0.1,5,10);
        tot1++;
        tot2++;
        tot3++;
        tot4++;
        tot5++;
        return 1.0;
    }

	double MoveAA(double tuning, int n, int nrep)	{
		double nacc = 0;
		double ntot = 0;
		double bk[Naa];
        for (int i=0; i<Nsite; i++) {
            vector<double>& aa = (*aafitnessarray)[i];
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

    //-------------------
    // Traces and Monitors
    // ------------------

	void TraceHeader(std::ostream& os) const {
		os << "#logprior\tlnL\tlength\t";
		os << "omega\t";
        os << "aaent\t";
        os << "aaconc\t";
        os << "aacenter\t";
		os << "statent\t";
		os << "rrent\n";
	}

	void Trace(ostream& os) const {	
		os << GetLogPrior() << '\t';
		os << GetLogLikelihood() << '\t';
        // 3x: per coding site (and not per nucleotide site)
        os << 3*branchlength->GetTotalLength() << '\t';
		os << omega << '\t';
        os << aafitnessarray->GetMeanEntropy() << '\t';
        os << aaconc << '\t';
        os << Random::GetEntropy(aacenter) << '\t';
		os << Random::GetEntropy(nucstat) << '\t';
		os << Random::GetEntropy(nucrelrate) << '\n';
	}

	void Monitor(ostream& os) const {
        os << "prop time in aa moves: " << aachrono.GetTime() / totchrono.GetTime() << '\n';
        os << "acceptance rates:\n";
        os << acc1/tot1 << '\n';
        os << acc2/tot2 << '\n';
        os << acc3/tot3 << '\n';
        os << acc4/tot4 << '\n';
        os << acc5/tot5 << '\n';
    }

	void FromStream(istream& is) {}
	void ToStream(ostream& os) const {}

};


