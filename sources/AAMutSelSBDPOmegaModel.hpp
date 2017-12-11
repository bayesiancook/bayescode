
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

class AAMutSelSBDPOmegaModel : public ProbModel {

	Tree* tree;
	FileSequenceAlignment* data;
	const TaxonSet* taxonset;
	CodonSequenceAlignment* codondata;

	int Nsite;
	int Ntaxa;
	int Nbranch;

    int Ncat;

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
	
    // mixture over amino-acid fitness profiles

    // (1) component weights (truncated stick breaking)
    double kappa;
    StickBreakingProcess* weight;

    OccupancySuffStat* occupancy;

    // (2) component values: independent aa fitness profiles from Dirichlet distribution
    vector<double> aacenter;
    double aaconc;
    // double aainvconc;
    IIDDirichlet* componentaafitnessarray;
    DirichletSuffStat aahypersuffstat;

    // (3) site allocations: multinomial given the weights
    // multinomial allocation of sites to components of aa fitness profile distribution
	MultinomialAllocationVector* sitealloc;

    // an array of codon matrices (one for each distinct aa fitness profile)
	AAMutSelOmegaCodonSubMatrixArray* componentcodonmatrixarray;

	// this one is used by PhyloProcess: has to be a Selector<SubMatrix>
	MixtureSelector<SubMatrix>* sitesubmatrixarray;

	PhyloProcess* phyloprocess;

	PathSuffStatArray* sitepathsuffstatarray;
	PathSuffStatArray* componentpathsuffstatarray;

    int fixbl;

    Chrono aachrono;
    Chrono totchrono;

	public:

    //-------------------
    // Construction and allocation
    // ------------------

	AAMutSelSBDPOmegaModel(string datafile, string treefile, int inNcat) : aahypersuffstat(20) {

        fixbl = 0;

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


		std::cerr << "-- Number of sites: " << Nsite << std::endl;

		taxonset = codondata->GetTaxonSet();

		// get tree from file (newick format)
		tree = new Tree(treefile);

		// check whether tree and data fits together
		tree->RegisterWith(taxonset);

		tree->SetIndices();
		Nbranch = tree->GetNbranch();

		// Allocate();
	}

    void Unfold()   {

		cerr << "-- unfold\n";
		phyloprocess->Unfold();
		cerr << phyloprocess->GetLogLikelihood() << '\n';
		std::cerr << "-- mapping substitutions\n";
		phyloprocess->ResampleSub();
		// Trace(cerr);
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

        kappa = 1.0;
        weight = new StickBreakingProcess(Ncat,kappa);
        occupancy = new OccupancySuffStat(Ncat);

        sitealloc = new MultinomialAllocationVector(Nsite,weight->GetArray());

        aacenter.assign(Naa,1.0/Naa);
        aaconc = ((double) Naa);
        // aainvconc = 1.0/Naa;
        componentaafitnessarray = new IIDDirichlet(Ncat,aacenter,aaconc);
        // componentaafitnessarray = new IIDDirichlet(Ncat,aacenter,1.0/aainvconc);

        omegahypermean = 1.0;
        omegahyperinvshape = 1.0;
		omega = 1.0;

        componentcodonmatrixarray = new AAMutSelOmegaCodonSubMatrixArray(GetCodonStateSpace(), nucmatrix, componentaafitnessarray, omega);

        sitesubmatrixarray = new MixtureSelector<SubMatrix>(componentcodonmatrixarray,sitealloc);

		phyloprocess = new PhyloProcess(tree,codondata,branchlength,0,sitesubmatrixarray);
		sitepathsuffstatarray = new PathSuffStatArray(Nsite);
        componentpathsuffstatarray = new PathSuffStatArray(Ncat);

        Update();
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

    const PoissonSuffStatBranchArray* GetLengthPathSuffStatArray() const {
        return lengthpathsuffstatarray;
    }

    //-------------------
    // Setting and updating
    // ------------------

    void SetFixBL(int infixbl)    {
        fixbl = infixbl;
    }

    void SetBranchLengths(const BranchSelector<double>& inbranchlength)    {
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
        componentcodonmatrixarray->SetOmega(omega);
		componentcodonmatrixarray->UpdateCodonMatrices();
	}

    void UpdateCodonMatrix(int k)    {
        (*componentcodonmatrixarray)[k].CorruptMatrix();
    }
		
    void UpdateMatrices()   {
        UpdateNucMatrix();
        UpdateCodonMatrices();
    }

    void NoUpdate() {}

    void Update()   {
        UpdateOccupancies();
        UpdateMatrices();
    }

    //-------------------
    // Priors and likelihood
    //-------------------

    double GetLogPrior() const {
        double total = 0;
        if (! fixbl)    {
            total += BranchLengthsHyperLogPrior();
            total += BranchLengthsLogPrior();
        }
        total += NucRatesLogPrior();
        total += StickBreakingHyperLogPrior();
        total += StickBreakingLogPrior();
        total += AAHyperLogPrior();
        total += AALogPrior();
        total += OmegaLogPrior();
        return total;
    }

	double GetLogLikelihood() const {
		return phyloprocess->GetLogLikelihood();
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

    double StickBreakingHyperLogPrior() const   {
        return -kappa/10;
    }

    double StickBreakingLogPrior() const    {
        return weight->GetLogProb(kappa);
    }

    double AAHyperLogPrior() const {
        return Random::logGammaDensity(aaconc,20.0,1.0);
        // return -aainvconc;
    }

    double AALogPrior() const {
        return componentaafitnessarray->GetLogProb();
    }

    double AALogPrior(int k) const {
        return componentaafitnessarray->GetLogProb(k);
    }

    double NucRatesLogPrior() const {
        return 0;
    }

    //-------------------
    // Suff Stat and suffstatlogprobs
    //-------------------

	double PathSuffStatLogProb() const {
        return componentpathsuffstatarray->GetLogProb(*componentcodonmatrixarray);
	}

    double PathSuffStatLogProb(int k) const {
        return componentpathsuffstatarray->GetVal(k).GetLogProb(componentcodonmatrixarray->GetVal(k));
    }

	double BranchLengthsHyperSuffStatLogProb() const {
		return hyperlengthsuffstat.GetLogProb(1.0,lambda);
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

    // for moving aa hyper params (aacenter and aaconc)
    double AAHyperLogProb() const   {
        return AAHyperLogPrior() + AAHyperSuffStatLogProb();
    }

    // for moving kappa
    double StickBreakingHyperLogProb() const {
        return StickBreakingHyperLogPrior() + StickBreakingLogPrior();
    }

    //-------------------
    //  Collecting Suff Stats
    //-------------------

    // per site
	void CollectSitePathSuffStat()	{
		sitepathsuffstatarray->Clear();
        sitepathsuffstatarray->AddSuffStat(*phyloprocess);
	}

    // per component of the mixture
	void CollectComponentPathSuffStat()	{
        componentpathsuffstatarray->Clear();
        componentpathsuffstatarray->Add(*sitepathsuffstatarray,*sitealloc);
    }

    void CollectLengthSuffStat()    {
		lengthpathsuffstatarray->Clear();
        lengthpathsuffstatarray->AddLengthPathSuffStat(*phyloprocess);
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
            if (! fixbl)    {
                ResampleBranchLengths();
                MoveBranchLengthsHyperParameter();
            }

			CollectSitePathSuffStat();
            CollectComponentPathSuffStat();

			MoveNucRates();
			MoveOmega();

            aachrono.Start();
            MoveAAMixture();
            aachrono.Stop();
            MoveAAHyperParameters();
            MoveKappa();
            totchrono.Stop();
		}
	}

	void ResampleBranchLengths()	{
        CollectLengthSuffStat();
		branchlength->GibbsResample(*lengthpathsuffstatarray);
	}

	void MoveBranchLengthsHyperParameter()	{
		hyperlengthsuffstat.Clear();
		hyperlengthsuffstat.AddSuffStat(*branchlength);
        ScalingMove(lambda,1.0,10,&AAMutSelSBDPOmegaModel::BranchLengthsHyperLogProb,&AAMutSelSBDPOmegaModel::NoUpdate,this);
        ScalingMove(lambda,0.3,10,&AAMutSelSBDPOmegaModel::BranchLengthsHyperLogProb,&AAMutSelSBDPOmegaModel::NoUpdate,this);
		branchlength->SetScale(lambda);
	}

	void MoveOmega()	{

		omegapathsuffstat.Clear();
		omegapathsuffstat.AddSuffStat(*componentcodonmatrixarray,*componentpathsuffstatarray);
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
		omega = Random::GammaSample(alpha + omegapathsuffstat.GetCount(), beta + omegapathsuffstat.GetBeta());
		UpdateCodonMatrices();
	}

	void MoveNucRates()	{

        ProfileMove(nucrelrate,0.1,1,3,&AAMutSelSBDPOmegaModel::NucRatesLogProb,&AAMutSelSBDPOmegaModel::UpdateMatrices,this);
        ProfileMove(nucrelrate,0.03,3,3,&AAMutSelSBDPOmegaModel::NucRatesLogProb,&AAMutSelSBDPOmegaModel::UpdateMatrices,this);
        ProfileMove(nucrelrate,0.01,3,3,&AAMutSelSBDPOmegaModel::NucRatesLogProb,&AAMutSelSBDPOmegaModel::UpdateMatrices,this);

        ProfileMove(nucstat,0.1,1,3,&AAMutSelSBDPOmegaModel::NucRatesLogProb,&AAMutSelSBDPOmegaModel::UpdateMatrices,this);
        ProfileMove(nucstat,0.01,1,3,&AAMutSelSBDPOmegaModel::NucRatesLogProb,&AAMutSelSBDPOmegaModel::UpdateMatrices,this);
	}

    void MoveAAHyperParameters()    {

        aahypersuffstat.Clear();
        componentaafitnessarray->AddSuffStat(aahypersuffstat,*occupancy);
        for (int rep=0; rep<10; rep++)  {
            ProfileMove(aacenter,0.1,1,10,&AAMutSelSBDPOmegaModel::AAHyperLogProb,&AAMutSelSBDPOmegaModel::NoUpdate,this);
            ProfileMove(aacenter,0.03,3,10,&AAMutSelSBDPOmegaModel::AAHyperLogProb,&AAMutSelSBDPOmegaModel::NoUpdate,this);
            ProfileMove(aacenter,0.01,3,10,&AAMutSelSBDPOmegaModel::AAHyperLogProb,&AAMutSelSBDPOmegaModel::NoUpdate,this);
            ScalingMove(aaconc,1.0,10,&AAMutSelSBDPOmegaModel::AAHyperLogProb,&AAMutSelSBDPOmegaModel::NoUpdate,this);
            ScalingMove(aaconc,0.3,10,&AAMutSelSBDPOmegaModel::AAHyperLogProb,&AAMutSelSBDPOmegaModel::NoUpdate,this);
            /*
            ScalingMove(aainvconc,1.0,10,&AAMutSelSBDPOmegaModel::AAHyperLogProb,&AAMutSelSBDPOmegaModel::NoUpdate,this);
            ScalingMove(aainvconc,0.3,10,&AAMutSelSBDPOmegaModel::AAHyperLogProb,&AAMutSelSBDPOmegaModel::NoUpdate,this);
            */
        }
        componentaafitnessarray->SetCenter(aacenter);
        componentaafitnessarray->SetConcentration(aaconc);
        // componentaafitnessarray->SetConcentration(1.0/aainvconc);
        componentaafitnessarray->PriorResample(*occupancy);
        componentcodonmatrixarray->UpdateCodonMatrices(*occupancy);
    }

    void MoveAAMixture()    {

        for (int rep=0; rep<3; rep++)  {
            MoveAAProfiles();
            ResampleAlloc();
            LabelSwitchingMove();
            ResampleWeights();
            UpdateCodonMatrices();
            CollectComponentPathSuffStat();
        }
    }

    double MoveAAProfiles() {
        MoveAAProfiles(1.0,1,3);
        MoveAAProfiles(1.0,3,3);
        MoveAAProfiles(0.3,3,3);
        MoveAAProfiles(0.1,5,3);
        componentaafitnessarray->PriorResample(*occupancy);
        componentcodonmatrixarray->UpdateCodonMatrices(*occupancy);
        return 1.0;
    }

    void ResampleAlloc()    {
        vector<double> postprob(Ncat,0);
        for (int i=0; i<Nsite; i++) {
            GetAllocPostProb(i,postprob);
            sitealloc->GibbsResample(i,postprob);
        }
        UpdateOccupancies();
    }

    void UpdateOccupancies()    {
        occupancy->Clear();
        occupancy->AddSuffStat(*sitealloc);
    }

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

    void LabelSwitchingMove()   {
        MoveOccupiedCompAlloc(5);
        MoveAdjacentCompAlloc(5);
    }

    double MoveOccupiedCompAlloc(int k0)	{

        const vector<double>& w = weight->GetArray();

        int nrep = (int) (k0 * kappa);
        ResampleWeights();
        double total = 0.0;
        int Nocc = GetNcluster();
        if (Nocc != 1)	{
            for (int i=0; i<nrep; i++)	{
                int occupiedComponentIndices[Nocc];
                int j=0;
                for (int k=0; k<Ncat; k++)	{
                    if ((*occupancy)[k] != 0)	{
                        occupiedComponentIndices[j] = k;
                        j++;
                    }
                }
                if (j != Nocc)	{
                    cerr << "error in MoveOccupiedCompAlloc.\n";
                    exit(1);
                }
                int indices[2];
                Random::DrawFromUrn(indices,2,Nocc);
                int cat1 = occupiedComponentIndices[indices[0]];
                int cat2 = occupiedComponentIndices[indices[1]];
                double logMetropolis = ((*occupancy)[cat2] - (*occupancy)[cat1]) * log(w[cat1] / w[cat2]);
                int accepted = (log(Random::Uniform()) < logMetropolis);
                if (accepted)	{
                    total += 1.0;
                    componentaafitnessarray->Swap(cat1,cat2);
                    sitealloc->SwapComponents(cat1,cat2);
                    occupancy->Swap(cat1,cat2);
                }
            }
            return total /= nrep;
        }
        return 0;
    }

    double MoveAdjacentCompAlloc(int k0)	{

        ResampleWeights();
        int nrep = (int) (k0 * kappa);
        
        double total = 0;

        const vector<double>& V = weight->GetBetaVariates();

        for (int i=0; i<nrep; i++)	{
            int cat1 = (int)(Random::Uniform() * (Ncat-2));  
            int cat2 = cat1 + 1;
            double logMetropolis = ((*occupancy)[cat1] * log(1 - V[cat2])) - ((*occupancy)[cat2] * log(1-V[cat1]));
            int accepted = (log(Random::Uniform()) < logMetropolis);
            if (accepted)	{
                total += 1.0;
                componentaafitnessarray->Swap(cat1,cat2);
                sitealloc->SwapComponents(cat1,cat2);
                weight->SwapComponents(cat1,cat2);
                occupancy->Swap(cat1,cat2);
            }
        }

        return total /= nrep;
    }

    void ResampleWeights()  {
        weight->GibbsResample(*occupancy);
    }

    void MoveKappa()    {
        ScalingMove(kappa,1.0,10,&AAMutSelSBDPOmegaModel::StickBreakingHyperLogProb,&AAMutSelSBDPOmegaModel::NoUpdate,this);
        ScalingMove(kappa,0.3,10,&AAMutSelSBDPOmegaModel::StickBreakingHyperLogProb,&AAMutSelSBDPOmegaModel::NoUpdate,this);
        weight->SetKappa(kappa);
    }

    int GetNcluster() const {

        int n = 0;
        for (int i=0; i<Ncat; i++)  {
            if (occupancy->GetVal(i))    {
                n++;
            }
        }
        return n;
    }

	double MoveAAProfiles(double tuning, int n, int nrep)	{
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

    //-------------------
    // Traces and Monitors
    // ------------------

    double GetMeanAAEntropy() const {
        return componentaafitnessarray->GetMeanEntropy();
    }

    double GetNucRREntropy() const  {
        return Random::GetEntropy(nucrelrate);
    }

    double GetNucStatEntropy() const    {
        return Random::GetEntropy(nucrelrate);
    }

	void TraceHeader(std::ostream& os) const {
		os << "#logprior\tlnL\tlength\t";
		os << "omega\t";
        os << "ncluster\t";
        os << "kappa\t";
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
        os << GetNcluster() << '\t';
        os << kappa << '\t';
        os << GetMeanAAEntropy() << '\t';
        os << aaconc << '\t';
        os << Random::GetEntropy(aacenter) << '\t';
		os << Random::GetEntropy(nucstat) << '\t';
		os << Random::GetEntropy(nucrelrate) << '\n';
	}

	void Monitor(ostream& os) const {
        os << totchrono.GetTime() << '\t' << aachrono.GetTime() << '\n';
        os << "prop time in aa moves: " << aachrono.GetTime() / totchrono.GetTime() << '\n';
    }

	void FromStream(istream& is) {}
	void ToStream(ostream& os) const {}

};


