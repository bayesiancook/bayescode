
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
#include "StickBreakingProcess.hpp"
#include "MultinomialAllocationVector.hpp"
#include "Chrono.hpp"
#include "Permutation.hpp"

class AAMutSelDSBDPOmegaModel : public ProbModel {

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
    ConstMixtureArray<vector<double> >* componentcenterarray;
    ConstMixtureArray<double>* componentconcentrationarray;

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

	// this one is used by PhyloProcess: has to be a ConstArray<SubMatrix>
	ConstMixtureArray<SubMatrix>* sitesubmatrixarray;

	PhyloProcess* phyloprocess;

	PathSuffStatArray* sitepathsuffstatarray;
	PathSuffStatArray* componentpathsuffstatarray;

    int fixbl;
    int fixbasemix;
    int fixomega;

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

	AAMutSelDSBDPOmegaModel(string datafile, string treefile, int inNcat, int inbaseNcat)   {

        fixbl = 0;
        fixbasemix = 0;

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
        if (baseNcat == -1) {
            baseNcat = 100;
        }

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

        acca1=acca2=acca3=acca4=0;
        tota1=tota2=tota3=tota4=0;
        accb1=accb2=accb3=accb4=0;
        totb1=totb2=totb3=totb4=0;

		// Allocate();
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
        componentcenterarray = new ConstMixtureArray<vector<double> >(basecenterarray,componentalloc);
        componentconcentrationarray = new ConstMixtureArray<double>(baseconcentrationarray,componentalloc);

        componentaafitnessarray = new MultiDirichlet(componentcenterarray,componentconcentrationarray);

        // mixture of aa fitness profiles
        kappa = 1.0;
        weight = new StickBreakingProcess(Ncat,kappa);
        occupancy = new OccupancySuffStat(Ncat);

        sitealloc = new MultinomialAllocationVector(Nsite,weight->GetArray());

        omegahypermean = 1.0;
        omegahyperinvshape = 1.0;
		omega = 1.0;

        componentcodonmatrixarray = new AAMutSelOmegaCodonSubMatrixArray(GetCodonStateSpace(), nucmatrix, componentaafitnessarray, omega);
        sitesubmatrixarray = new ConstMixtureArray<SubMatrix>(componentcodonmatrixarray,sitealloc);

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

    int GetNsite() const    {
        return Nsite;
    }

    double GetOmega() const {
        return omega;
    }

    const PoissonSuffStatBranchArray* GetLengthSuffStatArray() const {
        return lengthsuffstatarray;
    }

    const DirichletSuffStatArray* GetBaseSuffStatArray() const   {
        return basesuffstatarray;
    }

    const OccupancySuffStat* GetBaseOccupancies() const {
        return baseoccupancy;
    }

    //-------------------
    // Setting and updating
    // ------------------

    void SetFixBL(int infixbl)    {
        fixbl = infixbl;
    }

    void SetFixOmega(int inom)  {
        fixomega = inom;
    }

    void SetFixBaseMix(int inmix)    {
        fixbasemix = inmix;
    }

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

    void SetBaseMixture(const ConstArray<vector<double> >& inbasecenterarray, const ConstArray<double>& inbaseconcentrationarray, const ConstArray<double>& inbaseweight, const ConstArray<int>& inpermut) {

        basecenterarray->Copy(inbasecenterarray);
        baseconcentrationarray->Copy(inbaseconcentrationarray);
        baseweight->Copy(inbaseweight);
        componentalloc->Permute(inpermut);
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
        UpdateBaseOccupancies();
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
        total += BaseStickBreakingHyperLogPrior();
        total += BaseStickBreakingLogPrior();
        total += StickBreakingHyperLogPrior();
        total += StickBreakingLogPrior();
        total += BaseLogPrior();
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

    double NucRatesLogPrior() const {
        return 0;
    }

    double BaseStickBreakingHyperLogPrior() const   {
        return -basekappa/10;
    }

    double BaseStickBreakingLogPrior() const    {
        return baseweight->GetLogProb(basekappa);
    }

    double StickBreakingHyperLogPrior() const   {
        return -kappa/10;
    }

    double StickBreakingLogPrior() const    {
        return weight->GetLogProb(kappa);
    }

    double BaseLogPrior() const {
        double total = 0;
        total += basecenterarray->GetLogProb();
        total += baseconcentrationarray->GetLogProb();
        if (isinf(total))   {
            cerr << "in BaseLogPrior: inf\n";
            exit(1);
        }
        return total;
    }

    double BaseLogPrior(int k) const {
        double total = 0;
        total += basecenterarray->GetLogProb(k);
        total += baseconcentrationarray->GetLogProb(k);
        return total;
    }

    double AALogPrior() const {
        return componentaafitnessarray->GetLogProb();
    }

    double AALogPrior(int k) const {
        return componentaafitnessarray->GetLogProb(k);
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
		return lambdasuffstat.GetLogProb(1.0,lambda);
	}

    double BaseSuffStatLogProb(int k) const   {
        return basesuffstatarray->GetVal(k).GetLogProb(basecenterarray->GetVal(k),baseconcentrationarray->GetVal(k));
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

    // for moving aa hyper params (center and concentration)
    // for component k of the mixture
    double BaseLogProb(int k) const   {
        return BaseLogPrior(k) + BaseSuffStatLogProb(k);
    }

    // for moving basekappa
    double BaseStickBreakingHyperLogProb() const {
        return BaseStickBreakingHyperLogPrior() + BaseStickBreakingLogPrior();
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
		phyloprocess->AddPathSuffStat(*sitepathsuffstatarray);
	}

    // per component of the mixture
	void CollectComponentPathSuffStat()	{
        componentpathsuffstatarray->Clear();
        sitepathsuffstatarray->AddToComponents(*componentpathsuffstatarray,*sitealloc);
    }

    void CollectLengthSuffStat()    {
		lengthsuffstatarray->Clear();
		phyloprocess->AddLengthSuffStat(*lengthsuffstatarray);
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

            if (! fixomega) {
                MoveOmega();
            }

            aachrono.Start();
            MoveAAMixture(3);
            aachrono.Stop();

            basechrono.Start();
            if (baseNcat > 1)   {
                ResampleBaseAlloc();
            }
            if (! fixbasemix) {
                MoveBaseMixture(3);
            }
            basechrono.Stop();

            totchrono.Stop();
		}
	}

	void ResampleBranchLengths()	{
        CollectLengthSuffStat();
		branchlength->GibbsResample(*lengthsuffstatarray);
	}

	void MoveBranchLengthsHyperParameter()	{
		lambdasuffstat.Clear();
		branchlength->AddSuffStat(lambdasuffstat);
        ScalingMove(lambda,1.0,10,&AAMutSelDSBDPOmegaModel::BranchLengthsHyperLogProb,&AAMutSelDSBDPOmegaModel::NoUpdate,this);
        ScalingMove(lambda,0.3,10,&AAMutSelDSBDPOmegaModel::BranchLengthsHyperLogProb,&AAMutSelDSBDPOmegaModel::NoUpdate,this);
		branchlength->SetScale(lambda);
	}

	void MoveOmega()	{

		omegasuffstat.Clear();
		omegasuffstat.AddSuffStat(*componentcodonmatrixarray,*componentpathsuffstatarray);
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
		omega = Random::GammaSample(alpha + omegasuffstat.GetCount(), beta + omegasuffstat.GetBeta());
		UpdateCodonMatrices();
	}

	void MoveNucRates()	{

        ProfileMove(nucrelrate,0.1,1,3,&AAMutSelDSBDPOmegaModel::NucRatesLogProb,&AAMutSelDSBDPOmegaModel::UpdateMatrices,this);
        ProfileMove(nucrelrate,0.03,3,3,&AAMutSelDSBDPOmegaModel::NucRatesLogProb,&AAMutSelDSBDPOmegaModel::UpdateMatrices,this);
        ProfileMove(nucrelrate,0.01,3,3,&AAMutSelDSBDPOmegaModel::NucRatesLogProb,&AAMutSelDSBDPOmegaModel::UpdateMatrices,this);

        ProfileMove(nucstat,0.1,1,3,&AAMutSelDSBDPOmegaModel::NucRatesLogProb,&AAMutSelDSBDPOmegaModel::UpdateMatrices,this);
        ProfileMove(nucstat,0.01,1,3,&AAMutSelDSBDPOmegaModel::NucRatesLogProb,&AAMutSelDSBDPOmegaModel::UpdateMatrices,this);
	}

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

    void ResampleEmptyComponents()  {
        componentaafitnessarray->PriorResample(*occupancy);
        componentcodonmatrixarray->UpdateCodonMatrices(*occupancy);
    }

    void MoveAAProfiles()   {
        CompMoveAAProfiles(3);
        MulMoveAAProfiles(3);
    }

    double CompMoveAAProfiles(int nrep) {
        accb1 += MoveAA(1.0,1,nrep);
        accb2 += MoveAA(1.0,3,nrep);
        accb3 += MoveAA(0.3,3,nrep);
        accb4 += MoveAA(0.1,3,nrep);
        totb1++;
        totb2++;
        totb3++;
        totb4++;
        return 1.0;
    }

    double MulMoveAAProfiles(int nrep) {
        acca1 += MoveAAGamma(3.0,nrep);
        acca2 += MoveAAGamma(1.0,nrep);
        acca3 += MoveAAGamma(0.3,nrep);
        acca4 += MoveAAGamma(0.1,nrep);
        tota1++;
        tota2++;
        tota3++;
        tota4++;
        return 1.0;
    }

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

    double GammaAALogPrior(const vector<double>& x, const vector<double>& aacenter, double aaconc) {
        double total = 0;
        for (int l=0; l<Naa; l++)   {
            total += (aaconc*aacenter[l] -1)*log(x[l]) - x[l] - Random::logGamma(aaconc*aacenter[l]);
        }
        return total;
    }

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
        sitealloc->AddSuffStat(*occupancy);
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
        Permutation permut(Ncat);
        weight->LabelSwitchingMove(5,*occupancy,permut);
        sitealloc->Permute(permut);
        componentaafitnessarray->Permute(permut);
    }

    void ResampleWeights()  {
        weight->GibbsResample(*occupancy);
    }

    void MoveKappa()    {
        ScalingMove(kappa,1.0,10,&AAMutSelDSBDPOmegaModel::StickBreakingHyperLogProb,&AAMutSelDSBDPOmegaModel::NoUpdate,this);
        ScalingMove(kappa,0.3,10,&AAMutSelDSBDPOmegaModel::StickBreakingHyperLogProb,&AAMutSelDSBDPOmegaModel::NoUpdate,this);
        weight->SetKappa(kappa);
    }

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

    void CollectBaseSuffStat()   {
        basesuffstatarray->Clear();
        componentaafitnessarray->AddSuffStat(*basesuffstatarray,*componentalloc);
    }

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

    void ResampleBaseEmptyComponents()  {
        basecenterarray->PriorResample(*baseoccupancy);
        baseconcentrationarray->PriorResample(*baseoccupancy);
    }

    void ResampleBaseAlloc()    {
        vector<double> postprob(baseNcat,0);
        for (int i=0; i<Ncat; i++) {
            GetBaseAllocPostProb(i,postprob);
            componentalloc->GibbsResample(i,postprob);
        }
        UpdateBaseOccupancies();
    }

    void UpdateBaseOccupancies()    {
        baseoccupancy->Clear();
        componentalloc->AddSuffStat(*baseoccupancy);
    }

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

    void BaseLabelSwitchingMove()   {
        Permutation permut(baseNcat);
        baseweight->LabelSwitchingMove(5,*baseoccupancy,permut);
        componentalloc->Permute(permut);
        basecenterarray->Permute(permut);
        baseconcentrationarray->Permute(permut);
        basesuffstatarray->Permute(permut);
    }

    void ResampleBaseWeights()  {
        baseweight->GibbsResample(*baseoccupancy);
    }

    void MoveBaseKappa()    {
        ScalingMove(basekappa,1.0,10,&AAMutSelDSBDPOmegaModel::BaseStickBreakingHyperLogProb,&AAMutSelDSBDPOmegaModel::NoUpdate,this);
        ScalingMove(basekappa,0.3,10,&AAMutSelDSBDPOmegaModel::BaseStickBreakingHyperLogProb,&AAMutSelDSBDPOmegaModel::NoUpdate,this);
        baseweight->SetKappa(basekappa);
    }

    //-------------------
    // Traces and Monitors
    // ------------------

    int GetNcluster() const {

        int n = 0;
        for (int i=0; i<Ncat; i++)  {
            if (occupancy->GetVal(i))    {
                n++;
            }
        }
        return n;
    }

    int GetBaseNcluster() const {

        int n = 0;
        for (int i=0; i<baseNcat; i++)  {
            if (baseoccupancy->GetVal(i))    {
                n++;
            }
        }
        return n;
    }

    double GetMeanAAEntropy() const {
        return componentaafitnessarray->GetMeanEntropy();
    }

    double GetMeanComponentAAConcentration() const {

        double tot = 0;
        for (int i=0; i<baseNcat; i++)  {
            tot += baseoccupancy->GetVal(i) * baseconcentrationarray->GetVal(i);
        }
        return tot / Ncat;
    }

    double GetMeanComponentAAEntropy() const {

        double tot = 0;
        for (int i=0; i<baseNcat; i++)  {
            tot += baseoccupancy->GetVal(i) * Random::GetEntropy(basecenterarray->GetVal(i));
        }
        return tot / Ncat;
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
        os << "basencluster\t";
        os << "basekappa\t";
        os << "aaent\t";
		os << "meanaaconc\t";
        os << "aacenterent\t";
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
        os << GetBaseNcluster() << '\t';
        os << basekappa << '\t';
        os << GetMeanAAEntropy() << '\t';
		os << GetMeanComponentAAConcentration() << '\t';
        os << GetMeanComponentAAEntropy() << '\t';
		os << Random::GetEntropy(nucstat) << '\t';
		os << Random::GetEntropy(nucrelrate) << '\n';
	}

	void Monitor(ostream& os) const {
        os << totchrono.GetTime() << '\t' << aachrono.GetTime() << '\t' << basechrono.GetTime() << '\n';
        os << "prop time in aa moves  : " << aachrono.GetTime() / totchrono.GetTime() << '\n';
        os << "prop time in base moves: " << basechrono.GetTime() / totchrono.GetTime() << '\n';
    }

	void FromStream(istream& is) {}
	void ToStream(ostream& os) const {}

};


